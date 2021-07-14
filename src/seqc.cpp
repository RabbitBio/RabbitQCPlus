//
// Created by ylf9811 on 2021/7/6.
//
#include "seqc.h"

/**
 * @brief Construct function
 * @param cmd_info1 : cmd information
 */
SeQc::SeQc(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    int out_block_nums = int(1.0 * cmd_info1->in_file_size1_ / cmd_info1->out_block_size_);
    printf("out_block_nums %d\n", out_block_nums);
    out_queue_ = new moodycamel::ConcurrentQueue<std::pair<char *, int>>
            (out_block_nums + 1000);
    if (cmd_info1->write_data_) {
        printf("open stream %s\n", cmd_info1->out_file_name1_.c_str());
        out_stream_ = std::fstream(cmd_info1->out_file_name1_, std::ios::out | std::ios::binary);
    }
    duplicate_ = NULL;
    if (cmd_info1->state_duplicate_) {
        duplicate_ = new Duplicate(cmd_info1);
    }
}

SeQc::~SeQc() {}


/**
 * @brief get fastq data chunk from fastq_data_pool and put it into data queue
 * @param file : fastq file name, which is also the input file
 * @param fastq_data_pool : fastq data pool
 * @param dq : a data queue
 */

void SeQc::ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::fq::FastqFileReader *fqFileReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastq_data_pool);
    rabbit::int64 n_chunks = 0;
    while (true) {
        rabbit::fq::FastqDataChunk *fqdatachunk;// = new rabbit::fq::FastqDataChunk;
        fqdatachunk = fqFileReader->readNextChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        //std::cout << "readed chunk: " << n_chunks << std::endl;
        dq.Push(n_chunks, fqdatachunk);
    }
    dq.SetCompleted();
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
}
//
//void SeQc::PrintRead(neoReference &ref) {
//    std::cout << std::string((char *) ref.base + ref.pname, ref.lname) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pseq, ref.lseq) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pstrand, ref.lstrand) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pqual, ref.lqual) << std::endl;
//}

std::string SeQc::Read2String(neoReference &ref) {
    return std::string((char *) ref.base + ref.pname, ref.lname) + "\n" +
           std::string((char *) ref.base + ref.pseq, ref.lseq) + "\n" +
           std::string((char *) ref.base + ref.pstrand, ref.lstrand) + "\n" +
           std::string((char *) ref.base + ref.pqual, ref.lqual) + "\n";
}

void SeQc::Read2Chars(neoReference &ref, char *out_data, int &pos) {
    memcpy(out_data + pos, ref.base + ref.pname, ref.lname);
    pos += ref.lname;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pseq, ref.lseq);
    pos += ref.lseq;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pstrand, ref.lstrand);
    pos += ref.lstrand;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pqual, ref.lqual);
    pos += ref.lqual;
    out_data[pos++] = '\n';
}

/**
 * @brief get fastq data chunks from the data queue and do QC for them
 * @param thread_info : thread information
 * @param fastq_data_pool :a fastq data pool, it will be used to release data chunk
 * @param dq : data queue
 */
void SeQc::ConsumerSeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastq_data_pool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataChunk *fqdatachunk;// = new rabbit::fq::FastqDataChunk;
    while (dq.Pop(id, fqdatachunk)) {
        std::vector<neoReference> data;
        std::vector<neoReference> pass_data;
        rabbit::fq::chunkFormat(fqdatachunk, data, true);
        int out_len = 0;
        for (auto item:data) {
            thread_info->pre_state1_->StateInfo(item);
            if (cmd_info_->state_duplicate_) {
                duplicate_->statRead(item);
            }

            bool trim_res = filter_->TrimSeq(item, cmd_info_->trim_front1_, cmd_info_->trim_tail1_);

//            printf("gogo %s\n",
//                   std::string(reinterpret_cast<const char *>(item.base + item.pname), item.lname).c_str());
            if (trim_res && cmd_info_->trim_polyg_) {
                PolyX::trimPolyG(item, cmd_info_->trim_poly_len_);
            }

            if (trim_res && cmd_info_->trim_polyx_) {
                PolyX::trimPolyX(item, cmd_info_->trim_poly_len_);
            }


            if (trim_res && cmd_info_->trim_adapter_ && cmd_info_->detect_adapter1_) {
//                printf("Adapter::TrimAdapter ...\n");
                Adapter::TrimAdapter(item, cmd_info_->adapter_seq1_, false);
            }
            int filter_res = filter_->ReadFiltering(item, trim_res);
            if (filter_res == 0) {
                thread_info->aft_state1_->StateInfo(item);
                if (cmd_info_->write_data_) {
                    pass_data.push_back(item);
                    out_len += item.lname + item.lseq + item.lstrand + item.lqual + 4;
                }
            }
        }

        if (cmd_info_->write_data_) {
            if (pass_data.size() > 0) {
                char *out_data = new char[out_len];
                int pos = 0;
                for (auto item:pass_data) {
                    Read2Chars(item, out_data, pos);
                }
                ASSERT(pos == out_len);
                out_queue_->enqueue({out_data, out_len});
            }
        }

        fastq_data_pool->Release(fqdatachunk);
    }
    done_thread_number_++;

}

/**
 * @brief a function to write data from out_data queue to file
 */
void SeQc::WriteSeFastqTask() {
    int cnt = 0;
    while (true) {
        if (out_queue_->size_approx() == 0 && done_thread_number_ == cmd_info_->thread_number_) {
            break;
        }
        if (out_queue_->size_approx() == 0) {
            usleep(100);
        }
        std::pair<char *, int> now;
        while (out_queue_->size_approx()) {
            out_queue_->try_dequeue(now);
//            printf("write thread working, write %d %d\n", now.second, cnt++);
            out_stream_.write(now.first, now.second);
            delete now.first;
        }
    }
    out_stream_.close();
}

/**
 * @brief do QC for single-end data
 */

void SeQc::ProcessSeFastq() {
    auto *fastqPool = new rabbit::fq::FastqDataPool(32, 1 << 22);
    //TODO replace this queue
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(64, 1);

    auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_);
    }
    std::thread *write_thread;
    if (cmd_info_->write_data_) {
        write_thread = new std::thread(std::bind(&SeQc::WriteSeFastqTask, this));
    }
    //TODO bind ?
    std::thread producer(
            std::bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, std::ref(queue1)));
    auto **threads = new std::thread *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t] = new std::thread(
                std::bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info[t], fastqPool, std::ref(queue1)));
    }
    producer.join();
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t]->join();
    }
    if (cmd_info_->write_data_) {
        write_thread->join();
    }

    printf("all thrad done\n");
    printf("now merge thread info\n");

    std::vector<State *> pre_vec_state;
    std::vector<State *> aft_vec_state;

    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        pre_vec_state.push_back(p_thread_info[t]->pre_state1_);
        aft_vec_state.push_back(p_thread_info[t]->aft_state1_);
    }
    auto pre_state = State::MergeStates(pre_vec_state);
    auto aft_state = State::MergeStates(aft_vec_state);

    printf("merge done\n");
    printf("print pre state info :\n");
    State::PrintStates(pre_state);
    printf("print aft state info :\n");
    State::PrintStates(aft_state);

    int *dupHist = NULL;
    double *dupMeanGC = NULL;
    double dupRate = 0.0;
    int histSize = 32;
    if (cmd_info_->state_duplicate_) {
        dupHist = new int[histSize];
        memset(dupHist, 0, sizeof(int) * histSize);
        dupMeanGC = new double[histSize];
        memset(dupMeanGC, 0, sizeof(double) * histSize);
        dupRate = duplicate_->statAll(dupHist, dupMeanGC, histSize);
        printf("Duplication rate (may be overestimated since this is SE data): %.5f %%\n", dupRate * 100.0);
        delete[] dupHist;
        delete[] dupMeanGC;

    }


    delete fastqPool;
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        delete threads[t];
    }
    delete[] threads;
    delete[] p_thread_info;

}
