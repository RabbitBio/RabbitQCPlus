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
    out_queue_ = NULL;

    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != std::string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != std::string::npos;
    if (cmd_info1->write_data_) {
        out_queue_ = new moodycamel::ConcurrentQueue<std::pair<char *, int>>(out_block_nums + 1000);
        if (out_is_zip_) {
            if (cmd_info1->use_pigz_) {
                string out_name1 = cmd_info1->out_file_name1_;
                out_name1 = out_name1.substr(0, out_name1.find(".gz"));
                out_stream_ = std::fstream(out_name1, std::ios::out | std::ios::binary);
                out_stream_.close();
//                out_stream_ = std::fstream(out_name1 + "_tmp", std::ios::out | std::ios::binary);

                printf("now use pigz to compress output data\n");
            } else {
                printf("open gzip stream %s\n", cmd_info1->out_file_name1_.c_str());
                zip_out_stream = gzopen(cmd_info1->out_file_name1_.c_str(), "w");
                gzsetparams(zip_out_stream, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream, 1024 * 1024);
            }
        } else {
            printf("open stream %s\n", cmd_info1->out_file_name1_.c_str());
            out_stream_ = std::fstream(cmd_info1->out_file_name1_, std::ios::out | std::ios::binary);
        }
    }
    duplicate_ = NULL;
    if (cmd_info1->state_duplicate_) {
        duplicate_ = new Duplicate(cmd_info1);
    }
    umier_ = NULL;
    if (cmd_info1->add_umi_) {
        umier_ = new Umier(cmd_info1);
    }
    if (cmd_info1->use_pugz_) {
        pugzQueue = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 20);
    }
    if (cmd_info1->use_pigz_) {
        pigzQueue = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 20);
        pigzLast.first = new char[1 << 24];
        pigzLast.second = 0;
    }
    pugzDone = 0;
    producerDone = 0;
    writerDone = 0;
}

SeQc::~SeQc() {
    delete filter_;
    if (cmd_info_->write_data_) {
        delete out_queue_;
    }

    if (cmd_info_->state_duplicate_) {
        delete duplicate_;
    }
    if (cmd_info_->add_umi_) {
        delete umier_;
    }
}


/**
 * @brief get fastq data chunk from fastq_data_pool and put it into data queue
 * @param file : fastq file name, which is also the input file
 * @param fastq_data_pool : fastq data pool
 * @param dq : a data queue
 */

void SeQc::ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    double t0 = GetTime();

    rabbit::fq::FastqFileReader *fqFileReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, fastq_data_pool, "", in_is_zip_);
    int64_t n_chunks = 0;

    if (cmd_info_->use_pugz_) {
        pair<char *, int> last_info;
        last_info.first = new char[1 << 20];
        last_info.second = 0;
        while (true) {
            rabbit::fq::FastqDataChunk *fqdatachunk;
            fqdatachunk = fqFileReader->readNextChunk(pugzQueue, &pugzDone, last_info);
            if (fqdatachunk == NULL) break;
            n_chunks++;
            //std::cout << "readed chunk: " << n_chunks << std::endl;
            dq.Push(n_chunks, fqdatachunk);
        }
    } else {
        while (true) {
            rabbit::fq::FastqDataChunk *fqdatachunk;
            fqdatachunk = fqFileReader->readNextChunk();
            if (fqdatachunk == NULL) break;
            n_chunks++;
            //std::cout << "readed chunk: " << n_chunks << std::endl;
            dq.Push(n_chunks, fqdatachunk);
        }
    }


    dq.SetCompleted();
    delete fqFileReader;
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
    printf("producer cost %.3f\n", GetTime() - t0);
    producerDone = 1;

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
    rabbit::fq::FastqDataChunk *fqdatachunk;

    if (cmd_info_->is_TGS_) {
        while (dq.Pop(id, fqdatachunk)) {
            std::vector<neoReference> data;
            rabbit::fq::chunkFormat(fqdatachunk, data, true);
            for (auto item:data) {
                thread_info->TGS_state_->tgsStatRead(item, cmd_info_->isPhred64_);
            }
            fastq_data_pool->Release(fqdatachunk);
        }
    } else {
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

                if (cmd_info_->add_umi_) {
                    umier_->ProcessSe(item);
                }
                bool trim_res = filter_->TrimSeq(item, cmd_info_->trim_front1_, cmd_info_->trim_tail1_);

                if (trim_res && cmd_info_->trim_polyg_) {
                    PolyX::trimPolyG(item, cmd_info_->trim_poly_len_);
                }

                if (trim_res && cmd_info_->trim_polyx_) {
                    PolyX::trimPolyX(item, cmd_info_->trim_poly_len_);
                }


                if (trim_res && cmd_info_->trim_adapter_ && cmd_info_->detect_adapter1_) {
                    Adapter::TrimAdapter(item, cmd_info_->adapter_seq1_, false);
                }
                int filter_res = filter_->ReadFiltering(item, trim_res, cmd_info_->isPhred64_);
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
                        //TODO delete name
                        Read2Chars(item, out_data, pos);
                    }
                    ASSERT(pos == out_len);
                    out_queue_->enqueue({out_data, out_len});
                }
            }

            fastq_data_pool->Release(fqdatachunk);
        }

    }
    done_thread_number_++;

}

/**
 * @brief a function to write data from out_data queue to file
 */
void SeQc::WriteSeFastqTask() {
    double t0 = GetTime();
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
            if (out_is_zip_) {
                if (cmd_info_->use_pigz_) {
                    while (pigzQueue->try_enqueue(now) == 0) {
                        printf("waiting to push a chunk to pigz queue\n");
                        usleep(100);
                    }
//                    out_stream_.write(now.first, now.second);
//                    printf("writer data to pigzQueue and disk\n");

                } else {
                    int written = gzwrite(zip_out_stream, now.first, now.second);
                    if (written != now.second) {
                        printf("GG");
                        exit(0);
                    }
                    delete[] now.first;
                }
            } else {
                out_stream_.write(now.first, now.second);
                delete[] now.first;
            }
        }
    }
    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {
//            out_stream_.close();

        } else {
            if (zip_out_stream) {
                gzflush(zip_out_stream, Z_FINISH);
                gzclose(zip_out_stream);
                zip_out_stream = NULL;
            }
        }
    } else {
        out_stream_.close();
    }
    printf("write cost %.5f\n", GetTime() - t0);
}


/**
 * @brief do pugz
 */

void SeQc::PugzTask() {
    printf("pugz start\n");
    auto t0 = GetTime();
    main_pugz(cmd_info_->in_file_name1_, cmd_info_->pugz_threads_, pugzQueue, &producerDone);
    printf("pugz cost %.5f\n", GetTime() - t0);
    pugzDone = 1;
}


void SeQc::PigzTask() {
    /*
 argc 9
argv ./pigz
argv -p
argv 16
argv -k
argv -4
argv -f
argv -b
argv -4096
argv p.fq
 */
    int cnt = 9;

    char **infos = new char *[9];
    infos[0] = "./pigz";
    infos[1] = "-p";
    int th_num = cmd_info_->pigz_threads_;
//    printf("th num is %d\n", th_num);
    string th_num_s = to_string(th_num);
//    printf("th num s is %s\n", th_num_s.c_str());
//    printf("th num s len is %d\n", th_num_s.length());

    infos[2] = new char[th_num_s.length() + 1];
    memcpy(infos[2], th_num_s.c_str(), th_num_s.length());
    infos[2][th_num_s.length()] = '\0';
    infos[3] = "-k";
    infos[4] = "-2";
    infos[5] = "-f";
    infos[6] = "-b";
    infos[7] = "4096";
    string out_name1 = cmd_info_->out_file_name1_;
    string out_file = out_name1.substr(0, out_name1.find(".gz"));
//    printf("th out_file is %s\n", out_file.c_str());
//    printf("th out_file len is %d\n", out_file.length());
    infos[8] = new char[out_file.length() + 1];
    memcpy(infos[8], out_file.c_str(), out_file.length());
    infos[8][out_file.length()] = '\0';

    main_pigz(cnt, infos, pigzQueue, &writerDone, pigzLast);

    printf("pigz done\n");
}


/**
 * @brief do QC for single-end data
 */

void SeQc::ProcessSeFastq() {

    thread *pugzer;

    if (cmd_info_->use_pugz_) {
        pugzer = new thread(bind(&::SeQc::PugzTask, this));
    }


    auto *fastqPool = new rabbit::fq::FastqDataPool(128, 1 << 22);
    //TODO replace this queue
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(128, 1);

    auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, false);
    }
    std::thread *write_thread;
    if (cmd_info_->write_data_) {
        write_thread = new std::thread(std::bind(&SeQc::WriteSeFastqTask, this));
    }
    thread *pigzer;
    if (cmd_info_->use_pigz_) {
        pigzer = new thread(bind(&SeQc::PigzTask, this));
    }

    //TODO bind ?
    std::thread producer(
            std::bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, std::ref(queue1)));
    auto **threads = new std::thread *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t] = new std::thread(
                std::bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info[t], fastqPool, std::ref(queue1)));
    }
    if (cmd_info_->use_pugz_) {
        pugzer->join();
    }

    producer.join();
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t]->join();
    }
    if (cmd_info_->write_data_) {
        write_thread->join();
        writerDone = 1;
    }
    if (cmd_info_->use_pigz_) {
        pigzer->join();
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
    printf("\nprint pre state info :\n");
    State::PrintStates(pre_state);
    printf("\nprint aft state info :\n");
    State::PrintStates(aft_state);

    if (cmd_info_->do_overrepresentation_) {
        auto hash_graph1 = pre_state->GetHashGraph();
        int hash_num1 = pre_state->GetHashNum();

        auto hash_graph2 = aft_state->GetHashGraph();
        int hash_num2 = aft_state->GetHashNum();


        ofstream ofs;
        ofs.open("ORP2.log", ifstream::out);
        for (int i = 0; i < hash_num1; i++) {
            ofs << hash_graph1[i].seq << " " << hash_graph1[i].cnt << "\n";
        }
        ofs.close();
        ofs.open("ORP3.log", ifstream::out);
        for (int i = 0; i < hash_num2; i++) {
            ofs << hash_graph2[i].seq << " " << hash_graph2[i].cnt << "\n";
        }
        ofs.close();
    }


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

    Repoter::ReportHtmlSe(pre_state, aft_state, cmd_info_->in_file_name1_, dupRate * 100.0);


    delete pre_state;
    delete aft_state;

    delete fastqPool;


    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        delete threads[t];
        delete p_thread_info[t];
    }

    delete[] threads;
    delete[] p_thread_info;
    if (cmd_info_->write_data_) {
        delete write_thread;
    }

}

void SeQc::ProcessSeTGS() {

    thread *pugzer;
    if(cmd_info_->use_pugz_){
        pugzer=new thread(bind(&::SeQc::PugzTask, this));
    }


    auto *fastqPool = new rabbit::fq::FastqDataPool(128, 1 << 22);
    //TODO replace this queue
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(128, 1);

    auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, false);
    }
    //TODO bind ?
    std::thread producer(
            std::bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, std::ref(queue1)));
    auto **threads = new std::thread *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t] = new std::thread(
                std::bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info[t], fastqPool, std::ref(queue1)));
    }

    if (cmd_info_->use_pugz_){
        pugzer->join();
    }

    producer.join();
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t]->join();
    }

    printf("all thrad done\n");
    printf("now merge thread info\n");

    std::vector<TGSStats *> vec_state;

    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        vec_state.push_back(p_thread_info[t]->TGS_state_);
    }
    auto mer_state = TGSStats::merge(vec_state);

    printf("merge done\n");
    printf("\nprint TGS state info :\n");
    //mer_state->print();

//    report3(mer_state);
    mer_state->CalReadsLens();
    Repoter::ReportHtmlTGS(mer_state, cmd_info_->in_file_name1_);

    delete fastqPool;
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        delete threads[t];
        delete p_thread_info[t];
    }

    delete[] threads;
    delete[] p_thread_info;
}
