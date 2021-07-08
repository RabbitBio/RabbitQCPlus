//
// Created by ylf9811 on 2021/7/6.
//
#include <functional>
#include <cassert>
#include <cstring>
#include "seqc.h"


SeQc::SeQc(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;
    q20bases_ = 0;
    q30bases_ = 0;
    lines_ = 0;
    pass_lines_ = 0;
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    int out_block_nums = int(1.0 * cmd_info1->in_file_size1_ / cmd_info1->out_block_size_);
    printf("out_block_nums %d\n", out_block_nums);
    out_queue_ = new moodycamel::ConcurrentQueue<pair<char *, int>>(out_block_nums + 1000);
    if (cmd_info1->write_data_) {
        printf("open stream %s\n", cmd_info1->out_file_name1_.c_str());
        out_stream_ = std::fstream(cmd_info1->out_file_name1_, std::ios::out | std::ios::binary);
    }
}

SeQc::~SeQc() {}

void SeQc::StateInfo(ThreadInfo *thread_info, neoReference &ref) {
    int slen = ref.lseq;
    int qlen = ref.lqual;
    ASSERT(slen == qlen);
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);
    const char q20 = '5';
    const char q30 = '?';
    for (int i = 0; i < slen; i++) {
        char base = bases[i];
        char qual = quals[i];
        if (qual >= q30) {
            thread_info->q20bases_++;
            thread_info->q30bases_++;
        } else if (qual >= q20) {
            thread_info->q20bases_++;
        }
    }
}

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

void SeQc::PrintRead(neoReference &ref) {
    std::cout << std::string((char *) ref.base + ref.pname, ref.lname) << std::endl;
    std::cout << std::string((char *) ref.base + ref.pseq, ref.lseq) << std::endl;
    std::cout << std::string((char *) ref.base + ref.pstrand, ref.lstrand) << std::endl;
    std::cout << std::string((char *) ref.base + ref.pqual, ref.lqual) << std::endl;
}

std::string SeQc::Read2String(neoReference &ref) {
    return std::string((char *) ref.base + ref.pname, ref.lname) + "\n" +
           std::string((char *) ref.base + ref.pseq, ref.lseq) + "\n" +
           std::string((char *) ref.base + ref.pstrand, ref.lstrand) + "\n" +
           std::string((char *) ref.base + ref.pqual, ref.lqual) + "\n";
}

void SeQc::ConsumerSeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastq_data_pool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    long line_sum = 0;
    long pass_line_sum = 0;
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataChunk *fqdatachunk;// = new rabbit::fq::FastqDataChunk;
    std::string out_data;
    while (dq.Pop(id, fqdatachunk)) {
        std::vector<neoReference> data;
        int res = rabbit::fq::chunkFormat(fqdatachunk, data, true);
        line_sum += res;
        for (auto item:data) {
            StateInfo(thread_info, item);
            int filter_res = filter_->ReadFiltering(item);
            if (filter_res == 0) {
                out_data += Read2String(item);
                pass_line_sum++;
                if (out_data.size() >= cmd_info_->out_block_size_) {
                    char *ldata = new char[out_data.size()];
                    memcpy(ldata, out_data.c_str(), out_data.size());
                    out_queue_->enqueue({ldata, out_data.size()});
                    out_data.clear();
                }
            }
        }

        fastq_data_pool->Release(fqdatachunk);
    }
    if (out_data.size() > 0) {
        char *ldata = new char[out_data.size()];
        memcpy(ldata, out_data.c_str(), out_data.size());
        out_queue_->enqueue({ldata, out_data.size()});
        out_data.clear();
    }
    thread_info->lines_ = line_sum;
    thread_info->pass_lines_ = pass_line_sum;
    done_thread_number_++;

}

void SeQc::WriteSeFastqTask() {
    while (true) {
        if (out_queue_->size_approx() == 0 && done_thread_number_ == cmd_info_->thread_number_) {
            break;
        }
        if (out_queue_->size_approx() == 0) {
            usleep(100);
        }
        pair<char *, int> now;
        while (out_queue_->size_approx()) {
            out_queue_->try_dequeue(now);
//            printf("write thread working, write %d\n", now.second);
            out_stream_.write(now.first, now.second);
            delete now.first;
        }
    }
    out_stream_.close();
}

void SeQc::ProcessSeFastq() {
    rabbit::fq::FastqDataPool *fastqPool = new rabbit::fq::FastqDataPool(32, 1 << 22);
    //TODO replace this queue
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(64, 1);

    ThreadInfo **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        p_thread_info[t] = new ThreadInfo();
    }
    std::thread *write_thread;
    if (cmd_info_->write_data_) {
        write_thread = new std::thread(std::bind(&SeQc::WriteSeFastqTask, this));
    }
    //TODO bind ?
    std::thread producer(
            std::bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, std::ref(queue1)));
    std::thread **threads = new std::thread *[cmd_info_->thread_number_];
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
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        q20bases_ += p_thread_info[t]->q20bases_;
        q30bases_ += p_thread_info[t]->q30bases_;
        lines_ += p_thread_info[t]->lines_;
        pass_lines_ += p_thread_info[t]->pass_lines_;
    }

    printf("merge done\n");
    printf("read lines %lld\n", lines_);
    printf("pass lines %lld\n", pass_lines_);
    printf("q20bases %lld\n", q20bases_);
    printf("q30bases %lld\n", q30bases_);


    delete fastqPool;
    delete[] threads;
    delete[] p_thread_info;
}
