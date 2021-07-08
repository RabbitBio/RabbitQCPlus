//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RABBITQCPLUS_SEQC_H
#define RABBITQCPLUS_SEQC_H

#include <atomic>
#include <fstream>
#include "FastxChunk.h"
#include "FastxStream.h"
#include "Formater.h"
#include "cmdinfo.h"
#include "threadinfo.h"
#include "filter.h"
#include "concurrentqueue.h"

class SeQc {
public:
    SeQc(CmdInfo *cmd_info);

    SeQc();

    ~SeQc();

    void ProcessSeFastq();


private:
    void StateInfo(ThreadInfo *thread_info, neoReference &ref);

    void PrintRead(neoReference &ref);

    std::string Read2String(neoReference &ref);

    void ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void ConsumerSeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastq_data_pool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void WriteSeFastqTask();

private:
    CmdInfo *cmd_info_;
    int64_t q20bases_;
    int64_t q30bases_;
    int64_t lines_;
    int64_t pass_lines_;
    Filter *filter_;
    moodycamel::ConcurrentQueue<pair<char *, int>> *out_queue_;
//replace concurrentqueue with char*[]
//    char **out_queue_;
    atomic_int queue_head_;
    atomic_int queue_tail_;
    atomic_int done_thread_number_;
    fstream out_stream_;


};


#endif //RABBITQCPLUS_SEQC_H
