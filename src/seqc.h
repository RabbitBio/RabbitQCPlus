//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RABBITQCPLUS_SEQC_H
#define RABBITQCPLUS_SEQC_H

#include <atomic>
#include <fstream>
#include <functional>
#include <cstring>

#include "Globals.h"
#include "Formater.h"
#include "cmdinfo.h"
#include "threadinfo.h"
#include "filter.h"
#include "concurrentqueue.h"
#include "state.h"

class SeQc {
public:
    SeQc(CmdInfo *cmd_info);

    SeQc();

    ~SeQc();

    void ProcessSeFastq();


private:

//    void PrintRead(neoReference &ref);

    std::string Read2String(neoReference &ref);

    void Read2Chars(neoReference &ref, char *out_data, int &pos);

    void ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void ConsumerSeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastq_data_pool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void WriteSeFastqTask();

private:
    CmdInfo *cmd_info_;

    Filter *filter_;
    moodycamel::ConcurrentQueue<std::pair<char *, int>> *out_queue_;
//TODO replace concurrentqueue with char*[]
    std::atomic_int done_thread_number_;
    std::fstream out_stream_;


};


#endif //RABBITQCPLUS_SEQC_H
