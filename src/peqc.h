//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RABBITQCPLUS_PEQC_H
#define RABBITQCPLUS_PEQC_H

#include <atomic>
#include <fstream>
#include <functional>
#include <cstring>
#include <sys/time.h>

#include "Globals.h"
#include "Formater.h"
#include "cmdinfo.h"
#include "threadinfo.h"
#include "filter.h"
#include "concurrentqueue.h"
#include "state.h"
#include "adapter.h"
#include "duplicate.h"
#include "polyx.h"
#include "umier.h"


class PeQc {
public:
    PeQc(CmdInfo *cmd_info);

    PeQc();

    ~PeQc();

    void ProcessPeFastq();


private:

//    void PrintRead(neoReference &ref);

    std::string Read2String(neoReference &ref);

    void Read2Chars(neoReference &ref, char *out_data, int &pos);

//
//    void ProducerPeFastqTask(std::string file, std::string file2, rabbit::fq::FastqDataPool *fastqPool1,
//                             rabbit::fq::FastqDataPool *fastqPool2,
//                             rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq);
//
//    void ConsumerPeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool1,
//                             rabbit::fq::FastqDataPool *fastqPool2,
//                             rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq);

    void ProducerPeFastqTask(std::string file, std::string file2, rabbit::fq::FastqDataPool *fastqPool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq);

    void ConsumerPeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq);

    void WriteSeFastqTask1();

    void WriteSeFastqTask2();


private:
    CmdInfo *cmd_info_;
    Filter *filter_;
    moodycamel::ConcurrentQueue<std::pair<char *, int>> *out_queue1_;
    moodycamel::ConcurrentQueue<std::pair<char *, int>> *out_queue2_;
//TODO replace concurrentqueue with char*[]
    std::atomic_int done_thread_number_;
    std::fstream out_stream1_;
    std::fstream out_stream2_;
    Duplicate *duplicate_;
    Umier *umier_;
};


#endif //RABBITQCPLUS_PEQC_H
