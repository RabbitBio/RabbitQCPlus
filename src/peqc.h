//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RERABBITQC_PEQC_H
#define RERABBITQC_PEQC_H

#include <atomic>
#include <fstream>
#include <functional>
#include <cstring>
#include <sys/time.h>
#include <queue>

#include "Globals.h"
#include "Formater.h"
#include "cmdinfo.h"
#include "threadinfo.h"
#include "filter.h"
#include "state.h"
#include "adapter.h"
#include "duplicate.h"
#include "polyx.h"
#include "umier.h"
#include "pugz.h"
#include "pigz.h"

#define CIPair std::pair<char *, int>

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

    void ProducerPeInterFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                                  rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void ConsumerPeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq);

    void ConsumerPeInterFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
                                  rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void WriteSeFastqTask1();

    void WriteSeFastqTask2();

    void PugzTask1();

    void PugzTask2();

    void PigzTask1();

    void PigzTask2();

private:
    CmdInfo *cmd_info_;
    Filter *filter_;
    //moodycamel::ConcurrentQueue<std::pair<char *, int>> *out_queue1_;

    CIPair *out_queue1_;
    //moodycamel::ConcurrentQueue<std::pair<char *, int>> *out_queue2_;
    CIPair *out_queue2_;
//TODO replace concurrentqueue with char*[]
    std::atomic_int done_thread_number_;
    std::fstream out_stream1_;
    std::fstream out_stream2_;
    Duplicate *duplicate_;
    Umier *umier_;

    gzFile zip_out_stream1;
    gzFile zip_out_stream2;
    bool in_is_zip_;
    bool out_is_zip_;


    moodycamel::ReaderWriterQueue<std::pair < char * , int>> *
    pugzQueue1;
    moodycamel::ReaderWriterQueue<std::pair < char * , int>> *
    pugzQueue2;
    std::atomic_int pugzDone1;
    std::atomic_int pugzDone2;
    moodycamel::ReaderWriterQueue<std::pair < char * , int>> *
    pigzQueue1;
    pair<char *, int> pigzLast1;
    moodycamel::ReaderWriterQueue<std::pair < char * , int>> *
    pigzQueue2;
    pair<char *, int> pigzLast2;


    std::atomic_int producerDone;
    std::atomic_int writerDone1;
    std::atomic_int queue1P1;
    std::atomic_int queue1P2;
    std::atomic_int queueNumNow1;
    std::atomic_int queueSizeLim1;
    std::atomic_int pigzQueueNumNow1;
    std::atomic_int pigzQueueSizeLim1;

    std::atomic_int writerDone2;
    std::atomic_int queue2P1;
    std::atomic_int queue2P2;
    std::atomic_int queueNumNow2;
    std::atomic_int queueSizeLim2;
    std::atomic_int pigzQueueNumNow2;
    std::atomic_int pigzQueueSizeLim2;
    std::mutex mylock;


};


#endif //RERABBITQC_PEQC_H
