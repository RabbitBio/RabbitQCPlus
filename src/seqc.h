//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RERABBITQC_SEQC_H
#define RERABBITQC_SEQC_H

#include <atomic>
#include <cstring>
#include <fstream>
#include <functional>

#include "Formater.h"
#include "Globals.h"
#include "adapter.h"
#include "cmdinfo.h"
#include "concurrentqueue.h"
#include "duplicate.h"
#include "filter.h"
#include "pigz.h"
#include "polyx.h"
#include "pugz.h"
#include "state.h"
#include "threadinfo.h"
#include "umier.h"


class SeQc {
public:
    SeQc(CmdInfo *cmd_info, int my_rank = 0, int comm_size = 1);

    SeQc();

    ~SeQc();

    void ProcessSeFastq();

    void ProcessSeFastqOneThread();

    void ProcessSeTGS();

    void ProcessSeTGSOneThread();


private:
    //    void PrintRead(neoReference &ref);

    std::string Read2String(neoReference &ref);

    void Read2Chars(neoReference &ref, char *out_data, int &pos);

    void NGSTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool, ThreadInfo **thread_infos);

    void TGSTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool, ThreadInfo **thread_infos);

    void ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void ConsumerSeFastqTask(ThreadInfo **thread_infos, rabbit::fq::FastqDataPool *fastq_data_pool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void WriteSeFastqTask();

    void PugzTask();

    void PigzTask();

private:
    int my_rank, comm_size;
    CmdInfo *cmd_info_;

    Filter *filter_;
    moodycamel::ConcurrentQueue<std::pair<char *, int>> *
            out_queue_;
    //TODO replace concurrentqueue with char*[]
    std::atomic_int done_thread_number_;
    //std::ofstream out_stream_;
    FILE *out_stream_;
    Duplicate *duplicate_;
    Umier *umier_;
    long long now_pos_;

    gzFile zip_out_stream;
    bool in_is_zip_;
    bool out_is_zip_;

    moodycamel::ReaderWriterQueue<std::pair<char *, int>> *
            pugzQueue;
    std::atomic_int pugzDone;

    moodycamel::ReaderWriterQueue<std::pair<char *, int>> *
            pigzQueue;
    std::pair<char *, int> pigzLast;

    std::atomic_int producerDone;
    std::atomic_int writerDone;
    std::atomic_int queueNumNow;
    std::atomic_int pigzQueueNumNow;
    std::atomic_int queueSizeLim;
    std::atomic_int pigzQueueSizeLim;
    int64_t now_chunks;
    int64_t mx_chunks;
    std::atomic_int writeCommDone;
    std::atomic_int producerStop;
    std::mutex mylock;
};


#endif//RERABBITQC_SEQC_H
