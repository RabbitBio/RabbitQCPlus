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
#include "qcdata.h"

#include "globalMutex.h"


#define CIPair std::pair<char *, std::pair<int, long long>>
#define CRPair std::pair<int64_t, std::vector<rabbit::fq::FastqDataChunk *>>


class SeQc {
public:
    SeQc(CmdInfo *cmd_info, int my_rank = 0, int comm_size = 1);

    SeQc();

    ~SeQc();

    void ProcessSeFastq();

    void ProcessSeTGS();


private:
    //    void PrintRead(neoReference &ref);

    std::string Read2String(neoReference &ref);

    void Read2Chars(neoReference &ref, char *out_data, int &pos);

    void ProducerSeFastqTask64(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool);

    void ConsumerSeFastqTask64(ThreadInfo **thread_infos, rabbit::fq::FastqDataPool *fastq_data_pool);

    void ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void ConsumerSeFastqTask(ThreadInfo **thread_infos, rabbit::fq::FastqDataPool *fastq_data_pool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);


    void WriteSeFastqTask();

//    void ProcessNgsData(bool &proDone, std::vector <neoReference> &data, rabbit::fq::FastqDataChunk *fqdatachunk, qc_data *para, rabbit::fq::FastqDataPool *fastq_data_pool);
    void ProcessNgsData(bool &proDone, std::vector <neoReference> data[64], std::vector <rabbit::fq::FastqDataChunk *> fqdatachunks, qc_data *para, rabbit::fq::FastqDataPool *fastq_data_pool);


private:
    int my_rank, comm_size;
    CmdInfo *cmd_info_;

    Filter *filter_;
    CIPair *out_queue_;
    //QChunkItem *out_queue_;
    std::vector<rabbit::fq::FastqDataChunk *> *p_out_queue_;
    std::atomic_int done_thread_number_;
    FILE *out_stream_;
    MPI_File fh;
    MPI_Status status;
    Duplicate *duplicate_;
    //Duplicate *duplicate_[64];
    Umier *umier_;
    long long now_pos_;
    long long zip_now_pos_;

    gzFile zip_out_stream;
    std::ofstream off_idx;
    bool in_is_zip_;
    bool out_is_zip_;
    int start_line_;
    int end_line_;

    moodycamel::ReaderWriterQueue<std::pair<char *, int>> *
            pugzQueue;
    std::atomic_int pugzDone;

    moodycamel::ReaderWriterQueue<std::pair<char *, int>> *
            pigzQueue;
    std::pair<char *, int> pigzLast;

    std::atomic_int p_queueP1;
    std::atomic_int p_queueP2;
    std::atomic_int p_queueNumNow;
    std::atomic_int p_queueSizeLim;

    std::atomic_int producerDone;
    std::atomic_int writerDone;
    std::atomic_int queueP1;
    std::atomic_int queueP2;
    std::atomic_int queueNumNow;
    std::atomic_int pigzQueueNumNow;
    std::atomic_int queueSizeLim;
    std::atomic_int pigzQueueSizeLim;
    std::atomic_int now_chunks;
    int64_t mx_chunks;
    std::atomic_int consumerCommDone;
    std::atomic_int producerStop;
    std::mutex mylock;
    std::mutex p_mylock;

    int64_t *part_sizes;
};


#endif//RERABBITQC_SEQC_H
