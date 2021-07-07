//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RABBITQCPLUS_SEQC_H
#define RABBITQCPLUS_SEQC_H

#include "FastxChunk.h"
#include "FastxStream.h"
#include "Formater.h"
#include "cmdinfo.h"

class SeQc {
public:
    SeQc(CmdInfo *cmd_info);

    SeQc();

    ~SeQc();
    void ProcessSeFastq();


private:
    void StateInfo(neoReference &ref);

    void PrintRead(neoReference &ref);

    void ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastqPool,
                            rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);

    void ConsumerSeFastqTask(rabbit::fq::FastqDataPool *fastqPool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq);


private:
    CmdInfo *cmd_info;
    int64_t q20bases;
    int64_t q30bases;

};


#endif //RABBITQCPLUS_SEQC_H
