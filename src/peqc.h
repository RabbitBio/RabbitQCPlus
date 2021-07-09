//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RABBITQCPLUS_PEQC_H
#define RABBITQCPLUS_PEQC_H

#include <functional>
#include "Globals.h"
#include "Formater.h"
#include "cmdinfo.h"


class PeQc {
public:
    PeQc(CmdInfo *cmd_info);

    PeQc();

    ~PeQc();

    void ProcessPeFastq();


private:
    void StateInfo(neoReference &ref);

    void PrintRead(neoReference &ref);

    void ProducerPeFastqTask(std::string file, std::string file2, rabbit::fq::FastqDataPool *fastqPool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq);

    void ConsumerPeFastqTask(rabbit::fq::FastqDataPool *fastqPool,
                             rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq);


private:
    CmdInfo *cmd_info;

};


#endif //RABBITQCPLUS_PEQC_H
