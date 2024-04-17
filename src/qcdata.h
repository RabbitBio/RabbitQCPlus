

#ifndef QC_DATA_H
#define QC_DATA_H
#include "cmdinfo.h"
#include "threadinfo.h"
#include "Globals.h"
#include "duplicate.h"

struct dupInfo{
    uint32_t key;
    uint64_t kmer32;
    uint8_t gc;
};


struct qc_data_tgs {
    ThreadInfo **thread_info_;
    CmdInfo *cmd_info_;
    std::vector<neoReference> *data1_;
    int *cnt;
    int bit_len;
};

struct qc_data {
    ThreadInfo **thread_info_;
    CmdInfo *cmd_info_;
    Duplicate *duplicate_;
    //std::vector<dupInfo> *dups;
    std::vector<neoReference> *data1_[64];
    std::vector<neoReference> *data2_[64];
    std::vector<neoReference> *pass_data1_[64];
    std::vector<neoReference> *pass_data2_[64];
    int *cnt;
    int bit_len;
};

#endif // QC_DATA_H

