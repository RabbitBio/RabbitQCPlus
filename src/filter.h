//
// Created by ylf9811 on 2021/7/7.
//

#ifndef RERABBITQC_FILTER_H
#define RERABBITQC_FILTER_H

#include <cstring>

#include "Globals.h"
#include "Reference.h"
#include "cmdinfo.h"

class Filter {
public:
    Filter();

    Filter(CmdInfo *cmd_info);

    int ReadFiltering(neoReference &ref, bool trim_res, bool isPhred64);

    bool TrimSeq(neoReference &r, int front, int tail);


    //TODO
    void PrintResult();

private:
    CmdInfo *cmd_info;
    //TODO
    int64_t filter_res_cnt_[16];
};


#endif //RERABBITQC_FILTER_H
