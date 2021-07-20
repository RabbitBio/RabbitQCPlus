//
// Created by ylf9811 on 2021/7/11.
//

#ifndef RABBITQCPLUS_ADAPTER_H
#define RABBITQCPLUS_ADAPTER_H

#include "Reference.h"
#include "repoter.h"



struct OverlapRes {
    bool overlaped;
    int offset;
    int overlap_len;
    int diff_num;
};

class Adapter {
public:
    Adapter();

    static OverlapRes AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require);

    static bool TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len);

    static bool TrimAdapter(neoReference &ref, std::string &adapter_seq, bool isR2 = false);

    static int CorrectData(neoReference &r1, neoReference &r2, OverlapRes &overlap_res);


};


#endif //RABBITQCPLUS_ADAPTER_H
