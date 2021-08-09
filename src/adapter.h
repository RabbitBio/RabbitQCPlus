//
// Created by ylf9811 on 2021/7/11.
//

#ifndef RABBITQCPLUS_ADAPTER_H
#define RABBITQCPLUS_ADAPTER_H

#include <iostream>
#include <vector>
#include <cstring>
#include <map>
#include <algorithm>
#include "Reference.h"
#include "repoter.h"
#include "Formater.h"
#include "nucleotidetree.h"

struct OverlapRes {
    bool overlapped;
    int offset;
    int overlap_len;
    int diff_num;
};

class Adapter {
public:
    Adapter();

    static std::string int2seq(unsigned int val, int seqlen);

    static int seq2int(const char *seq, int pos, int keylen, int lastVal);
//    static int seq2int(std::string &seq, int pos, int keylen, int lastVal);

    static std::string matchKnownAdapter(std::string seq);

    static std::string AutoDetect(std::string file_name, int trim_tail);

    static std::string
    getAdapterWithSeed(int seed, std::vector<neoReference> loadedReads, long records, int keylen, int trim_tail);

    static OverlapRes AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require);

    static bool TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len);

    static bool TrimAdapter(neoReference &ref, std::string &adapter_seq, bool isR2 = false);

    static int CorrectData(neoReference &r1, neoReference &r2, OverlapRes &overlap_res);


};


#endif //RABBITQCPLUS_ADAPTER_H
