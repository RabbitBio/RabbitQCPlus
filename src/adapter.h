//
// Created by ylf9811 on 2021/7/11.
//

#ifndef RERABBITQC_ADAPTER_H
#define RERABBITQC_ADAPTER_H

#include <iostream>
#include <vector>
#include <cstring>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <omp.h>
#include <set>

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

    static void PreOverAnalyze(std::string file_name, std::vector<std::string> &ho_seqs, int &eva_len);

    static std::string AutoDetect(std::string file_name, int trim_tail);
    
    static int EvalMaxLen(std::string file_name);

    static std::string
    getAdapterWithSeed(int seed, std::vector<neoReference> loadedReads, long records, int keylen, int trim_tail);

    static OverlapRes AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require);

    static int TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len, std::unordered_map<std::string,int> &mp1, std::unordered_map<std::string,int> &mp2, int adapter_len_lim);

    static int TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len);

    static int TrimAdapter(neoReference &ref, std::string &adapter_seq, std::unordered_map<std::string,int> &mp, int adapter_len_lim, bool isR2);

    static int TrimAdapter(neoReference &ref, std::string &adapter_seq, bool isR2 = false);

    static int CorrectData(neoReference &r1, neoReference &r2, OverlapRes &overlap_res, bool isPhred64);


};


#endif //RERABBITQC_ADAPTER_H
