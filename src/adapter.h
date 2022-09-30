//
// Created by ylf9811 on 2021/7/11.
//

#ifndef RERABBITQC_ADAPTER_H
#define RERABBITQC_ADAPTER_H

#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <omp.h>
#include <set>
#include <unordered_map>
#include <vector>

#include "Formater.h"
#include "Reference.h"
#include "nucleotidetree.h"
#include "repoter.h"

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

    static std::vector<std::string> LoadAdaptersFromFasta(std::string file_name);

    static std::string
    getAdapterWithSeed(int seed, std::vector<neoReference> loadedReads, long records, int keylen, int trim_tail);

    static OverlapRes AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require);

    static int TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len,
                           std::unordered_map<std::string, int> &mp1, std::unordered_map<std::string, int> &mp2,
                           int adapter_len_lim);

    static int TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len);

    static int TrimAdapter(neoReference &ref, std::string &adapter_seq, std::unordered_map<std::string, int> &mp,
                           int adapter_len_lim, bool isR2, int matchReq = 4);

    static int TrimAdapter(neoReference &ref, std::string &adapter_seq, bool isR2 = false, int matchReq = 4);

    static int TrimAdapters(neoReference &ref, std::vector<std::string> &adapter_seqs, std::unordered_map<std::string, int> &mp,
                           int adapter_len_lim, bool isR2);

    static int TrimAdapters(neoReference &ref, std::vector<std::string> &adapter_seqs, bool isR2 = false);


    static int CorrectData(neoReference &r1, neoReference &r2, OverlapRes &overlap_res, bool isPhred64);
};


#endif//RERABBITQC_ADAPTER_H
