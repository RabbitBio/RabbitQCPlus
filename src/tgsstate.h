//
// Created by ylf9811 on 2021/8/11.
//

#ifndef RERABBITQC_TGSSTATE_H
#define RERABBITQC_TGSSTATE_H

#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "Reference.h"

class TGSStats {
public:
    TGSStats(int minLen);

    ~TGSStats();

    static TGSStats *merge(std::vector<TGSStats *> &list);

    void print();

    void tgsStatRead(neoReference &ref);

    // a port of HTML report
    void reportHtml(std::ofstream &ofs, std::string filteringType, std::string readName);

    void reportHtmlQuality(std::ofstream &ofs, std::string seqFileName, bool isTail, std::string xAxisName,
                           std::string yAxisName, double *statsData);

    void reportHtmlContents(std::ofstream &ofs, std::string seqFileName, bool isTail, std::string xAxisName,
                            std::string yAxisName, double **statsData);

    void reportHistogram(std::ofstream &ofs);

    bool isLongRead();

    static std::string list2string(double *list, int size);

    static std::string list2string(double *list, int size, int64_t *coords);

    static std::string list2string(int64_t *list, int size);

    static std::string list2stringReversedOrder(int64_t *list, int size);

    int base2num(std::string base);

private:
    int mMinlen;
    int mHalfMinlen;
    std::vector<int> mLengths;
    std::vector<int> mTotalReadsLen;
    int64_t *head_seq_pos_count[4];
    int64_t *tail_seq_pos_count[4];
    int64_t *head_qual_sum;
public:
    int64_t *const *GetHeadSeqPosCount() const;

    int64_t *const *GetTailSeqPosCount() const;

    int64_t *GetHeadQualSum() const;

    int64_t *GetTailQualSum() const;

private:
    int64_t *tail_qual_sum;
};


#endif //RERABBITQC_TGSSTATE_H
