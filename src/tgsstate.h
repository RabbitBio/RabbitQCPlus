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


    void reportJson(std::ofstream &ofs, std::string padding);

    // a port of HTML report
    void reportHtml(std::ofstream &ofs, std::string filteringType, std::string readName);

    void reportHtmlQuality(std::ofstream &ofs, std::string seqFileName, bool isTail, std::string xAxisName,
                           std::string yAxisName, double *statsData);

    void reportHtmlContents(std::ofstream &ofs, std::string seqFileName, bool isTail, std::string xAxisName,
                            std::string yAxisName, double **statsData);

    void reportHistogram(std::ofstream &ofs);

    bool isLongRead();

    static std::string list2string(double *list, int size);

    static std::string list2string(double *list, int size, long *coords);

    static std::string list2string(long *list, int size);

    static std::string list2stringReversedOrder(long *list, int size);

    int base2num(std::string base);

private:
    int mMinlen;
    int mHalfMinlen;
    std::vector<int> mLengths;
    std::vector<int> mTotalReadsLen;
    long *head_seq_pos_count[4];
    long *tail_seq_pos_count[4];
    long *head_qual_sum;
    long *tail_qual_sum;
};


#endif //RERABBITQC_TGSSTATE_H
