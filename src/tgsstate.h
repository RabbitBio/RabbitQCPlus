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

        void tgsStatRead(neoReference &ref, bool isPhred64);

        void CalReadsLens();

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

        int64_t *tail_qual_sum;

        int *readsLens;

        int mMaxReadsLen;

        double mAvgReadsLen;

        int64_t mReadsNum;

        int64_t mBasesNum;

        int64_t *mBases51015Num;

        std::pair<double,int>* mTop5QualReads;

        std::pair<int,double>* mTop5LengReads;

        void updateTop5(int rlen,double avgQual);

    public:
        int64_t *const *GetHeadSeqPosCount() const;

        int64_t *const *GetTailSeqPosCount() const;

        const int64_t *GetHeadQualSum() const;

        const int64_t *GetTailQualSum() const;

        const int GetMaxReadsLen() const;

        const double GetAvgReadsLen() const;

        const int64_t GetReadsNum() const;

        const int64_t GetBasesNum() const;

        const int64_t* GetBases51015Num() const;

        const std::pair<double,int>* GetTop5QualReads() const;

        const std::pair<int,double>* GetTop5LengReads() const;

        int *GetReadsLens() const;

};


#endif //RERABBITQC_TGSSTATE_H
