//
// Created by ylf9811 on 2021/8/11.
//

#include "tgsstate.h"
using namespace std;

TGSStats::TGSStats(int minLen, bool isMerge) {
    if(isMerge) mTotalReadsLen = new int[64 << 20];
    else mTotalReadsLen = new int[1 << 20];
    countLenNum = 0;
    mMinlen = minLen;
    mHalfMinlen = mMinlen >> 1;
    int size_range = mMinlen >> 1;
    head_qual_sum = new int64_t[size_range];
    tail_qual_sum = new int64_t[size_range];
    int i;
    for (i = 0; i < 4; i++) {
        head_seq_pos_count[i] = new int64_t[size_range];
        tail_seq_pos_count[i] = new int64_t[size_range];
        memset(head_seq_pos_count[i], 0, sizeof(int64_t) * size_range);
        memset(tail_seq_pos_count[i], 0, sizeof(int64_t) * size_range);
    }
    //init
    memset(head_qual_sum, 0, size_range * sizeof(int64_t));
    memset(tail_qual_sum, 0, size_range * sizeof(int64_t));
    mReadsNum = 0;
    mBasesNum = 0;
    mBases51015Num = new int64_t[3];
    mBases51015Num[0] = 0;
    mBases51015Num[1] = 0;
    mBases51015Num[2] = 0;
    mTop5QualReads = new pair<double, int>[5];
    mTop5LengReads = new pair<int, double>[5];
    for (int i = 0; i < 5; i++) {
        mTop5QualReads[i] = {0, 0};
        mTop5LengReads[i] = {0, 0};
    }
}

TGSStats::~TGSStats() {
    delete[] head_qual_sum;
    delete[] tail_qual_sum;
    int i;
    for (i = 0; i < 4; i++) {
        delete[] head_seq_pos_count[i];
        delete[] tail_seq_pos_count[i];
    }
    delete mBases51015Num;
    delete mTop5QualReads;
    delete mTop5LengReads;
}

void TGSStats::updateTop5(int rlen, double avgQual) {
    int pos = 5;
    for (int i = 0; i < 5; i++) {
        if (mTop5QualReads[i].first < avgQual) {
            for (int j = 4; j > i; j--) {
                mTop5QualReads[j] = mTop5QualReads[j - 1];
            }
            mTop5QualReads[i] = {avgQual, rlen};
            break;
        }
    }
    for (int i = 0; i < 5; i++) {
        if (mTop5LengReads[i].first < rlen) {
            for (int j = 4; j > i; j--) {
                mTop5LengReads[j] = mTop5LengReads[j - 1];
            }
            mTop5LengReads[i] = {rlen, avgQual};
            break;
        }
    }
}


void TGSStats::print() {


    auto res3 = this->GetBases51015Num();
    printf("total bases %lld\n", this->GetBasesNum());
    printf("q5 bases %lld\n", res3[0]);
    printf("q10 bases %lld\n", res3[1]);
    printf("q15 bases %lld\n", res3[2]);
    printf("read number %lld\n", this->GetReadsNum());
    printf("average read length %.3f\n", this->GetAvgReadsLen());
    printf("\n");
    /*
    //cerr << "nothing here" << endl;
    int i;
    cout << "head A    T   C   G" << endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        cout << head_seq_pos_count[0][i] << " ";
        cout << head_seq_pos_count[1][i] << " ";
        cout << head_seq_pos_count[2][i] << " ";
        cout << head_seq_pos_count[3][i] << " ";
        cout << endl;
    }
    cout << "tail A    T   C   G" << endl;
    //cout << endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        cout << tail_seq_pos_count[0][i] << " ";
        cout << tail_seq_pos_count[1][i] << " ";
        cout << tail_seq_pos_count[2][i] << " ";
        cout << tail_seq_pos_count[3][i] << " ";
        cout << endl;
    }
    cout << "head and tail quality" << endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        cout << head_qual_sum[i] << " " << tail_qual_sum[i];
        cout << endl;
    }
    */
}

TGSStats *TGSStats::merge(vector<TGSStats *> &list) {
    const int minLen = list[0]->mMinlen;
    TGSStats *s = new TGSStats(minLen, true);
    for (TGSStats *ds: list) {
        //printf("ooo1\n");
        for(int i = 0; i < ds->countLenNum; i++){
          s->mTotalReadsLen[s->countLenNum++] = ds->mTotalReadsLen[i];
        }
        //printf("ooo2\n");
        for (int i = 0; i < (minLen >> 1); ++i) {
            //head_seq
            s->head_seq_pos_count[0][i] += ds->head_seq_pos_count[0][i];
            s->head_seq_pos_count[1][i] += ds->head_seq_pos_count[1][i];
            s->head_seq_pos_count[2][i] += ds->head_seq_pos_count[2][i];
            s->head_seq_pos_count[3][i] += ds->head_seq_pos_count[3][i];
            //tail_seq
            s->tail_seq_pos_count[0][i] += ds->tail_seq_pos_count[0][i];
            s->tail_seq_pos_count[1][i] += ds->tail_seq_pos_count[1][i];
            s->tail_seq_pos_count[2][i] += ds->tail_seq_pos_count[2][i];
            s->tail_seq_pos_count[3][i] += ds->tail_seq_pos_count[3][i];
        }
        //printf("ooo3\n");
        for (int i = 0; i < (minLen >> 1); ++i) {
            //head_qual
            s->head_qual_sum[i] += ds->head_qual_sum[i];
            //tail_qual
            s->tail_qual_sum[i] += ds->tail_qual_sum[i];
        }
        //printf("ooo4\n");
        s->mReadsNum += ds->mReadsNum;
        s->mBasesNum += ds->mBasesNum;
        s->mBases51015Num[0] += ds->mBases51015Num[0];
        s->mBases51015Num[1] += ds->mBases51015Num[1];
        s->mBases51015Num[2] += ds->mBases51015Num[2];
        //printf("ooo5\n");
        for (int i = 0; i < 5; i++) {
            s->updateTop5(ds->mTop5QualReads[i].second, ds->mTop5QualReads[i].first);
            s->updateTop5(ds->mTop5LengReads[i].first, ds->mTop5LengReads[i].second);
        }
        //printf("ooo6\n");
    }


    return s;
}

void TGSStats::CalReadsLens() {
    int Ppre = 10;
    mMaxReadsLen = 0;
    int64_t readLenSum = 0;
    for (int i = 0; i < countLenNum; i++) {
        mMaxReadsLen = max(mMaxReadsLen, mTotalReadsLen[i]);
        readLenSum += mTotalReadsLen[i];
    }
    mAvgReadsLen = 1.0 * readLenSum / countLenNum;
    mMaxReadsLen /= Ppre;
    mMaxReadsLen += Ppre;
    readsLens = new int[mMaxReadsLen];
    for (int i = 0; i < mMaxReadsLen; i++) readsLens[i] = 0;
    for (int i = 0; i < countLenNum; i++) {
        readsLens[(mTotalReadsLen[i] + Ppre - 1) / Ppre]++;
    }
}


//generate html data

bool TGSStats::isLongRead() {
    return mHalfMinlen > 300;
}

int TGSStats::base2num(string base) {
    if (base == "A")
        return 0;
    else if (base == "T")
        return 1;
    else if (base == "C")
        return 2;
    else if (base == "G")
        return 3;

    return -1;//fix warning (no ACGT)
}


string TGSStats::list2string(double *list, int size) {
    stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

string TGSStats::list2string(double *list, int size, int64_t *coords) {
    stringstream ss;
    for (int i = 0; i < size; i++) {
        // coords is 1,2,3,...
        int64_t start = 0;
        if (i > 0)
            start = coords[i - 1];
        int64_t end = coords[i];

        double total = 0.0;
        for (int k = start; k < end; k++)
            total += list[k];

        // get average
        if (end == start)
            ss << "0.0";
        else
            ss << total / (end - start);
        //ss << list[coords[i]-1];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

string TGSStats::list2string(int64_t *list, int size) {
    stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

string TGSStats::list2stringReversedOrder(int64_t *list, int size) {
    stringstream ss;
    for (int i = size - 1; i >= 0; --i) {
        ss << list[i];
        if (i > 0)
            ss << ",";
    }
    return ss.str();
}


int64_t *const *TGSStats::GetHeadSeqPosCount() const {
    return head_seq_pos_count;
}

int64_t *const *TGSStats::GetTailSeqPosCount() const {
    return tail_seq_pos_count;
}

const int64_t *TGSStats::GetHeadQualSum() const {
    return head_qual_sum;
}

const int64_t *TGSStats::GetTailQualSum() const {
    return tail_qual_sum;
}

const int TGSStats::GetMaxReadsLen() const {
    return mMaxReadsLen;
}

int *TGSStats::GetReadsLens() const {
    return readsLens;
}


const double TGSStats::GetAvgReadsLen() const {
    return mAvgReadsLen;
}

const int64_t TGSStats::GetReadsNum() const {
    return mReadsNum;
}

const int64_t TGSStats::GetBasesNum() const {
    return mBasesNum;
}

const int64_t *TGSStats::GetBases51015Num() const {
    return mBases51015Num;
}

const pair<double, int> *TGSStats::GetTop5QualReads() const {
    return mTop5QualReads;
}

const pair<int, double> *TGSStats::GetTop5LengReads() const {
    return mTop5LengReads;
}
