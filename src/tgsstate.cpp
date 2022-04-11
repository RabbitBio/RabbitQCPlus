//
// Created by ylf9811 on 2021/8/11.
//

#include "tgsstate.h"


TGSStats::TGSStats(int minLen) {
    mMinlen = minLen;
    mHalfMinlen = mMinlen >> 1;
    ///mLengths;
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
    mReadsNum=0;
    mBasesNum=0;
    mBases51015Num=new int64_t[3];
    mBases51015Num[0]=0;
    mBases51015Num[1]=0;
    mBases51015Num[2]=0;
    mTop5QualReads=new std::pair<double,int>[5];
    mTop5LengReads=new std::pair<int,double>[5];
    for(int i=0;i<5;i++){
        mTop5QualReads[i]={0,0};
        mTop5LengReads[i]={0,0};
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

void TGSStats::updateTop5(int rlen,double avgQual){
    int pos=5;
    for(int i=0;i<5;i++){
        if(mTop5QualReads[i].first<avgQual){
            for(int j=4;j>i;j--){
                mTop5QualReads[j]=mTop5QualReads[j-1];
            }
            mTop5QualReads[i]={avgQual,rlen};
            break;
        }
    }
    for(int i=0;i<5;i++){
        if(mTop5LengReads[i].first<rlen){
            for(int j=4;j>i;j--){
                mTop5LengReads[j]=mTop5LengReads[j-1];
            }
            mTop5LengReads[i]={rlen,avgQual};
            break;
        }
    }
}

void TGSStats::tgsStatRead(neoReference &ref, bool isPhred64) {
    const int rlen = ref.lseq;
    const char *seq = reinterpret_cast<const char *>(ref.base + ref.pseq);
    const char *quality = reinterpret_cast<const char *>(ref.base + ref.pqual);
    int size_range = mMinlen >> 1;
    size_range=min(size_range,rlen);
    int i;
    mReadsNum++;
    mBasesNum+=rlen;
    char c1, c2;
    mTotalReadsLen.push_back(rlen);
    int phredSub = 33;
    if (isPhred64)phredSub = 64;
    int sumQual=0;
    for(int i=0;i<rlen;i++){
        int qual=(quality[i] - phredSub);
        sumQual+=qual;
        if(qual>=5)mBases51015Num[0]++;
        if(qual>=10)mBases51015Num[1]++;
        if(qual>=15)mBases51015Num[2]++;
    }
    updateTop5(rlen,1.0*sumQual/rlen);
    //if (rlen > mMinlen) {
        //[1] stats lengths
        mLengths.push_back(rlen);
        //--head
        for (i = 0; i < size_range; ++i) {
            c1 = seq[i];
            //[2] quality sum
            head_qual_sum[i] += (quality[i] - phredSub);
            //[3] sequence count
            if (c1 == 'A') {
                head_seq_pos_count[0][i]++;
            } else if (c1 == 'T') {
                head_seq_pos_count[1][i]++;
            } else if (c1 == 'C') {
                head_seq_pos_count[2][i]++;
            } else if (c1 == 'G') {
                head_seq_pos_count[3][i]++;
            }

        }
        //--tail
        for (i = rlen - 1; i >= rlen - size_range; --i) {
            c2 = seq[i];
            tail_qual_sum[i - (rlen - size_range)] += (quality[i] - phredSub);

            if (c2 == 'A') {
                tail_seq_pos_count[0][i - (rlen - size_range)]++;
            } else if (c2 == 'T') {
                tail_seq_pos_count[1][i - (rlen - size_range)]++;
            } else if (c2 == 'C') {
                tail_seq_pos_count[2][i - (rlen - size_range)]++;
            } else if (c2 == 'G') {
                tail_seq_pos_count[3][i - (rlen - size_range)]++;
            }
        }
    //}
}

void TGSStats::print() {
    //cerr << "nothing here" << std::endl;
    int i;
    std::cout << "head A    T   C   G" << std::endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        std::cout << head_seq_pos_count[0][i] << " ";
        std::cout << head_seq_pos_count[1][i] << " ";
        std::cout << head_seq_pos_count[2][i] << " ";
        std::cout << head_seq_pos_count[3][i] << " ";
        std::cout << std::endl;
    }
    std::cout << "tail A    T   C   G" << std::endl;
    //std::cout << std::endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        std::cout << tail_seq_pos_count[0][i] << " ";
        std::cout << tail_seq_pos_count[1][i] << " ";
        std::cout << tail_seq_pos_count[2][i] << " ";
        std::cout << tail_seq_pos_count[3][i] << " ";
        std::cout << std::endl;
    }
    std::cout << "head and tail quality" << std::endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        std::cout << head_qual_sum[i] << " " << tail_qual_sum[i];
        std::cout << std::endl;
    }

}

TGSStats *TGSStats::merge(std::vector<TGSStats *> &list) {
    int i;
    const int minLen = list[0]->mMinlen;
    TGSStats *s = new TGSStats(minLen);
    for (TGSStats *ds : list) {
        //mLengths
        //s->mLengths.push_back(ds->mLengths);
        //ds->print();
        s->mLengths.insert(s->mLengths.end(), ds->mLengths.begin(), ds->mLengths.end());
        s->mTotalReadsLen.insert(s->mTotalReadsLen.end(), ds->mTotalReadsLen.begin(), ds->mTotalReadsLen.end());
        for (i = 0; i < (minLen >> 1); ++i) {
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
        for (i = 0; i < (minLen >> 1); ++i) {
            //head_qual
            s->head_qual_sum[i] += ds->head_qual_sum[i];
            //tail_qual
            s->tail_qual_sum[i] += ds->tail_qual_sum[i];
        }
        s->mReadsNum+=ds->mReadsNum;
        s->mBasesNum+=ds->mBasesNum;
        s->mBases51015Num[0]+=ds->mBases51015Num[0];
        s->mBases51015Num[1]+=ds->mBases51015Num[1];
        s->mBases51015Num[2]+=ds->mBases51015Num[2];
        for(int i=0;i<5;i++){
            s->updateTop5(ds->mTop5QualReads[i].second,ds->mTop5QualReads[i].first);
            s->updateTop5(ds->mTop5LengReads[i].first,ds->mTop5LengReads[i].second);
        }
        
    }

    
    //s->head_qual_mean = s->head_qual_sum/s->mLengths.size()
    return s;

}
void TGSStats::CalReadsLens(){
    int Ppre=10;
    mMaxReadsLen=0;
    int64_t readLenSum=0; 
    for(auto item:mTotalReadsLen){
        mMaxReadsLen=max(mMaxReadsLen,item);
        readLenSum+=item;
    }
    mAvgReadsLen=1.0*readLenSum/mTotalReadsLen.size();
    mMaxReadsLen/=Ppre;
    mMaxReadsLen+=Ppre;
    readsLens=new int[mMaxReadsLen];
    for(int i=0;i<mMaxReadsLen;i++)readsLens[i]=0;
    for(auto item:mTotalReadsLen){
        readsLens[(item+Ppre)/Ppre]++;
    }
}


//generate html data

bool TGSStats::isLongRead() {
    return mHalfMinlen > 300;
}

int TGSStats::base2num(std::string base) {
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


std::string TGSStats::list2string(double *list, int size) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

std::string TGSStats::list2string(double *list, int size, int64_t *coords) {
    std::stringstream ss;
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

std::string TGSStats::list2string(int64_t *list, int size) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

std::string TGSStats::list2stringReversedOrder(int64_t *list, int size) {
    std::stringstream ss;
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

const int TGSStats::GetMaxReadsLen() const{
    return mMaxReadsLen;
}

int *TGSStats::GetReadsLens() const {
    return readsLens;
}


const double TGSStats::GetAvgReadsLen() const{
    return mAvgReadsLen;
}

const int64_t TGSStats::GetReadsNum() const{
    return mReadsNum;
}

const int64_t TGSStats::GetBasesNum() const{
    return mBasesNum;
}

const int64_t* TGSStats::GetBases51015Num() const{
    return mBases51015Num;
}

const std::pair<double,int>* TGSStats::GetTop5QualReads() const{
    return mTop5QualReads;
}

const std::pair<int,double>* TGSStats::GetTop5LengReads() const{
    return mTop5LengReads;
}
