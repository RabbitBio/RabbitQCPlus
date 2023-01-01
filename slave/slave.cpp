#include <cstring>
#include <vector>
#include <unordered_map>

#include "tuple_spawn.hpp"
#include "slave.h"
#include "Reference.h"
#include "cmdinfo.h"
#include "state.h"
#include "tgsstate.h"
#include "threadinfo.h"
#define gcMax 100
#define LOCALSIZE 1 * (1 << 10)
#define BATCH_SIZE 2
__thread_local unsigned char local_data1[LOCALSIZE];
__thread_local unsigned char local_data2[LOCALSIZE];


struct OverlapRes {
    bool overlapped;
    int offset;
    int overlap_len;
    int diff_num;
};


struct dupInfo{
    uint32_t key;
    uint kmer32;
    uint8_t gc;
};

__thread_local dupInfo local_dup_info[BATCH_SIZE];

__thread_local int valAGCT[8]  = {-1, 0, -1, 2, 1, -1, -1, 3};

__thread_local char reMap[8] = {0, 'T', 0, 'G', 'A', 0, 0, 'C'};
int max(int x,int y) {
    return x > y ? x : y;
}

int min(int x,int y) {
    return x < y ? x : y;
}

void rtc_(unsigned long *counter){

    unsigned long rpcc;

    asm volatile("rcsr %0, 4":"=r"(rpcc));

    *counter=rpcc;

}


extern "C" void* ldm_malloc(size_t size);
extern "C" void ldm_free(void *addr, size_t size);
extern "C" size_t get_allocatable_size();

struct qc_data {
    ThreadInfo **thread_info_;
    CmdInfo *cmd_info_;
    std::vector <dupInfo> *dups;
    std::vector <neoReference> *data1_;
    std::vector <neoReference> *data2_;
    std::vector <neoReference> *pass_data1_;
    std::vector <neoReference> *pass_data2_;
    int *cnt;
    int bit_len;
};




void updateTop5(int rlen, double avgQual, TGSStats *tgs_state) {
    int pos = 5;
    for (int i = 0; i < 5; i++) {
        if (tgs_state->mTop5QualReads[i].first < avgQual) {
            for (int j = 4; j > i; j--) {
                tgs_state->mTop5QualReads[j] = tgs_state->mTop5QualReads[j - 1];
            }
            tgs_state->mTop5QualReads[i] = {avgQual, rlen};
            break;
        }
    }
    for (int i = 0; i < 5; i++) {
        if (tgs_state->mTop5LengReads[i].first < rlen) {
            for (int j = 4; j > i; j--) {
                tgs_state->mTop5LengReads[j] = tgs_state->mTop5LengReads[j - 1];
            }
            tgs_state->mTop5LengReads[i] = {rlen, avgQual};
            break;
        }
    }
}

void tgsfunc(qc_data *para) {
    std::vector <neoReference> *data = para->data1_;
    ThreadInfo *thread_info = para->thread_info_[_PEN];
    CmdInfo *cmd_info_ = para->cmd_info_;
    int data_num = data->size();
    int base_num = 0;
    // copy data to ldm
    int size_range_ = thread_info->TGS_state_->mMinlen >> 1;
    int64_t *head_qual_sum = (int64_t *)ldm_malloc(size_range_ * sizeof(int64_t));
    int64_t *tail_qual_sum = (int64_t *)ldm_malloc(size_range_ * sizeof(int64_t));
    dma_getn(thread_info->TGS_state_->head_qual_sum, head_qual_sum, size_range_);
    dma_getn(thread_info->TGS_state_->tail_qual_sum, tail_qual_sum, size_range_);

    int64_t *head_seq_pos_count[4];
    int64_t *tail_seq_pos_count[4];
    for(int i = 0; i < 4; i++){
        head_seq_pos_count[i] = (int64_t*)ldm_malloc(size_range_ * sizeof(int64_t));
        dma_getn(thread_info->TGS_state_->head_seq_pos_count[i], head_seq_pos_count[i], size_range_);
        tail_seq_pos_count[i] = (int64_t*)ldm_malloc(size_range_ * sizeof(int64_t));
        dma_getn(thread_info->TGS_state_->tail_seq_pos_count[i], tail_seq_pos_count[i], size_range_);
    }
    int64_t mBases51015Num[3];
    mBases51015Num[0] = thread_info->TGS_state_->mBases51015Num[0];
    mBases51015Num[1] = thread_info->TGS_state_->mBases51015Num[1];
    mBases51015Num[2] = thread_info->TGS_state_->mBases51015Num[2];
    int isPhred64_ = cmd_info_->isPhred64_;
    for(int id = _PEN * BATCH_SIZE; id < data_num; id += 64 * BATCH_SIZE) {
        int start_pos = id;
        int end_pos = id + BATCH_SIZE;
        if(end_pos > data_num) end_pos = data_num;
        for(int ids = start_pos; ids < end_pos; ids++) {
            auto ref = (*data)[ids];
            const int rlen = ref.lseq;
            const char *seq = reinterpret_cast<const char *>(ref.base + ref.pseq);
            const char *quality = reinterpret_cast<const char *>(ref.base + ref.pqual);
            int size_range = size_range_;
            if(rlen < size_range) size_range = rlen;
            if(thread_info->TGS_state_->countLenNum < (1 << 20))
                thread_info->TGS_state_->mTotalReadsLen[thread_info->TGS_state_->countLenNum++] = rlen;
            thread_info->TGS_state_->mReadsNum++;
            thread_info->TGS_state_->mBasesNum += rlen;
            char c1, c2;
            int phredSub = 33;
            if (isPhred64_) phredSub = 64;
            int sumQual = 0;
            for (int i = 0; i < rlen; i++) {
                int qual = (quality[i] - phredSub);
                sumQual += qual;
                if (qual >= 5) mBases51015Num[0]++;
                if (qual >= 10) mBases51015Num[1]++;
                if (qual >= 15) mBases51015Num[2]++;
            }
            updateTop5(rlen, 1.0 * sumQual / rlen, thread_info->TGS_state_);
            for (int i = 0; i < size_range; ++i) {
                c1 = seq[i];
                head_qual_sum[i] += (quality[i] - phredSub);
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
            for (int i = rlen - 1; i >= rlen - size_range; --i) {
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
        }

    }
    dma_putn(thread_info->TGS_state_->head_qual_sum, head_qual_sum, size_range_);
    dma_putn(thread_info->TGS_state_->tail_qual_sum, tail_qual_sum, size_range_);
    ldm_free(head_qual_sum, size_range_ * sizeof(int64_t));
    ldm_free(tail_qual_sum, size_range_ * sizeof(int64_t));
    for(int i = 0; i < 4; i++){
        dma_putn(thread_info->TGS_state_->head_seq_pos_count[i], head_seq_pos_count[i], size_range_);
        dma_putn(thread_info->TGS_state_->tail_seq_pos_count[i], tail_seq_pos_count[i], size_range_);
        ldm_free(head_seq_pos_count[i], size_range_ * sizeof(int64_t));
        ldm_free(tail_seq_pos_count[i], size_range_ * sizeof(int64_t));
    }

    dma_putn(thread_info->TGS_state_->mBases51015Num, mBases51015Num, 3);
}



void StateInfo(State *state, neoReference &ref, CmdInfo *cmd_info_, int *pos_cnt_, int *pos_qul_, int *len_cnt_, int *gc_cnt_, int *qul_cnt_, int &q20bases_, int &q30bases_, int &tot_bases_, int &gc_bases_, int &real_seq_len_, int &lines_){
    int slen = ref.lseq;
    int qlen = ref.lqual;
    //ASSERT(slen == qlen);
    //if (slen > malloc_seq_len_) {
    //    ExtendBuffer(malloc_seq_len_, max(slen + 100, slen * 2));
    //}
    //TODO
    if(slen > real_seq_len_) real_seq_len_ = slen;
    len_cnt_[slen - 1]++;
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);
    int flag = 4;
    int gc_cnt = 0;
    int qul_tot = 0;
    int kmer = 0;
    int phredSub = 33;
    if (cmd_info_->isPhred64_) phredSub = 64;
    for (int i = 0; i < slen; i++) {
        char b = bases[i] & 0x07;
        if (b == 3 || b == 7) gc_cnt++;
        int q = quals[i] - phredSub;
        if(q < 0) q = 0;
        qul_tot += q;
        if (q >= 30) {
            q20bases_++;
            q30bases_++;
        } else if (q >= 20) {
            q20bases_++;
        }
        pos_cnt_[i * 4 + valAGCT[b]]++;
        pos_qul_[i] += q;
        //if (bases[i] == 'N') flag = 5;
        //int val = valAGCT[b];
        //kmer = ((kmer << 2) & 0x3FC) | val;
        //if (flag <= 0) state->kmer_[kmer]++;
        //flag--;
    }
    tot_bases_ += slen;
    gc_bases_ += gc_cnt;
    gc_cnt_[int(1.0 * gcMax * gc_cnt / slen)]++;
    qul_cnt_[int(1.0 * qul_tot / slen)]++;
    // do overrepresentation analysis for 1 of every 20 reads
    //if (do_over_represent_analyze_) {
    //    if (lines_ % over_representation_sampling_ == 0) {
    //        //StateORP(ref);
    //        //TODO
    //    }
    //}
    lines_++;
}

bool TrimSeq(neoReference &ref, neoReference &ref2, int front, int tail, CmdInfo *cmd_info_) {
    int new_seq_len = ref.lseq - front - tail;
    if (new_seq_len <= 0) return false;

    int w = cmd_info_->cut_window_size_;
    int l = ref.lseq;
    const char *seq = reinterpret_cast<const char *>(ref.base + ref.pseq);
    const char *qualstr = reinterpret_cast<const char *>(ref.base + ref.pqual);
    int phredSub = 33;
    int cut_mean_quality_ = cmd_info_->cut_mean_quality_;
    if (cmd_info_->isPhred64_) phredSub = 64;
    // quality cutting forward
    if (cmd_info_->trim_5end_) {
        int s = front;
        if (l - front - tail - w <= 0) {
            return 0;
        }

        int totalQual = 0;

        // preparing rolling
        for (int i = 0; i < w - 1; i++)
            totalQual += qualstr[s + i];

        for (s = front; s + w < l - tail; s++) {
            totalQual += qualstr[s + w - 1];
            // rolling
            if (s > front) {
                totalQual -= qualstr[s - 1];
            }
            // add 33 for phred33 transforming
            if ((double) totalQual / (double) w >= phredSub + cut_mean_quality_)
                break;
        }

        // the trimming in front is forwarded and rlen is recalculated
        if (s > 0)
            s = s + w - 1;
        while (s < l && seq[s] == 'N')
            s++;
        front = s;
        new_seq_len = l - front - tail;
    }

    // quality cutting backward
    if (cmd_info_->trim_3end_) {
        if (l - front - tail - w <= 0) {
            return 0;
        }

        int totalQual = 0;
        int t = l - tail - 1;

        // preparing rolling
        for (int i = 0; i < w - 1; i++)
            totalQual += qualstr[t - i];

        for (t = l - tail - 1; t - w >= front; t--) {
            totalQual += qualstr[t - w + 1];
            // rolling
            if (t < l - tail - 1) {
                totalQual -= qualstr[t + 1];
            }
            // add 33 for phred33 transforming
            if ((double) totalQual / (double) w >= phredSub + cut_mean_quality_)
                break;
        }

        if (t < l - 1)
            t = t - w + 1;
        while (t >= 0 && seq[t] == 'N')
            t--;
        new_seq_len = t - front + 1;
    }

    if (new_seq_len <= 0 || front >= l - 1) {
        return 0;
    }


    ref.pseq += front;
    ref.pqual += front;
    ref.lseq = new_seq_len;
    ref.lqual = new_seq_len;

    ref2.pseq += front;
    ref2.pqual += front;
    ref2.lseq = new_seq_len;
    ref2.lqual = new_seq_len;

    return true;

}

int ReadFiltering(neoReference &ref, bool trim_res, bool isPhred64, CmdInfo *cmd_info_) {
    if (!trim_res) return 2;
    int seq_len = ref.lseq;
    int qul_len = ref.lqual;
    //ASSERT(seq_len == qul_len);
    if (seq_len >= cmd_info_->length_limit_) {
        return 3;
    }
    if (seq_len == 0 || seq_len < cmd_info_->length_required_) {
        return 2;
    }
    int n_number = 0;
    int low_qual_number = 0;// < q15
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);

    int phredSub = 33;
    if (isPhred64) phredSub = 64;
    for (int i = 0; i < seq_len; i++) {
        int q = std::max(0, quals[i] - phredSub);
        //int q = quals[i] - phredSub;
        if (q < 15) {
            low_qual_number++;
        }
        if (bases[i] == 'N') {
            n_number++;
        }
    }
    if (low_qual_number > cmd_info_->low_qual_perc_limit_ * seq_len / 100.0) return 4;
    if (n_number > cmd_info_->n_number_limit_) return 1;
    return 0;
}

int TrimAdapterPE(neoReference &r1, neoReference &r2, neoReference &r1p, neoReference &r2p, int offset, int overlap_len) {
    //    if(ov.diff<=5 && ov.overlapped && ov.offset < 0 && ol > r1->length()/3)

    if (overlap_len > 0 && offset < 0) {
        //string adapter1 = string(reinterpret_cast<const char *>(r1.base + r1.pseq + overlap_len),
        //        r1.lseq - overlap_len);
        //string adapter2 = string(reinterpret_cast<const char *>(r2.base + r2.pseq + overlap_len),
        //        r2.lseq - overlap_len);
        int res_len = r1.lseq - overlap_len + r2.lseq - overlap_len;
        r1.lseq = overlap_len;
        r1p.lseq = overlap_len;
        r1.lqual = overlap_len;
        r1p.lqual = overlap_len;
        r2.lseq = overlap_len;
        r2p.lseq = overlap_len;
        r2.lqual = overlap_len;
        r2p.lqual = overlap_len;

        return res_len;
        //return adapter1.length() + adapter2.length();
    }
    return false;
}




int TrimAdapter(neoReference &ref, neoReference &refp, std::string &adapter_seq, bool isR2) {
    int matchReq = 4;
    const int allowOneMismatchForEach = 8;

    int rlen = ref.lseq;
    int alen = adapter_seq.length();

    const char *adata = adapter_seq.c_str();
    const char *rdata = reinterpret_cast<const char *>(ref.base + ref.pseq);

    if (alen < matchReq)
        return false;

    int pos = 0;
    bool found = false;
    int start = 0;
    if (alen >= 16)
        start = -4;
    else if (alen >= 12)
        start = -3;
    else if (alen >= 8)
        start = -2;

    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for (int i = max(0, -pos); i < cmplen; i++) {
            mismatch += adata[i] != rdata[i + pos];
            if (mismatch > allowedMismatch) break;
        }
        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }

    if (found) {
        int res_len = 0;
        if (pos < 0) {
            //std::string adapter = adapter_seq.substr(0, alen + pos);
            ref.lseq = 0;
            refp.lseq = 0;
            ref.lqual = 0;
            refp.lqual = 0;
            res_len = alen + pos;

        } else {
            //std::string new_adapter_seq = std::string(reinterpret_cast<const char *>(ref.base + ref.pseq + pos),
            //        rlen - pos);
            ref.lseq = pos;
            refp.lseq = pos;
            ref.lqual = pos;
            refp.lqual = pos;
            res_len = rlen - pos;
        }
        return res_len;
    }

    return false;
}


void trimPolyG(neoReference &ref, neoReference &refp, int compareReq) {
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;

    const char *data = reinterpret_cast<const char *>(ref.base + ref.pseq);

    int rlen = ref.lseq;

    int mismatch = 0;
    int i = 0;
    int firstGPos = rlen - 1;
    for (i = 0; i < rlen; i++) {
        if (data[rlen - i - 1] != 'G') {
            mismatch++;
        } else {
            firstGPos = rlen - i - 1;
        }

        int allowedMismatch = (i + 1) / allowOneMismatchForEach;
        if (mismatch > maxMismatch || (mismatch > allowedMismatch && i >= compareReq - 1))
            break;
    }

    if (i >= compareReq) {
        ref.lseq = firstGPos;
        refp.lseq = firstGPos;
        ref.lqual = firstGPos;
        refp.lqual = firstGPos;
    }
}
void trimPolyX(neoReference &ref, neoReference &refp, int compareReq) {
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;

    const char *data = reinterpret_cast<const char *>(ref.base + ref.pseq);

    int rlen = ref.lseq;

    int atcgNumbers[4] = {0, 0, 0, 0};
    char atcgBases[4] = {'A', 'T', 'C', 'G'};
    int pos = 0;
    for (pos = 0; pos < rlen; pos++) {
        switch (data[rlen - pos - 1]) {
            case 'A':
                atcgNumbers[0]++;
                break;
            case 'T':
                atcgNumbers[1]++;
                break;
            case 'C':
                atcgNumbers[2]++;
                break;
            case 'G':
                atcgNumbers[3]++;
                break;
            case 'N':
                atcgNumbers[0]++;
                atcgNumbers[1]++;
                atcgNumbers[2]++;
                atcgNumbers[3]++;
                break;
            default:
                break;
        }

        int cmp = (pos + 1);
        int allowedMismatch = min(maxMismatch, cmp / allowOneMismatchForEach);

        bool needToBreak = true;
        for (int b = 0; b < 4; b++) {
            if (cmp - atcgNumbers[b] <= allowedMismatch)
                needToBreak = false;
        }
        if (needToBreak && (pos >= allowOneMismatchForEach || pos + 1 >= compareReq - 1)) {
            break;
        }
    }


    // has polyX
    if (pos + 1 >= compareReq && pos < rlen) {
        // find the poly
        int poly;
        int maxCount = -1;
        for (int b = 0; b < 4; b++) {
            if (atcgNumbers[b] > maxCount) {
                maxCount = atcgNumbers[b];
                poly = b;
            }
        }
        char polyBase = atcgBases[poly];
        while (data[rlen - pos - 1] != polyBase && pos >= 0)
            pos--;
        ref.lseq = rlen - pos - 1;
        refp.lseq = rlen - pos - 1;
        ref.lqual = rlen - pos - 1;
        refp.lqual = rlen - pos - 1;
    }
}

void trimPolyGPE(neoReference &r1, neoReference &r2, neoReference &r1p, neoReference &r2p, int compareReq) {
    trimPolyG(r1, r1p, compareReq);
    trimPolyG(r2, r2p, compareReq);
}

void trimPolyXPE(neoReference &r1, neoReference &r2, neoReference &r1p, neoReference &r2p, int compareReq) {
    trimPolyX(r1, r1p, compareReq);
    trimPolyX(r2, r2p, compareReq);
}



uint seq2int(const char *data, int start, int key_len, bool &valid) {
    uint ret = 0;
    for (int i = 0; i < key_len; i++) {
        ret <<= 2;
        if (data[start + i] == 'N') {
            valid = false;
            return 0;
        }
        ret += valAGCT[data[start + i] & 0x07];
        // if it's not the last one, shift it by 2 bits
    }
    return ret;
}

void StateDup(neoReference &ref, int ids, int bit_len) {
    if (ref.lseq < 32)
        return;
    int start1 = 0;
    uint zr = 0;
    int start2 = max(zr, ref.lseq - 32 - 5);

    const char *data = reinterpret_cast<const char *>(ref.base + ref.pseq);
    bool valid = true;

    uint ret = 0;
    ret = seq2int(data, start1, bit_len, valid);
    uint32_t key = (uint32_t) ret;
    if (!valid)
        return;

    uint kmer32 = 0;
    kmer32 = seq2int(data, start2, 32, valid);
    if (!valid)
        return;

    int gc = 0;
    for (int i = 0; i < ref.lseq; i++) {
        if (data[i] == 'C' || data[i] == 'G')
            gc++;
    }
    gc = int(255.0 * gc / ref.lseq + 0.5);
    local_dup_info[ids & (BATCH_SIZE - 1)].key = key;
    local_dup_info[ids & (BATCH_SIZE - 1)].kmer32 = kmer32;
    local_dup_info[ids & (BATCH_SIZE - 1)].gc = (uint8_t)gc;
}

void ngsfunc(qc_data *para){
    unsigned long c_trim = 0;
    unsigned long c_dup = 0;
    unsigned long c_loop = 0;
    unsigned long c_dma_all = 0;
    unsigned long c_dma2 = 0;
    unsigned long c_dma3 = 0;
    unsigned long c_dma4 = 0;
    unsigned long c_dma5 = 0;
    unsigned long c_state1 = 0;
    unsigned long c_state2 = 0;
    unsigned long c_write1 = 0;
    unsigned long c_write2 = 0;
    unsigned long c_dma_dup = 0;
    unsigned long c_filter = 0;
    unsigned long start_t, end_t;
    unsigned long start_tt, end_tt;
    rtc_(&start_tt);
    std::vector <neoReference> *data = para->data1_;
    std::vector <neoReference> *pass_data = para->pass_data1_;
    ThreadInfo *thread_info = para->thread_info_[_PEN];
    CmdInfo *cmd_info_ = para->cmd_info_;
    int bit_len = para->bit_len;
    int data_num = data->size();

    int pre_pos_seq_len = thread_info->pre_state1_->malloc_seq_len_ * 4;
    int *pre_pos_cnt_ = (int*)ldm_malloc(pre_pos_seq_len * sizeof(int));
    dma_getn(thread_info->pre_state1_->pos_cnt_, pre_pos_cnt_, pre_pos_seq_len);

    int aft_pos_seq_len = thread_info->aft_state1_->malloc_seq_len_ * 4;
    int *aft_pos_cnt_ = (int*)ldm_malloc(aft_pos_seq_len * sizeof(int));
    dma_getn(thread_info->aft_state1_->pos_cnt_, aft_pos_cnt_, aft_pos_seq_len);


    int pre_malloc_seq_len_ = thread_info->pre_state1_->malloc_seq_len_;
    int *pre_len_cnt_ = (int*)ldm_malloc(pre_malloc_seq_len_ * sizeof(int));
    dma_getn(thread_info->pre_state1_->len_cnt_, pre_len_cnt_, pre_malloc_seq_len_);
    int *pre_pos_qul_ = (int*)ldm_malloc(pre_malloc_seq_len_ * sizeof(int));
    dma_getn(thread_info->pre_state1_->pos_qul_, pre_pos_qul_, pre_malloc_seq_len_);



    int aft_malloc_seq_len_ = thread_info->aft_state1_->malloc_seq_len_;
    int *aft_len_cnt_ = (int*)ldm_malloc(aft_malloc_seq_len_ * sizeof(int));
    dma_getn(thread_info->aft_state1_->len_cnt_, aft_len_cnt_, aft_malloc_seq_len_);
    int *aft_pos_qul_ = (int*)ldm_malloc(aft_malloc_seq_len_ * sizeof(int));
    dma_getn(thread_info->aft_state1_->pos_qul_, aft_pos_qul_, aft_malloc_seq_len_);

    int *pre_gc_cnt_ = (int*)ldm_malloc((gcMax + 1) * sizeof(int));
    dma_getn(thread_info->pre_state1_->gc_cnt_, pre_gc_cnt_, gcMax + 1);

    int *aft_gc_cnt_ = (int*)ldm_malloc((gcMax + 1) * sizeof(int));
    dma_getn(thread_info->aft_state1_->gc_cnt_, aft_gc_cnt_, gcMax + 1);

    int pre_qul_range_ = thread_info->pre_state1_->qul_range_;
    int *pre_qul_cnt_ = (int*)ldm_malloc(pre_qul_range_ * sizeof(int));
    dma_getn(thread_info->pre_state1_->qul_cnt_, pre_qul_cnt_, pre_qul_range_);

    int aft_qul_range_ = thread_info->aft_state1_->qul_range_;
    int *aft_qul_cnt_ = (int*)ldm_malloc(aft_qul_range_ * sizeof(int));
    dma_getn(thread_info->aft_state1_->qul_cnt_, aft_qul_cnt_, aft_qul_range_);
    rtc_(&end_tt);
    c_dma_all += end_tt - start_tt;

    rtc_(&start_tt);
    int pre_q20bases_ = 0;
    int pre_q30bases_ = 0;
    int pre_tot_bases_ = 0;
    int pre_gc_bases_ = 0;
    int pre_real_seq_len_ = 0;
    int pre_lines_ = 0;

    int aft_q20bases_ = 0;
    int aft_q30bases_ = 0;
    int aft_tot_bases_ = 0;
    int aft_gc_bases_ = 0;
    int aft_real_seq_len_ = 0;
    int aft_lines_ = 0;
    for(int id = _PEN * BATCH_SIZE; id < data_num; id += 64 * BATCH_SIZE) {


        rtc_(&start_t);
        int start_pos = id;
        int end_pos = id + BATCH_SIZE;
        if(end_pos > data_num) end_pos = data_num;
        int cal_size = (*data)[end_pos - 1].pqual + (*data)[end_pos - 1].lqual + 1 - (*data)[start_pos].pname;
        auto item0 = (*data)[start_pos];
        dma_getn((unsigned char*)(item0.base + item0.pname), local_data1, cal_size);
        rtc_(&end_t);
        c_dma2 += end_t - start_t;



        for(int ids = start_pos; ids < end_pos; ids++) {



            rtc_(&start_t);
            neoReference item;
            dma_getn(&((*data)[ids]), &item, 1);
            neoReference item2;
            dma_getn(&((*data)[ids]), &item2, 1);
            int now_pos = item.pname - item0.pname;
            item.base = &(local_data1[now_pos]) - item.pname;
            rtc_(&end_t);
            c_dma3 += end_t - start_t;




            rtc_(&start_t);
            StateInfo(thread_info->pre_state1_, item, cmd_info_, pre_pos_cnt_, pre_pos_qul_, pre_len_cnt_, pre_gc_cnt_, pre_qul_cnt_, pre_q20bases_, pre_q30bases_, pre_tot_bases_, pre_gc_bases_, pre_real_seq_len_, pre_lines_); 
            rtc_(&end_t);
            c_state1 += end_t - start_t;




            rtc_(&start_t);
            if (cmd_info_->state_duplicate_) {
                StateDup(item, ids, bit_len);
            }
            if (cmd_info_->add_umi_) {
                //umier_->ProcessSe(item);
                //TODO
            }
            rtc_(&end_t);
            c_dup += end_t - start_t;

            rtc_(&start_t);
            bool trim_res = 0;
            trim_res = TrimSeq(item, item2, cmd_info_->trim_front2_, cmd_info_->trim_tail1_, cmd_info_);
            if (trim_res && cmd_info_->trim_polyg_) {
                trimPolyG(item, item2, cmd_info_->trim_poly_len_);
            }
            if (trim_res && cmd_info_->trim_polyx_) {
                trimPolyX(item, item2, cmd_info_->trim_poly_len_);
            }
            if (trim_res && cmd_info_->trim_adapter_) {
                int res = 0;
                bool is_trimmed = false;
                if (cmd_info_->detect_adapter1_) {
                    if (cmd_info_->print_what_trimmed_) {
                        //res = TrimAdapter(item, cmd_info_->adapter_seq1_, thread_info->aft_state1_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                        //TODO
                    } else {
                        res = TrimAdapter(item, item2, cmd_info_->adapter_seq1_, false);
                    }
                    if (res) {
                        is_trimmed = true;
                        thread_info->aft_state1_->trim_adapter_++;
                        thread_info->aft_state1_->trim_adapter_bases_ += res;
                    }
                }
                if (cmd_info_->adapter_from_fasta_.size() > 0) {
                    if (cmd_info_->print_what_trimmed_) {
                        //res = TrimAdapters(item, cmd_info_->adapter_from_fasta_, thread_info->aft_state1_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                        //TODO
                    } else {
                        //res = TrimAdapters(item, cmd_info_->adapter_from_fasta_, false);
                        //TODO
                    }
                    if (res) {
                        if (!is_trimmed) thread_info->aft_state1_->trim_adapter_++;
                        thread_info->aft_state1_->trim_adapter_bases_ += res;
                    }
                }
            }
            rtc_(&end_t);
            c_trim += end_t - start_t;



            int filter_res = 0;
            rtc_(&start_t);
            filter_res = ReadFiltering(item, trim_res, cmd_info_->isPhred64_, cmd_info_);
            rtc_(&end_t);
            c_filter += end_t - start_t;
            if (filter_res == 0) {
                rtc_(&start_t);
                StateInfo(thread_info->aft_state1_, item, cmd_info_, aft_pos_cnt_, aft_pos_qul_, aft_len_cnt_, aft_gc_cnt_, aft_qul_cnt_, aft_q20bases_, aft_q30bases_, aft_tot_bases_, aft_gc_bases_, aft_real_seq_len_, aft_lines_); 
                rtc_(&end_t);
                c_state2 += end_t - start_t;



                rtc_(&start_t);
                thread_info->aft_state1_->pass_reads_++;
                if (cmd_info_->write_data_) {
                    dma_putn(&((*pass_data)[ids]), &item2, 1);
                    //(*pass_data)[ids] = item2;
                }
                rtc_(&end_t);
                c_dma4 += end_t - start_t;




            } else {
                rtc_(&start_t);
                if (cmd_info_->write_data_) {
                    item2.lname = 0;
                    item2.lseq = 0;
                    item2.lstrand = 0;
                    item2.lqual = 0;
                    dma_putn(&((*pass_data)[ids]), &item2, 1);
                    //(*pass_data)[ids] = item2;
                }
                rtc_(&end_t);
                c_write1 += end_t - start_t;




                rtc_(&start_t);
                if (filter_res == 1) {
                    thread_info->aft_state1_->fail_N_++;
                } else if (filter_res == 2) {
                    thread_info->aft_state1_->fail_short_++;
                } else if (filter_res == 3) {
                    thread_info->aft_state1_->fail_long_++;
                } else if (filter_res == 4) {
                    thread_info->aft_state1_->fail_lowq_++;
                }
                rtc_(&end_t);
                c_write2 += end_t - start_t;
            }

        }
        rtc_(&start_t);
        if(cmd_info_->state_duplicate_) {
            dma_putn(&((*(para->dups))[start_pos]), local_dup_info, BATCH_SIZE);
        }
        rtc_(&end_t);
        c_dma_dup += end_t - start_t;
    }
    rtc_(&end_tt);
    c_loop += end_tt - start_tt;
    rtc_(&start_tt);
    dma_putn(thread_info->pre_state1_->pos_cnt_, pre_pos_cnt_, pre_pos_seq_len);
    dma_putn(thread_info->aft_state1_->pos_cnt_, aft_pos_cnt_, aft_pos_seq_len);
    ldm_free(pre_pos_cnt_, pre_pos_seq_len * sizeof(int));
    ldm_free(aft_pos_cnt_, aft_pos_seq_len * sizeof(int));

    dma_putn(thread_info->pre_state1_->pos_qul_, pre_pos_qul_, pre_malloc_seq_len_);
    dma_putn(thread_info->aft_state1_->pos_qul_, aft_pos_qul_, aft_malloc_seq_len_);
    dma_putn(thread_info->pre_state1_->len_cnt_, pre_len_cnt_, pre_malloc_seq_len_);
    dma_putn(thread_info->aft_state1_->len_cnt_, aft_len_cnt_, aft_malloc_seq_len_);
    dma_putn(thread_info->pre_state1_->gc_cnt_, pre_gc_cnt_, gcMax + 1);
    dma_putn(thread_info->aft_state1_->gc_cnt_, aft_gc_cnt_, gcMax + 1);
    dma_putn(thread_info->pre_state1_->qul_cnt_, pre_qul_cnt_, pre_qul_range_);
    dma_putn(thread_info->aft_state1_->qul_cnt_, aft_qul_cnt_, aft_qul_range_);
    ldm_free(pre_pos_qul_, pre_malloc_seq_len_ * sizeof(int));
    ldm_free(aft_pos_qul_, aft_malloc_seq_len_ * sizeof(int));
    ldm_free(pre_len_cnt_, pre_malloc_seq_len_ * sizeof(int));
    ldm_free(aft_len_cnt_, aft_malloc_seq_len_ * sizeof(int));
    ldm_free(pre_gc_cnt_, (gcMax + 1) * sizeof(int));
    ldm_free(aft_gc_cnt_, (gcMax + 1) * sizeof(int));
    ldm_free(pre_qul_cnt_, pre_qul_range_ * sizeof(int));
    ldm_free(aft_qul_cnt_, aft_qul_range_ * sizeof(int));

    thread_info->pre_state1_->q20bases_ += pre_q20bases_;
    thread_info->pre_state1_->q30bases_ += pre_q30bases_;
    thread_info->pre_state1_->tot_bases_ += pre_tot_bases_;
    thread_info->pre_state1_->gc_bases_ += pre_gc_bases_;
    thread_info->pre_state1_->real_seq_len_ = pre_real_seq_len_;
    thread_info->pre_state1_->lines_ += pre_lines_;

    thread_info->aft_state1_->q20bases_ += aft_q20bases_;
    thread_info->aft_state1_->q30bases_ += aft_q30bases_;
    thread_info->aft_state1_->tot_bases_ += aft_tot_bases_;
    thread_info->aft_state1_->gc_bases_ += aft_gc_bases_;
    thread_info->aft_state1_->real_seq_len_ = aft_real_seq_len_;
    thread_info->aft_state1_->lines_ += aft_lines_;
    rtc_(&end_tt);
    c_dma5 += end_tt - start_tt;
    if(_PEN == 0) {
        printf("== dma all %lf\n", c_dma_all * 1e-6);
        printf("== loop %lf\n", c_loop * 1e-6);
        printf("== dma5 %lf\n", c_dma5 * 1e-6);

        printf("== dma2 %lf\n", c_dma2 * 1e-6);
        printf("== dma3 %lf\n", c_dma3 * 1e-6);
        printf("== state1 %lf\n", c_state1 * 1e-6);
        printf("== dup %lf\n", c_dup * 1e-6);
        printf("== trim %lf\n", c_trim * 1e-6);
        printf("== filter %lf\n", c_filter * 1e-6);
        printf("== state2 %lf\n", c_state2 * 1e-6);
        printf("== dma4 %lf\n", c_dma4 * 1e-6);
        printf("== write1 %lf\n", c_write1 * 1e-6);
        printf("== write2 %lf\n", c_write2 * 1e-6);
        printf("== dma_dup %lf\n", c_dma_dup * 1e-6);
    }
}

OverlapRes AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require) {


    int len1 = r1.lseq;
    int len2 = r2.lseq;
    const char *str1 = reinterpret_cast<const char *>(r1.base + r1.pseq);
    const char *str2 = reinterpret_cast<const char *>(r2.base + r2.pseq);


    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff_num = 0;

    // forward
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1 - overlap_require) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);
        int leng = min(complete_compare_require, overlap_len);
        diff_num = 0;
        for (int i = 0; i < leng; i++) {
            diff_num += str1[offset + i] != reMap[str2[len2 - 1 - i] & 0x07];
            if (diff_num > overlap_diff_limit)
                break;
        }

        if (diff_num <= overlap_diff_limit) {
            return {true, offset, overlap_len, diff_num};
        }

        offset += 1;
    }


    // reverse
    // in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    // check if distance can get smaller if offset goes negative
    // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    // we go reversely
    offset = 0;
    while (offset > -(len2 - overlap_require)) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1, len2 - abs(offset));
        int leng = min(complete_compare_require, overlap_len);
        diff_num = 0;
        for (int i = 0; i < leng; i++) {
            diff_num += str1[i] != reMap[str2[len2 - 1 + offset - i] & 0x07];
            if (diff_num > overlap_diff_limit)
                break;
        }

        if (diff_num <= overlap_diff_limit) {
            return {true, offset, overlap_len, diff_num};
        }

        offset -= 1;
    }

    return {false, 0, 0, 0};
}


void PrintString(char *pos, int len) {
    for(int i = 0; i < len; i++) {
        printf("%c", pos[i]);
    }
    printf("\n");
}

void UpdateIterm(neoReference &up, neoReference now) {
    up.pname = now.pname;
    up.lname = now.lname;
    up.pseq = now.pseq;
    up.lseq = now.lseq;
    up.pstrand = now.pstrand;
    up.lstrand = now.lstrand;
    up.pqual = now.pqual;
    up.lqual = now.lqual;
}

void ngspefunc(qc_data *para){
    unsigned long c_trim = 0;
    unsigned long c_dup = 0;
    unsigned long c_loop = 0;
    unsigned long c_dma_all = 0;
    unsigned long c_dma2 = 0;
    unsigned long c_dma3 = 0;
    unsigned long c_dma4 = 0;
    unsigned long c_dma5 = 0;
    unsigned long c_state1 = 0;
    unsigned long c_state2 = 0;
    unsigned long c_write1 = 0;
    unsigned long c_write2 = 0;
    unsigned long c_dma_dup = 0;
    unsigned long c_filter = 0;
    unsigned long start_t, end_t;
    unsigned long start_tt, end_tt;
    rtc_(&start_tt);
    std::vector <neoReference> *data1 = para->data1_;
    std::vector <neoReference> *data2 = para->data2_;
    std::vector <neoReference> *pass_data1 = para->pass_data2_;
    std::vector <neoReference> *pass_data2 = para->pass_data1_;
    ThreadInfo *thread_info = para->thread_info_[_PEN];
    CmdInfo *cmd_info_ = para->cmd_info_;
    int bit_len = para->bit_len;

    int data_num1 = data1->size();
    int data_num2 = data2->size();


    int pre_pos_seq_len1 = thread_info->pre_state1_->malloc_seq_len_ * 4;
    int *pre_pos_cnt1_ = (int*)ldm_malloc(pre_pos_seq_len1 * sizeof(int));
    dma_getn(thread_info->pre_state1_->pos_cnt_, pre_pos_cnt1_, pre_pos_seq_len1);

    int aft_pos_seq_len1 = thread_info->aft_state1_->malloc_seq_len_ * 4;
    int *aft_pos_cnt1_ = (int*)ldm_malloc(aft_pos_seq_len1 * sizeof(int));
    dma_getn(thread_info->aft_state1_->pos_cnt_, aft_pos_cnt1_, aft_pos_seq_len1);


    int pre_malloc_seq_len1_ = thread_info->pre_state1_->malloc_seq_len_;
    int *pre_pos_qul1_ = (int*)ldm_malloc(pre_malloc_seq_len1_ * sizeof(int));
    dma_getn(thread_info->pre_state1_->pos_qul_, pre_pos_qul1_, pre_malloc_seq_len1_);
    int *pre_len_cnt1_ = (int*)ldm_malloc(pre_malloc_seq_len1_ * sizeof(int));
    dma_getn(thread_info->pre_state1_->len_cnt_, pre_len_cnt1_, pre_malloc_seq_len1_);

    int aft_malloc_seq_len1_ = thread_info->aft_state1_->malloc_seq_len_;
    int *aft_pos_qul1_ = (int*)ldm_malloc(aft_malloc_seq_len1_ * sizeof(int));
    dma_getn(thread_info->aft_state1_->pos_qul_, aft_pos_qul1_, aft_malloc_seq_len1_);
    int *aft_len_cnt1_ = (int*)ldm_malloc(aft_malloc_seq_len1_ * sizeof(int));
    dma_getn(thread_info->aft_state1_->len_cnt_, aft_len_cnt1_, aft_malloc_seq_len1_);

    int *pre_gc_cnt1_ = (int*)ldm_malloc((gcMax + 1) * sizeof(int));
    dma_getn(thread_info->pre_state1_->gc_cnt_, pre_gc_cnt1_, gcMax + 1);

    int *aft_gc_cnt1_ = (int*)ldm_malloc((gcMax + 1) * sizeof(int));
    dma_getn(thread_info->aft_state1_->gc_cnt_, aft_gc_cnt1_, gcMax + 1);

    int pre_qul_range1_ = thread_info->pre_state1_->qul_range_;
    int *pre_qul_cnt1_ = (int*)ldm_malloc(pre_qul_range1_ * sizeof(int));
    dma_getn(thread_info->pre_state1_->qul_cnt_, pre_qul_cnt1_, pre_qul_range1_);

    int aft_qul_range1_ = thread_info->aft_state1_->qul_range_;
    int *aft_qul_cnt1_ = (int*)ldm_malloc(aft_qul_range1_ * sizeof(int));
    dma_getn(thread_info->aft_state1_->qul_cnt_, aft_qul_cnt1_, aft_qul_range1_);







    int pre_pos_seq_len2 = thread_info->pre_state2_->malloc_seq_len_ * 4;
    int *pre_pos_cnt2_ = (int*)ldm_malloc(pre_pos_seq_len2 * sizeof(int));
    dma_getn(thread_info->pre_state2_->pos_cnt_, pre_pos_cnt2_, pre_pos_seq_len2);

    int aft_pos_seq_len2 = thread_info->aft_state2_->malloc_seq_len_ * 4;
    int *aft_pos_cnt2_ = (int*)ldm_malloc(aft_pos_seq_len2 * sizeof(int));
    dma_getn(thread_info->aft_state2_->pos_cnt_, aft_pos_cnt2_, aft_pos_seq_len2);

    int pre_malloc_seq_len2_ = thread_info->pre_state2_->malloc_seq_len_;
    int *pre_pos_qul2_ = (int*)ldm_malloc(pre_malloc_seq_len2_ * sizeof(int));
    dma_getn(thread_info->pre_state2_->pos_qul_, pre_pos_qul2_, pre_malloc_seq_len2_);
    int *pre_len_cnt2_ = (int*)ldm_malloc(pre_malloc_seq_len2_ * sizeof(int));
    dma_getn(thread_info->pre_state2_->len_cnt_, pre_len_cnt2_, pre_malloc_seq_len2_);


    int aft_malloc_seq_len2_ = thread_info->aft_state2_->malloc_seq_len_;
    int *aft_pos_qul2_ = (int*)ldm_malloc(aft_malloc_seq_len2_ * sizeof(int));
    dma_getn(thread_info->aft_state2_->pos_qul_, aft_pos_qul2_, aft_malloc_seq_len2_);
    int *aft_len_cnt2_ = (int*)ldm_malloc(aft_malloc_seq_len2_ * sizeof(int));
    dma_getn(thread_info->aft_state2_->len_cnt_, aft_len_cnt2_, aft_malloc_seq_len2_);


    int *pre_gc_cnt2_ = (int*)ldm_malloc((gcMax + 1) * sizeof(int));
    dma_getn(thread_info->pre_state2_->gc_cnt_, pre_gc_cnt2_, gcMax + 1);


    int *aft_gc_cnt2_ = (int*)ldm_malloc((gcMax + 1) * sizeof(int));
    dma_getn(thread_info->aft_state2_->gc_cnt_, aft_gc_cnt2_, gcMax + 1);


    int pre_qul_range2_ = thread_info->pre_state2_->qul_range_;
    int *pre_qul_cnt2_ = (int*)ldm_malloc(pre_qul_range2_ * sizeof(int));
    dma_getn(thread_info->pre_state2_->qul_cnt_, pre_qul_cnt2_, pre_qul_range2_);


    int aft_qul_range2_ = thread_info->aft_state2_->qul_range_;
    int *aft_qul_cnt2_ = (int*)ldm_malloc(aft_qul_range2_ * sizeof(int));
    dma_getn(thread_info->aft_state2_->qul_cnt_, aft_qul_cnt2_, aft_qul_range2_);


    rtc_(&end_tt);
    c_dma_all += end_tt - start_tt;



    rtc_(&start_tt);
    int pre_q20bases1_ = 0;
    int pre_q30bases1_ = 0;
    int pre_tot_bases1_ = 0;
    int pre_gc_bases1_ = 0;
    int pre_real_seq_len1_ = 0;
    int pre_lines1_ = 0;
    int aft_q20bases1_ = 0;
    int aft_q30bases1_ = 0;
    int aft_tot_bases1_ = 0;
    int aft_gc_bases1_ = 0;
    int aft_real_seq_len1_ = 0;
    int aft_lines1_ = 0;


    int pre_q20bases2_ = 0;
    int pre_q30bases2_ = 0;
    int pre_tot_bases2_ = 0;
    int pre_gc_bases2_ = 0;
    int pre_real_seq_len2_ = 0;
    int pre_lines2_ = 0;
    int aft_q20bases2_ = 0;
    int aft_q30bases2_ = 0;
    int aft_tot_bases2_ = 0;
    int aft_gc_bases2_ = 0;
    int aft_real_seq_len2_ = 0;
    int aft_lines2_ = 0;



    for(int id = _PEN * BATCH_SIZE; id < data_num1; id += 64 * BATCH_SIZE) {


        rtc_(&start_t);
        int start_pos = id;
        int end_pos = id + BATCH_SIZE;
        if(end_pos > data_num1) end_pos = data_num1;
        int cal_size1 = (*data1)[end_pos - 1].pqual + (*data1)[end_pos - 1].lqual + 1 - (*data1)[start_pos].pname;
        auto item_s1 = (*data1)[start_pos];
        dma_getn((unsigned char*)(item_s1.base + item_s1.pname), local_data1, cal_size1);

        int cal_size2 = (*data2)[end_pos - 1].pqual + (*data2)[end_pos - 1].lqual + 1 - (*data2)[start_pos].pname;
        auto item_s2 = (*data2)[start_pos];
        dma_getn((unsigned char*)(item_s2.base + item_s2.pname), local_data2, cal_size2);
        rtc_(&end_t);
        c_dma2 += end_t - start_t;



        for(int ids = start_pos; ids < end_pos; ids++) {



            rtc_(&start_t);
            neoReference item1;
            dma_getn(&((*data1)[ids]), &item1, 1);
            neoReference item1p;
            dma_getn(&((*data1)[ids]), &item1p, 1);
            int now_pos1 = item1.pname - item_s1.pname;
            item1.base = &(local_data1[now_pos1]) - item1.pname;

            neoReference item2;
            dma_getn(&((*data2)[ids]), &item2, 1);
            neoReference item2p;
            dma_getn(&((*data2)[ids]), &item2p, 1);
            int now_pos2 = item2.pname - item_s2.pname;
            item2.base = &(local_data2[now_pos2]) - item2.pname;
            rtc_(&end_t);
            c_dma3 += end_t - start_t;


            rtc_(&start_t);
            StateInfo(thread_info->pre_state1_, item1, cmd_info_, pre_pos_cnt1_, pre_pos_qul1_, pre_len_cnt1_, pre_gc_cnt1_, pre_qul_cnt1_, pre_q20bases1_, pre_q30bases1_, pre_tot_bases1_, pre_gc_bases1_, pre_real_seq_len1_, pre_lines1_); 
            StateInfo(thread_info->pre_state2_, item2, cmd_info_, pre_pos_cnt2_, pre_pos_qul2_, pre_len_cnt2_, pre_gc_cnt2_, pre_qul_cnt2_, pre_q20bases2_, pre_q30bases2_, pre_tot_bases2_, pre_gc_bases2_, pre_real_seq_len2_, pre_lines2_); 
            rtc_(&end_t);
            c_state1 += end_t - start_t;






            rtc_(&start_t);
            if (cmd_info_->state_duplicate_) {
                //TODO
            }
            if (cmd_info_->add_umi_) {
                //TODO
            }
            rtc_(&end_t);
            c_dup += end_t - start_t;

            rtc_(&start_t);
            bool trim_res1 = 0;
            bool trim_res2 = 0;
            trim_res1 = TrimSeq(item1, item1p, cmd_info_->trim_front2_, cmd_info_->trim_tail1_, cmd_info_);
            trim_res2 = TrimSeq(item2, item2p, cmd_info_->trim_front2_, cmd_info_->trim_tail1_, cmd_info_);
            if (trim_res1 && trim_res2 && cmd_info_->trim_polyg_) {
                trimPolyGPE(item1, item2, item1p, item2p, cmd_info_->trim_poly_len_);
            }
            if (trim_res1 && trim_res2 && cmd_info_->trim_polyx_) {
                trimPolyXPE(item1, item2, item1p, item2p, cmd_info_->trim_poly_len_);
            }
            //TODO 
            rtc_(&end_t);
            c_trim += end_t - start_t;


            OverlapRes overlap_res;
            if (trim_res1 && trim_res2 && cmd_info_->analyze_overlap_) {
                overlap_res = AnalyzeOverlap(item1, item2, cmd_info_->overlap_diff_limit_, cmd_info_->overlap_require_);
                int now_size = cmd_info_->max_insert_size_;
                if (overlap_res.overlapped) {
                    if (overlap_res.offset > 0)
                        now_size = item1.lseq + item2.lseq - overlap_res.overlap_len;
                    else
                        now_size = overlap_res.overlap_len;
                }
                now_size = min(now_size, cmd_info_->max_insert_size_);
                thread_info->insert_size_dist_[now_size]++;
            }
            if (trim_res1 && trim_res2 && cmd_info_->correct_data_) {
                //Adapter::CorrectData(item1, item2, overlap_res, cmd_info_->isPhred64_);
                //TODO
            }
            if (trim_res1 && trim_res2 && cmd_info_->trim_adapter_) {
                int trimmed= false;
                if (cmd_info_->print_what_trimmed_) {
                    //TODO
                } else {
                    trimmed = TrimAdapterPE(item1, item2, item1p, item2p, overlap_res.offset, overlap_res.overlap_len);
                }
                if (trimmed) {
                    thread_info->aft_state1_->trim_adapter_++;
                    thread_info->aft_state1_->trim_adapter_++;
                    thread_info->aft_state1_->trim_adapter_bases_ +=trimmed;
                }
                int res1 = 0, res2 = 0;
                bool is_trimmed1 = trimmed;
                bool is_trimmed2 = trimmed;

                if (!trimmed) {
                    if (cmd_info_->detect_adapter1_) {
                        int res1;
                        if (cmd_info_->print_what_trimmed_) {
                            //TODO
                        } else {
                            res1 = TrimAdapter(item1, item1p, cmd_info_->adapter_seq1_, false);
                        }
                        if (res1) {
                            is_trimmed1 = true;
                            thread_info->aft_state1_->trim_adapter_++;
                            thread_info->aft_state1_->trim_adapter_bases_ += res1;
                        }
                    }
                    if (cmd_info_->detect_adapter2_) {
                        int res2;
                        if (cmd_info_->print_what_trimmed_) {
                            //TODO
                        } else {
                            res2 = TrimAdapter(item2, item2p, cmd_info_->adapter_seq2_, true);
                        }
                        if (res2) {
                            is_trimmed2 = true;
                            thread_info->aft_state1_->trim_adapter_++;
                            thread_info->aft_state1_->trim_adapter_bases_ += res2;
                        }
                    }
                }

                if(cmd_info_->adapter_from_fasta_.size() > 0) {
                    if (cmd_info_->print_what_trimmed_) {
                        //TODO
                    } else {
                        //TODO
                    }
                    if (res1) {
                        if(!is_trimmed1) thread_info->aft_state1_->trim_adapter_++;
                        thread_info->aft_state1_->trim_adapter_bases_ += res1;
                    }

                    if (cmd_info_->print_what_trimmed_) {
                        //TODO
                    } else {
                        //TODO
                    }
                    if (res2) {
                        if(!is_trimmed2) thread_info->aft_state1_->trim_adapter_++;
                        thread_info->aft_state1_->trim_adapter_bases_ += res2;
                    }
                }

            }

            rtc_(&start_t);
            int filter_res1 = 0;
            int filter_res2 = 0;
            filter_res1 = ReadFiltering(item1, trim_res1, cmd_info_->isPhred64_, cmd_info_);
            filter_res2 = ReadFiltering(item2, trim_res2, cmd_info_->isPhred64_, cmd_info_);
            int filter_res = max(filter_res1, filter_res2);
            rtc_(&end_t);
            c_filter += end_t - start_t;


            if (filter_res == 0) {
                rtc_(&start_t);
                StateInfo(thread_info->aft_state1_, item1, cmd_info_, aft_pos_cnt1_, aft_pos_qul1_, aft_len_cnt1_, aft_gc_cnt1_, aft_qul_cnt1_, aft_q20bases1_, aft_q30bases1_, aft_tot_bases1_, aft_gc_bases1_, aft_real_seq_len1_, aft_lines1_); 
                StateInfo(thread_info->aft_state2_, item2, cmd_info_, aft_pos_cnt2_, aft_pos_qul2_, aft_len_cnt2_, aft_gc_cnt2_, aft_qul_cnt2_, aft_q20bases2_, aft_q30bases2_, aft_tot_bases2_, aft_gc_bases2_, aft_real_seq_len2_, aft_lines2_); 
                rtc_(&end_t);
                c_state2 += end_t - start_t;




                rtc_(&start_t);
                thread_info->aft_state1_->pass_reads_++;
                thread_info->aft_state1_->pass_reads_++;
                if (cmd_info_->write_data_) {

                    UpdateIterm((*pass_data1)[ids], item1p);
                    UpdateIterm((*pass_data2)[ids], item2p);
                }
                rtc_(&end_t);
                c_dma4 += end_t - start_t;




            } else {
                rtc_(&start_t);
                if (cmd_info_->write_data_) {
                    item1p.lname = 0;
                    item1p.lseq = 0;
                    item1p.lstrand = 0;
                    item1p.lqual = 0;
                    UpdateIterm((*pass_data1)[ids], item1p);

                    item2p.lname = 0;
                    item2p.lseq = 0;
                    item2p.lstrand = 0;
                    item2p.lqual = 0;
                    UpdateIterm((*pass_data2)[ids], item2p);
                }
                rtc_(&end_t);
                c_write1 += end_t - start_t;

                rtc_(&start_t);
                if (filter_res == 1) {
                    thread_info->aft_state1_->fail_N_++;
                    thread_info->aft_state1_->fail_N_++;
                } else if (filter_res == 2) {
                    thread_info->aft_state1_->fail_short_++;
                    thread_info->aft_state1_->fail_short_++;
                } else if (filter_res == 3) {
                    thread_info->aft_state1_->fail_long_++;
                    thread_info->aft_state1_->fail_long_++;
                } else if (filter_res == 4) {
                    thread_info->aft_state1_->fail_lowq_++;
                    thread_info->aft_state1_->fail_lowq_++;
                }
                rtc_(&end_t);
                c_write2 += end_t - start_t;
            }

        }
        rtc_(&start_t);
        if(cmd_info_->state_duplicate_) {
            dma_putn(&((*(para->dups))[start_pos]), local_dup_info, BATCH_SIZE);
        }
        rtc_(&end_t);
        c_dma_dup += end_t - start_t;
    }
    rtc_(&end_tt);
    c_loop += end_tt - start_tt;



    rtc_(&start_tt);
    dma_putn(thread_info->pre_state1_->pos_cnt_, pre_pos_cnt1_, pre_pos_seq_len1);
    dma_putn(thread_info->aft_state1_->pos_cnt_, aft_pos_cnt1_, aft_pos_seq_len1);
    ldm_free(pre_pos_cnt1_, pre_pos_seq_len1 * sizeof(int));
    ldm_free(aft_pos_cnt1_, aft_pos_seq_len1 * sizeof(int));

    dma_putn(thread_info->pre_state1_->pos_qul_, pre_pos_qul1_, pre_malloc_seq_len1_);
    dma_putn(thread_info->aft_state1_->pos_qul_, aft_pos_qul1_, aft_malloc_seq_len1_);
    dma_putn(thread_info->pre_state1_->len_cnt_, pre_len_cnt1_, pre_malloc_seq_len1_);
    dma_putn(thread_info->aft_state1_->len_cnt_, aft_len_cnt1_, aft_malloc_seq_len1_);
    dma_putn(thread_info->pre_state1_->gc_cnt_, pre_gc_cnt1_, gcMax + 1);
    dma_putn(thread_info->aft_state1_->gc_cnt_, aft_gc_cnt1_, gcMax + 1);
    dma_putn(thread_info->pre_state1_->qul_cnt_, pre_qul_cnt1_, pre_qul_range1_);
    dma_putn(thread_info->aft_state1_->qul_cnt_, aft_qul_cnt1_, aft_qul_range1_);
    ldm_free(pre_pos_qul1_, pre_malloc_seq_len1_ * sizeof(int));
    ldm_free(aft_pos_qul1_, aft_malloc_seq_len1_ * sizeof(int));
    ldm_free(pre_len_cnt1_, pre_malloc_seq_len1_ * sizeof(int));
    ldm_free(aft_len_cnt1_, aft_malloc_seq_len1_ * sizeof(int));
    ldm_free(pre_gc_cnt1_, (gcMax + 1) * sizeof(int));
    ldm_free(aft_gc_cnt1_, (gcMax + 1) * sizeof(int));
    ldm_free(pre_qul_cnt1_, pre_qul_range1_ * sizeof(int));
    ldm_free(aft_qul_cnt1_, aft_qul_range1_ * sizeof(int));

    thread_info->pre_state1_->q20bases_ += pre_q20bases1_;
    thread_info->pre_state1_->q30bases_ += pre_q30bases1_;
    thread_info->pre_state1_->tot_bases_ += pre_tot_bases1_;
    thread_info->pre_state1_->gc_bases_ += pre_gc_bases1_;
    thread_info->pre_state1_->real_seq_len_ = pre_real_seq_len1_;
    thread_info->pre_state1_->lines_ += pre_lines1_;

    thread_info->aft_state1_->q20bases_ += aft_q20bases1_;
    thread_info->aft_state1_->q30bases_ += aft_q30bases1_;
    thread_info->aft_state1_->tot_bases_ += aft_tot_bases1_;
    thread_info->aft_state1_->gc_bases_ += aft_gc_bases1_;
    thread_info->aft_state1_->real_seq_len_ = aft_real_seq_len1_;
    thread_info->aft_state1_->lines_ += aft_lines1_;




    dma_putn(thread_info->pre_state2_->pos_cnt_, pre_pos_cnt2_, pre_pos_seq_len2);
    dma_putn(thread_info->aft_state2_->pos_cnt_, aft_pos_cnt2_, aft_pos_seq_len2);
    ldm_free(pre_pos_cnt2_, pre_pos_seq_len2 * sizeof(int));
    ldm_free(aft_pos_cnt2_, aft_pos_seq_len2 * sizeof(int));

    dma_putn(thread_info->pre_state2_->pos_qul_, pre_pos_qul2_, pre_malloc_seq_len2_);
    dma_putn(thread_info->aft_state2_->pos_qul_, aft_pos_qul2_, aft_malloc_seq_len2_);
    dma_putn(thread_info->pre_state2_->len_cnt_, pre_len_cnt2_, pre_malloc_seq_len2_);
    dma_putn(thread_info->aft_state2_->len_cnt_, aft_len_cnt2_, aft_malloc_seq_len2_);
    dma_putn(thread_info->pre_state2_->gc_cnt_, pre_gc_cnt2_, gcMax + 1);
    dma_putn(thread_info->aft_state2_->gc_cnt_, aft_gc_cnt2_, gcMax + 1);
    dma_putn(thread_info->pre_state2_->qul_cnt_, pre_qul_cnt2_, pre_qul_range2_);
    dma_putn(thread_info->aft_state2_->qul_cnt_, aft_qul_cnt2_, aft_qul_range2_);
    ldm_free(pre_pos_qul2_, pre_malloc_seq_len2_ * sizeof(int));
    ldm_free(aft_pos_qul2_, aft_malloc_seq_len2_ * sizeof(int));
    ldm_free(pre_len_cnt2_, pre_malloc_seq_len2_ * sizeof(int));
    ldm_free(aft_len_cnt2_, aft_malloc_seq_len2_ * sizeof(int));
    ldm_free(pre_gc_cnt2_, (gcMax + 1) * sizeof(int));
    ldm_free(aft_gc_cnt2_, (gcMax + 1) * sizeof(int));
    ldm_free(pre_qul_cnt2_, pre_qul_range2_ * sizeof(int));
    ldm_free(aft_qul_cnt2_, aft_qul_range2_ * sizeof(int));

    thread_info->pre_state2_->q20bases_ += pre_q20bases2_;
    thread_info->pre_state2_->q30bases_ += pre_q30bases2_;
    thread_info->pre_state2_->tot_bases_ += pre_tot_bases2_;
    thread_info->pre_state2_->gc_bases_ += pre_gc_bases2_;
    thread_info->pre_state2_->real_seq_len_ = pre_real_seq_len2_;
    thread_info->pre_state2_->lines_ += pre_lines2_;

    thread_info->aft_state2_->q20bases_ += aft_q20bases2_;
    thread_info->aft_state2_->q30bases_ += aft_q30bases2_;
    thread_info->aft_state2_->tot_bases_ += aft_tot_bases2_;
    thread_info->aft_state2_->gc_bases_ += aft_gc_bases2_;
    thread_info->aft_state2_->real_seq_len_ = aft_real_seq_len2_;
    thread_info->aft_state2_->lines_ += aft_lines2_;

    rtc_(&end_tt);
    c_dma5 += end_tt - start_tt;
    if(_PEN == 0) {
        printf("== dma all %lf\n", c_dma_all * 1e-6);
        printf("== loop %lf\n", c_loop * 1e-6);
        printf("== dma5 %lf\n", c_dma5 * 1e-6);

        printf("== dma2 %lf\n", c_dma2 * 1e-6);
        printf("== dma3 %lf\n", c_dma3 * 1e-6);
        printf("== state1 %lf\n", c_state1 * 1e-6);
        printf("== dup %lf\n", c_dup * 1e-6);
        printf("== trim %lf\n", c_trim * 1e-6);
        printf("== filter %lf\n", c_filter * 1e-6);
        printf("== state2 %lf\n", c_state2 * 1e-6);
        printf("== dma4 %lf\n", c_dma4 * 1e-6);
        printf("== write1 %lf\n", c_write1 * 1e-6);
        printf("== write2 %lf\n", c_write2 * 1e-6);
        printf("== dma_dup %lf\n", c_dma_dup * 1e-6);
    }
}



ATHREAD_VISIBLE(tgsfunc);
//ATHREAD_VISIBLE(ngsfunc);



