
#include <cstring>
#include <vector>
#include <unordered_map>
#include <unistd.h>

#include "tuple_spawn.hpp"
#include "slave.h"
#include "Reference.h"
#include "cmdinfo.h"
#include "state.h"
#include "tgsstate.h"
#include "threadinfo.h"

extern "C" void* ldm_malloc(size_t size);
extern "C" void ldm_free(void *addr, size_t size);
struct qc_data {
    ThreadInfo **thread_info_;
    CmdInfo *cmd_info_;
    std::vector <neoReference> *data_;
    std::vector <neoReference> *pass_data_;
    int *cnt;
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
    std::vector <neoReference> *data = para->data_;
    ThreadInfo *thread_info = para->thread_info_[_PEN];
    CmdInfo *cmd_info = para->cmd_info_;
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
    for (int id = _PEN; id < data_num; id += 64) {
        auto ref = (*data)[id];
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
        if (cmd_info->isPhred64_) phredSub = 64;
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
    //printf("=== %d base number %d\n", _PEN, base_num);
}

// template void tmpl_entrance(decltype(test__)*);
ATHREAD_VISIBLE(tgsfunc);
