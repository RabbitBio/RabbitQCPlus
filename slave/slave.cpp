
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
    //if(_PEN != 0) return;
    std::vector <neoReference> *data = para->data_;
    ThreadInfo *thread_info = para->thread_info_[_PEN];
    CmdInfo *cmd_info = para->cmd_info_;
    int data_num = data->size();
    int base_num = 0;
    for (int id = _PEN; id < data_num; id += 64) {
    //for (int id = 0; id < data_num; id++) {
        para->cnt[_PEN]++;
        auto ref = (*data)[id];
        //base_num += ref.lqual;
        const int rlen = ref.lseq;
        const char *seq = reinterpret_cast<const char *>(ref.base + ref.pseq);
        const char *quality = reinterpret_cast<const char *>(ref.base + ref.pqual);
        int size_range = thread_info->TGS_state_->mMinlen >> 1;
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
            if (qual >= 5) thread_info->TGS_state_->mBases51015Num[0]++;
            if (qual >= 10) thread_info->TGS_state_->mBases51015Num[1]++;
            if (qual >= 15) thread_info->TGS_state_->mBases51015Num[2]++;
        }
        updateTop5(rlen, 1.0 * sumQual / rlen, thread_info->TGS_state_);
        for (int i = 0; i < size_range; ++i) {
            c1 = seq[i];
            thread_info->TGS_state_->head_qual_sum[i] += (quality[i] - phredSub);
            if (c1 == 'A') {
                thread_info->TGS_state_->head_seq_pos_count[0][i]++;
            } else if (c1 == 'T') {
                thread_info->TGS_state_->head_seq_pos_count[1][i]++;
            } else if (c1 == 'C') {
                thread_info->TGS_state_->head_seq_pos_count[2][i]++;
            } else if (c1 == 'G') {
                thread_info->TGS_state_->head_seq_pos_count[3][i]++;
            }
        }
        for (int i = rlen - 1; i >= rlen - size_range; --i) {
            c2 = seq[i];
            thread_info->TGS_state_->tail_qual_sum[i - (rlen - size_range)] += (quality[i] - phredSub);

            if (c2 == 'A') {
                thread_info->TGS_state_->tail_seq_pos_count[0][i - (rlen - size_range)]++;
            } else if (c2 == 'T') {
                thread_info->TGS_state_->tail_seq_pos_count[1][i - (rlen - size_range)]++;
            } else if (c2 == 'C') {
                thread_info->TGS_state_->tail_seq_pos_count[2][i - (rlen - size_range)]++;
            } else if (c2 == 'G') {
                thread_info->TGS_state_->tail_seq_pos_count[3][i - (rlen - size_range)]++;
            }
        }

    }
    //printf("=== %d base number %d\n", _PEN, base_num);
}

// template void tmpl_entrance(decltype(test__)*);
ATHREAD_VISIBLE(tgsfunc);
