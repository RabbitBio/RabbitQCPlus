//
// Created by ylf9811 on 2021/7/8.
//

#include "state.h"

State::State(int seq_len, int qul_range) {
    q20bases_ = 0;
    q30bases_ = 0;
    lines_ = 0;
    malloc_seq_len_ = seq_len;
    qul_range_ = qul_range;
    real_seq_len_ = 0;
    int pos_seq_len = seq_len * 8;
    pos_cnt_ = new int64_t[pos_seq_len];
    memset(pos_cnt_, 0, pos_seq_len * sizeof(int64_t));
    pos_qul_ = new int64_t[pos_seq_len];
    memset(pos_qul_, 0, pos_seq_len * sizeof(int64_t));
    len_cnt_ = new int64_t[seq_len];
    memset(len_cnt_, 0, seq_len * sizeof(int64_t));
    gc_cnt_ = new int64_t[101];
    memset(gc_cnt_, 0, 101 * sizeof(int64_t));
    qul_cnt_ = new int64_t[qul_range];
    memset(qul_cnt_, 0, qul_range * sizeof(int64_t));

}

int State::get_seq_len() {
    return real_seq_len_;
}

/**
 * @brief State reference information
 * @param ref
 */
void State::StateInfo(neoReference &ref) {
    int slen = ref.lseq;
    int qlen = ref.lqual;
    ASSERT(slen == qlen);
    if (slen > malloc_seq_len_) {
        printf("exit because sequence length is too long\n");
        exit(0);
    }
    lines_++;
    real_seq_len_ = std::max(real_seq_len_, slen);
    len_cnt_[slen - 1]++;
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);
    const char q20 = '5';
    const char q30 = '?';
    int gc_cnt = 0;
    int qul_tot = 0;
    for (int i = 0; i < slen; i++) {
        char b = bases[i] & 0x07;
        if (b == 3 || b == 7)gc_cnt++;
        qul_tot += quals[i] - 33;
        if (quals[i] >= q30) {
            q20bases_++;
            q30bases_++;
        } else if (quals[i] >= q20) {
            q20bases_++;
        }
        pos_cnt_[i * 8 + b]++;
        pos_qul_[i * 8 + b] += quals[i] - 33;
    }
    gc_cnt_[int(100.0 * gc_cnt / slen)]++;
    qul_cnt_[int(1.0 * qul_tot / slen)]++;
}

/**
 * @brief Merge some states to one state
 * @param states
 * @return
 */
State *State::MergeStates(const std::vector<State *> &states) {
    int now_seq_len = 0;
    for (auto item:states) {
        now_seq_len = std::max(now_seq_len, item->real_seq_len_);
    }
    auto *res_state = new State(now_seq_len, states[0]->qul_range_);
    res_state->real_seq_len_ = now_seq_len;
    for (auto item:states) {
        res_state->q20bases_ += item->q20bases_;
        res_state->q30bases_ += item->q30bases_;
        res_state->lines_ += item->lines_;
        for (int i = 0; i < now_seq_len; i++) {
            for (int j = 0; j < 8; j++) {
                res_state->pos_cnt_[i * 8 + j] += item->pos_cnt_[i * 8 + j];
                res_state->pos_qul_[i * 8 + j] += item->pos_qul_[i * 8 + j];
            }
            res_state->len_cnt_[i] += item->len_cnt_[i];
        }
        for (int i = 0; i < res_state->qul_range_; i++) {
            res_state->qul_cnt_[i] += item->qul_cnt_[i];
        }
        for (int i = 0; i <= 100; i++) {
            res_state->gc_cnt_[i] += item->gc_cnt_[i];
        }
    }
    return res_state;
}

/**
 * @brief Print state information to screen
 * @param state
 */
void State::PrintStates(const State *state) {
    printf("q20bases %lld\n", state->q20bases_);
    printf("q30bases %lld\n", state->q30bases_);
    printf("lines %lld\n", state->lines_);
    int now_seq_len = state->real_seq_len_;
//    printf("position--quality :\n");
//    for (int i = 0; i < now_seq_len; i++) {
//        int64_t tot_cnt = 0;
//        int64_t tot_qul = 0;
//        for (int j = 0; j < 8; j++) {
//            tot_cnt += state->pos_cnt_[i * 8 + j];
//            tot_qul += state->pos_qul_[i * 8 + j];
//        }
//        printf("pos %d, quality %.5f\n", i, 1.0 * tot_qul / tot_cnt);
//    }
//    printf("mean_quality--ref_number :\n");
//    for (int i = 0; i < state->qul_range_; i++) {
//        printf("quality %d, ref_number %lld\n", i, state->qul_cnt_[i]);
//    }
//    printf("gc%%--ref_number :\n");
//    for (int i = 0; i <= 100; i++) {
//        printf("gc%% %d, ref_number %lld\n", i, state->gc_cnt_[i]);
//    }
//    printf("seq_len--ref_number :\n");
//    for (int i = 0; i < state->real_seq_len_; i++) {
//        printf("seq_len %d, ref_number %lld\n", i + 1, state->len_cnt_[i]);
//    }


}