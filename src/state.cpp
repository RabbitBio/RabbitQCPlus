//
// Created by ylf9811 on 2021/7/8.
//

#include "state.h"

#define mod 9999991
#ifdef Vec512

#include <immintrin.h>

#endif
#ifdef Vec256

#include <immintrin.h>

#endif

State::State(CmdInfo *cmd_info, int seq_len, int qul_range, bool is_read2) {
    is_read2_ = is_read2;
    cmd_info_ = cmd_info;
    q20bases_ = 0;
    q30bases_ = 0;
    lines_ = 0;
    malloc_seq_len_ = seq_len;
    qul_range_ = qul_range;
    real_seq_len_ = 0;
    has_summarize_ = false;
    int pos_seq_len = malloc_seq_len_ * 8;
    pos_cnt_ = new int64_t[pos_seq_len];
    memset(pos_cnt_, 0, pos_seq_len * sizeof(int64_t));
    pos_qul_ = new int64_t[pos_seq_len];
    memset(pos_qul_, 0, pos_seq_len * sizeof(int64_t));
    len_cnt_ = new int64_t[malloc_seq_len_];
    memset(len_cnt_, 0, malloc_seq_len_ * sizeof(int64_t));
    gc_cnt_ = new int64_t[101];
    memset(gc_cnt_, 0, 101 * sizeof(int64_t));
    qul_cnt_ = new int64_t[qul_range_];
    memset(qul_cnt_, 0, qul_range_ * sizeof(int64_t));

    kmer_buf_len_ = 2 << (5 * 2);
    kmer_ = new int64_t[kmer_buf_len_];
    memset(kmer_, 0, sizeof(int64_t) * kmer_buf_len_);

    tot_bases_ = 0;
    gc_bases_ = 0;

    do_over_represent_analyze_ = cmd_info->do_overrepresentation_;
    over_representation_sampling_ = cmd_info->overrepresentation_sampling_;

    head_hash_graph_ = NULL;
    hash_graph_ = NULL;
    hash_num_ = 0;

    if (do_over_represent_analyze_) {
        int maxn = 1e7;
        head_hash_graph_ = new int[maxn];
        for (int i = 0; i < maxn; i++)head_hash_graph_[i] = -1;
        if (is_read2_) {
            hash_graph_ = new node[cmd_info->hot_seqs2_.size()];
            for (auto item:cmd_info->hot_seqs2_) {
                HashInsert(item.c_str(), item.length(), cmd_info->eva_len2_);
            }
        } else {
            hash_graph_ = new node[cmd_info->hot_seqs_.size()];
            for (auto item:cmd_info->hot_seqs_) {
                HashInsert(item.c_str(), item.length(), cmd_info->eva_len_);
            }
        }

    }

}


State::~State() {
    delete[] pos_cnt_;
    delete[] pos_qul_;
    delete[] len_cnt_;
    delete[] gc_cnt_;
    delete[] qul_cnt_;
    delete[] kmer_;
    if (do_over_represent_analyze_) {
        for (int i = 0; i < hash_num_; i++)
            delete[] hash_graph_[i].dist;
        delete[] hash_graph_;
        delete[] head_hash_graph_;
    }

}


void State::HashInsert(const char *seq, int len, int eva_len) {
    int64_t now = 0;
    for (int i = 0; i < len; i++) {
        now = now * 5;
        now += valAGCT2[seq[i] & 0x07];
    }
    int ha = (now % mod + mod) % mod;
    hash_graph_[hash_num_].v = now;
    hash_graph_[hash_num_].pre = head_hash_graph_[ha];
    head_hash_graph_[ha] = hash_num_;
    hash_graph_[hash_num_].cnt = 0;
    hash_graph_[hash_num_].seq = std::string(seq, len);
    hash_graph_[hash_num_].dist = new int64_t[eva_len];
    memset(hash_graph_[hash_num_].dist, 0, sizeof(int64_t) * eva_len);
    hash_num_++;
}


void State::HashQueryAndAdd(const char *seq, int offset, int len, int eva_len) {
    int64_t now = 0;
    for (int i = 0; i < len; i++) {
        now = now * 5;
        now += valAGCT2[seq[offset + i] & 0x07];
    }
    int ha = (now % mod + mod) % mod;
    for (int i = head_hash_graph_[ha]; i != -1; i = hash_graph_[i].pre)
        if (hash_graph_[i].v == now) {
            for (int p = offset; p < offset + len && p < eva_len; p++) {
                hash_graph_[i].dist[p]++;
            }
            hash_graph_[i].cnt++;
            return;
        }
}


void State::ExtendBuffer(int old_len, int new_len) {

    malloc_seq_len_ = new_len;
    int64_t *newBuf = NULL;

    newBuf = new int64_t[new_len * 8];
    memset(newBuf, 0, sizeof(int64_t) * new_len * 8);
    memcpy(newBuf, pos_cnt_, sizeof(int64_t) * old_len * 8);
    delete pos_cnt_;
    pos_cnt_ = newBuf;

    newBuf = new int64_t[new_len * 8];
    memset(newBuf, 0, sizeof(int64_t) * new_len * 8);
    memcpy(newBuf, pos_qul_, sizeof(int64_t) * old_len * 8);
    delete pos_qul_;
    pos_qul_ = newBuf;


    newBuf = new int64_t[new_len];
    memset(newBuf, 0, sizeof(int64_t) * new_len);
    memcpy(newBuf, len_cnt_, sizeof(int64_t) * old_len);
    delete len_cnt_;
    len_cnt_ = newBuf;

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

        ExtendBuffer(malloc_seq_len_, std::max(slen + 100, slen * 2));
//        printf("exit because sequence length is too long, malloc_seq_len_ is %d, slen is %d\n", malloc_seq_len_, slen);
//        exit(0);
    }
    real_seq_len_ = std::max(real_seq_len_, slen);
    len_cnt_[slen - 1]++;
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);
    const char q20 = '5';
    const char q30 = '?';
    int flag = 4;
    int gc_cnt = 0;
    int qul_tot = 0;
    int kmer = 0;

#ifdef Vec512
    int i = 0;
    int64_t *p1, *p2, *p3, *p4, *p5, *p6;
    __m512i ad0, ad1, ad2, ad3, ad4, v1, v2, v3, v4, v5, v6, sub33, quamm;
    __m256i bse, and7, add8, idx;
    __m128i ide;
    bse = _mm256_set_epi32(7 * 8, 6 * 8, 5 * 8, 4 * 8, 3 * 8, 2 * 8, 1 * 8, 0 * 8);
    ad1 = _mm512_set1_epi64(1);
    and7 = _mm256_set1_epi32(0x07);
    add8 = _mm256_set1_epi32(64);
    sub33 = _mm512_set1_epi64(33);
    __m512i q20_vec = _mm512_set1_epi64((int64_t) '5');
    __m512i q30_vec = _mm512_set1_epi64((int64_t) '?');
    __m256i con3 = _mm256_set1_epi32(3);
    __m256i con7 = _mm256_set1_epi32(7);


    __m512i qul_tot_vec = _mm512_set1_epi64(0);


    for (; i + 8 <= slen; i += 8) {

        ide = _mm_maskz_loadu_epi8(0xFF, quals + i);
        quamm = _mm512_cvtepi8_epi64(ide);
        ad2 = _mm512_sub_epi64(quamm, sub33);
        __mmask8 q30_mask = _mm512_cmp_epi64_mask(quamm, q30_vec, _MM_CMPINT_NLT);
        __mmask8 q20_mask = _mm512_cmp_epi64_mask(quamm, q20_vec, _MM_CMPINT_NLT);
        ide = _mm_maskz_loadu_epi8(0xFF, bases + i);
        idx = _mm256_cvtepi8_epi32(ide);
        idx = _mm256_and_si256(idx, and7);

        __mmask8 gc_mask1 = _mm256_cmp_epi32_mask(idx, con3, _MM_CMPINT_EQ);
        __mmask8 gc_mask2 = _mm256_cmp_epi32_mask(idx, con7, _MM_CMPINT_EQ);
        gc_cnt += _mm_popcnt_u32(gc_mask1 | gc_mask2);

        idx = _mm256_add_epi32(bse, idx);
        bse = _mm256_add_epi32(bse, add8);

        q20bases_ += _mm_popcnt_u32(q20_mask);
        q30bases_ += _mm_popcnt_u32(q30_mask);

        qul_tot_vec = _mm512_add_epi64(qul_tot_vec, ad2);


        p3 = (int64_t *) pos_cnt_;
        v3 = _mm512_i32gather_epi64(idx, p3, 8);
        v3 = _mm512_add_epi64(v3, ad1);
        _mm512_i32scatter_epi64(p3, idx, v3, 8);

        p4 = (int64_t *) pos_qul_;
        v4 = _mm512_i32gather_epi64(idx, p4, 8);
        v4 = _mm512_add_epi64(v4, ad2);
        _mm512_i32scatter_epi64(p4, idx, v4, 8);

        for (int j = i; j < i + 8; j++) {
            if (bases[j] == 'N')flag = 5;
            int val = valAGCT[bases[j] & 0x07];
            kmer = ((kmer << 2) & 0x3FC) | val;
            if (flag <= 0)kmer_[kmer]++;
            flag--;
        }
    }
    int64_t tmp[8];
    _mm512_store_epi64(tmp, qul_tot_vec);
    for (int ii = 0; ii < 8; ii++)
        qul_tot += tmp[ii];
    for (; i < slen; i++) {
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
        if (bases[i] == 'N')flag = 5;
        int val = valAGCT[bases[i] & 0x07];
        kmer = ((kmer << 2) & 0x3FC) | val;
        if (flag <= 0)kmer_[kmer]++;
        flag--;
    }
//
//    int gc_cnt2 = 0;
//    int qul_tot2 = 0;
//    for (int i = 0; i < slen; i++) {
//        char b = bases[i] & 0x07;
//        if (b == 3 || b == 7)gc_cnt2++;
//        qul_tot2 += quals[i] - 33;
//    }
//    if (gc_cnt2 != gc_cnt) {
//        printf("GG1\n");
//        printf("%d %d\n", gc_cnt2, gc_cnt);
//        exit(0);
//    }
//    if (qul_tot2 != qul_tot) {
//        printf("GG2\n");
//        printf("%d %d\n", qul_tot2, qul_tot);
//        exit(0);
//    }

#elif Vec256
    print("pending...");
#else
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
        if (bases[i] == 'N')flag = 5;
        int val = valAGCT[b];
        kmer = ((kmer << 2) & 0x3FC) | val;
        if (flag <= 0)kmer_[kmer]++;
        flag--;

    }
#endif

    tot_bases_ += slen;
    gc_bases_ += gc_cnt;
    gc_cnt_[int(100.0 * gc_cnt / slen)]++;
    qul_cnt_[int(1.0 * qul_tot / slen)]++;


    // do overrepresentation analysis for 1 of every 20 reads
    if (do_over_represent_analyze_) {
        if (lines_ % over_representation_sampling_ == 0) {
            int eva_len = 0;
            if (is_read2_)eva_len = cmd_info_->eva_len2_;
            else eva_len = cmd_info_->eva_len_;
            const int steps[5] = {10, 20, 40, 100, std::min(150, eva_len - 2)};
            for (int i = 0; i < slen; i++) {
                char *seq = reinterpret_cast<char *>(ref.base + ref.pseq);
                for (int s = 0; s < 5; s++) {
                    if (i + steps[s] < slen)HashQueryAndAdd(seq, i, steps[s], eva_len);
                }
            }

        }
    }
    lines_++;
}


void State::Summarize() {
    if (has_summarize_)return;

    kmer_min_ = kmer_[0];
    kmer_max_ = kmer_[0];
    for (int i = 0; i < kmer_buf_len_; i++) {
//        printf("%d ", kmer_[i]);
        if (kmer_[i] > kmer_max_)
            kmer_max_ = kmer_[i];
        if (kmer_[i] < kmer_min_)
            kmer_min_ = kmer_[i];
    }
//    printf("\n");

    has_summarize_ = true;
}

/**
 * @brief Merge some states to one state
 * @param states
 * @return
 */
State *State::MergeStates(const std::vector<State *> &states) {
    int now_seq_len = 0;
    for (auto item:states) {
        item->Summarize();
        now_seq_len = std::max(now_seq_len, item->real_seq_len_);
    }
    auto *res_state = new State(states[0]->cmd_info_, now_seq_len, states[0]->qul_range_, states[0]->is_read2_);
    res_state->real_seq_len_ = now_seq_len;
    int eva_len = 0;
    if (states[0]->is_read2_)eva_len = states[0]->cmd_info_->eva_len2_;
    else eva_len = states[0]->cmd_info_->eva_len_;

    for (auto item:states) {
        res_state->q20bases_ += item->q20bases_;
        res_state->q30bases_ += item->q30bases_;
        res_state->lines_ += item->lines_;
        res_state->tot_bases_ += item->tot_bases_;
        res_state->gc_bases_ += item->gc_bases_;
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
        for (int i = 0; i < res_state->kmer_buf_len_; i++) {
            res_state->kmer_[i] += item->kmer_[i];
        }

        // merge over rep seq
        int hash_num = res_state->hash_num_;
        for (int i = 0; i < hash_num; i++) {
            res_state->hash_graph_[i].cnt += item->hash_graph_[i].cnt;
            for (int j = 0; j < eva_len; j++) {
                res_state->hash_graph_[i].dist[j] += item->hash_graph_[i].dist[j];
            }
        }
    }


    res_state->Summarize();
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
    printf("kmer max is %lld\n", state->kmer_max_);
    printf("kmer min is %lld\n", state->kmer_min_);

//    int now_seq_len = state->real_seq_len_;
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

int64_t State::GetQ20Bases() const {
    return q20bases_;
}

int64_t State::GetQ30Bases() const {
    return q30bases_;
}

int64_t State::GetLines() const {
    return lines_;
}

int State::GetMallocSeqLen() const {
    return malloc_seq_len_;
}

int State::GetQulRange() const {
    return qul_range_;
}

int State::GetRealSeqLen() const {
    return real_seq_len_;
}

int State::GetKmerBufLen() const {
    return kmer_buf_len_;
}

int64_t *State::GetPosQul() const {
    return pos_qul_;
}

int64_t *State::GetPosCnt() const {
    return pos_cnt_;
}

int64_t *State::GetLenCnt() const {
    return len_cnt_;
}

int64_t *State::GetGcCnt() const {
    return gc_cnt_;
}

int64_t *State::GetQulCnt() const {
    return qul_cnt_;
}

int64_t *State::GetKmer() const {
    return kmer_;
}

int64_t State::GetKmerMin() const {
    return kmer_min_;
}

int64_t State::GetKmerMax() const {
    return kmer_max_;
}

bool State::IsHasSummarize() const {
    return has_summarize_;
}

int64_t State::GetTotBases() const {
    return tot_bases_;
}

int64_t State::GetGcBases() const {
    return gc_bases_;
}

//const std::unordered_map<std::string, int64_t> &State::GetHotSeqsInfo() const {
//    return hot_seqs_info_;
//}

CmdInfo *State::GetCmdInfo() const {
    return cmd_info_;
}

int *State::GetHeadHashGraph() const {
    return head_hash_graph_;
}

node *State::GetHashGraph() const {
    return hash_graph_;
}

int State::GetHashNum() const {
    return hash_num_;
}

std::string State::list2string(int64_t *list, int size) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

std::string State::list2string(double *list, int size) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}
