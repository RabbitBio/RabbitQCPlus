//
// Created by ylf9811 on 2021/7/8.
//

#include "state.h"

#ifdef Vec512

#include <immintrin.h>

#endif
#ifdef Vec256

#include <immintrin.h>

#endif

#define MODB 20
#define gcMax 100000
using namespace std;


State::State(CmdInfo *cmd_info, int seq_len, int qul_range, bool is_read2) {
    is_read2_ = is_read2;
    cmd_info_ = cmd_info;
    orpCost = 0;
    q20bases_ = 0;
    q30bases_ = 0;
    lines_ = 0;
    malloc_seq_len_ = seq_len;
    avg_len = 0;
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
    gc_cnt_ = new int64_t[gcMax + 100];
    memset(gc_cnt_, 0, (gcMax + 1) * sizeof(int64_t));
    qul_cnt_ = new int64_t[qul_range_];
    memset(qul_cnt_, 0, qul_range_ * sizeof(int64_t));

    kmer_buf_len_ = 2 << (5 * 2);
    kmer_ = new int64_t[kmer_buf_len_];
    memset(kmer_, 0, sizeof(int64_t) * kmer_buf_len_);

    tot_bases_ = 0;
    gc_bases_ = 0;


    pass_reads_ = 0;
    fail_short_ = 0;
    fail_long_ = 0;
    fail_N_ = 0;
    fail_lowq_ = 0;
    trim_adapter_ = 0;
    trim_adapter_bases_ = 0;

    do_over_represent_analyze_ = cmd_info->do_overrepresentation_;
    over_representation_sampling_ = cmd_info->overrepresentation_sampling_;

    head_hash_graph_ = NULL;
    hash_graph_ = NULL;
    hash_num_ = 0;
    over_representation_qcnt_ = 0;
    over_representation_pcnt_ = 0;
    if (do_over_represent_analyze_) {
        head_hash_graph_ = new int[(1 << MODB) + 100];
        bf_zone_ = new int64_t[(1 << MODB) / 64 + 100];
        for (int i = 0; i < (1 << MODB) + 10; i++) head_hash_graph_[i] = -1;
        for (int i = 0; i < (1 << MODB) / 64 + 10; i++) bf_zone_[i] = 0;
        if (is_read2_) {
            hash_graph_ = new node[cmd_info->hot_seqs2_.size()];
            for (auto item: cmd_info->hot_seqs2_) {
                HashInsert(item.c_str(), item.length(), cmd_info->eva_len2_);
            }
        } else {
            hash_graph_ = new node[cmd_info->hot_seqs_.size()];
            for (auto item: cmd_info->hot_seqs_) {
                HashInsert(item.c_str(), item.length(), cmd_info->eva_len2_);
            }
        }

        //	HashState();
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

#define BB 233ll
#define MM 1000000000000000003ll

void State::HashInsert(const char *seq, int len, int eva_len) {
    uint64_t now = NT64(seq, len);
    int ha = now & ((1 << MODB) - 1);
    int ha_big = now & ((1 << MODB) - 1);
    bf_zone_[ha_big >> 6] |= (1ll << (ha_big & (0x3f)));
    hash_graph_[hash_num_].v = now;
    hash_graph_[hash_num_].pre = head_hash_graph_[ha];
    head_hash_graph_[ha] = hash_num_;
    hash_graph_[hash_num_].cnt = 0;
    hash_graph_[hash_num_].seq = string(seq, len);
    hash_graph_[hash_num_].dist = new int64_t[eva_len];
    memset(hash_graph_[hash_num_].dist, 0, sizeof(int64_t) * eva_len);
    hash_num_++;
}

void State::HashState() {
    int cntTotal[100];
    for (int i = 0; i < 100; i++) cntTotal[i] = 0;
    for (int ha = 0; ha < (1 << MODB); ha++) {
        int cnt = 0;
        for (int i = head_hash_graph_[ha]; i != -1; i = hash_graph_[i].pre) {
            cnt++;
        }
        cntTotal[cnt]++;
    }
    fprintf(stderr, "print hash table state info ===========================\n");
    for (int i = 1; i < 100; i++) {
        if (cntTotal[i]) fprintf(stderr, "%d %d\n", i, cntTotal[i]);
    }
    fprintf(stderr, "=======================================================\n");
}

inline bool State::HashQueryAndAdd(uint64_t now, int offset, int len, int eva_len) {
    over_representation_qcnt_++;
    int ha_big = now & ((1 << MODB) - 1);
    bool pass_bf = bf_zone_[ha_big >> 6] & (1ll << (ha_big & 0x3f));
    if (!pass_bf) return 0;
    else
        over_representation_pcnt_++;
    int ha = now & ((1 << MODB) - 1);
    for (int i = head_hash_graph_[ha]; i != -1; i = hash_graph_[i].pre) {
        //over_representation_pcnt_++;
        if (hash_graph_[i].v == now && hash_graph_[i].seq.length() == len) {
            //over_representation_pcnt_++;
            for (int p = offset; p < offset + len && p < eva_len; p++) {
                hash_graph_[i].dist[p]++;
            }
            hash_graph_[i].cnt++;
            return 1;
        }
    }
    return 0;
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
/*
void print128_8(__m128i var){
    char val[16];
    memcpy(val, &var, sizeof(val));
    fprintf(stderr, "Numerical: %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c\n", 
            val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7], val[8], val[9], val[10], val[11], val[12], val[13], val[14], val[15]);
}
void print128_32(__m128i var){
    int val[4];
    memcpy(val, &var, sizeof(val));
    fprintf(stderr, "Numerical: %d %d %d %d\n", 
            val[0], val[1], val[2], val[3]);
}
void print128_64(__m128i var) {
    int64_t v64val[2];
    memcpy(v64val, &var, sizeof(v64val));
    fprintf(stderr, "Numerical: %lld %lld\n", v64val[0], v64val[1]);
}


void print256_32(__m256i var){
    int val[8];
    memcpy(val, &var, sizeof(val));
    fprintf(stderr, "Numerical: %d %d %d %d %d %d %d %d\n", 
            val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
}
void print256_64(__m256i var) {
    int64_t v64val[4];
    memcpy(v64val, &var, sizeof(v64val));
    fprintf(stderr, "Numerical: %lld %lld %lld %lld\n", v64val[0], v64val[1], v64val[2], v64val[3]);
}

*/
/**
 * @brief State reference information
 * @param ref
 */
void State::StateInfo(neoReference &ref) {
    int slen = ref.lseq;
    int qlen = ref.lqual;
    ASSERT(slen == qlen);
    if (slen > malloc_seq_len_) {
        ExtendBuffer(malloc_seq_len_, max(slen + 100, slen * 2));
    }
    real_seq_len_ = max(real_seq_len_, slen);
    len_cnt_[slen - 1]++;
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);
    int flag = 4;
    int gc_cnt = 0;
    int qul_tot = 0;
    int kmer = 0;
    int phredSub = 33;
    if (cmd_info_->isPhred64_) phredSub = 64;

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
    sub33 = _mm512_set1_epi64(phredSub);
    __m512i q20_vec = _mm512_set1_epi64(20);
    __m512i q30_vec = _mm512_set1_epi64(30);
    __m512i q0_vec = _mm512_set1_epi64(0);
    __m256i con3 = _mm256_set1_epi32(3);
    __m256i con7 = _mm256_set1_epi32(7);


    __m512i qul_tot_vec = _mm512_set1_epi64(0);


    for (; i + 8 <= slen; i += 8) {


        //int q = max(0, quals[i] - phredSub);
        ide = _mm_maskz_loadu_epi8(0xFF, quals + i);
        quamm = _mm512_cvtepi8_epi64(ide);
        ad2 = _mm512_sub_epi64(quamm, sub33);
        __mmask8 mmsk = _mm512_cmp_epi64_mask(ad2, q0_vec, 1);
        _mm512_mask_set1_epi64(ad2, mmsk, 0);


        //if (q >= 30) {
        //    q20bases_++;
        //    q30bases_++;       
        //} else if (q >= 20) {
        //    q20bases_++;        
        //}
        __mmask8 q30_mask = _mm512_cmp_epi64_mask(ad2, q30_vec, _MM_CMPINT_NLT);
        __mmask8 q20_mask = _mm512_cmp_epi64_mask(ad2, q20_vec, _MM_CMPINT_NLT);
        q20bases_ += _mm_popcnt_u32(q20_mask);
        q30bases_ += _mm_popcnt_u32(q30_mask);


        //qul_tot += q;
        qul_tot_vec = _mm512_add_epi64(qul_tot_vec, ad2);


        //char b = bases[i] & 0x07;
        ide = _mm_maskz_loadu_epi8(0xFF, bases + i);
        idx = _mm256_cvtepi8_epi32(ide);
        idx = _mm256_and_si256(idx, and7);


        //if (b == 3 || b == 7) gc_cnt++;
        __mmask8 gc_mask1 = _mm256_cmp_epi32_mask(idx, con3, _MM_CMPINT_EQ);
        __mmask8 gc_mask2 = _mm256_cmp_epi32_mask(idx, con7, _MM_CMPINT_EQ);
        gc_cnt += _mm_popcnt_u32(gc_mask1 | gc_mask2);



        //i * 8 + b
        idx = _mm256_add_epi32(bse, idx);
        //i++
        bse = _mm256_add_epi32(bse, add8);



        //pos_cnt_[i * 8 + b]++
        p3 = (int64_t *) pos_cnt_;
        v3 = _mm512_i32gather_epi64(idx, p3, 8);
        v3 = _mm512_add_epi64(v3, ad1);
        _mm512_i32scatter_epi64(p3, idx, v3, 8);

        //pos_qul_[i * 8 + b] += q
        p4 = (int64_t *) pos_qul_;
        v4 = _mm512_i32gather_epi64(idx, p4, 8);
        v4 = _mm512_add_epi64(v4, ad2);
        _mm512_i32scatter_epi64(p4, idx, v4, 8);

        for (int j = i; j < i + 8; j++) {
            if (bases[j] == 'N') flag = 5;
            int val = valAGCT[bases[j] & 0x07];
            kmer = ((kmer << 2) & 0x3FC) | val;
            if (flag <= 0) kmer_[kmer]++;
            flag--;
        }
    }
    int64_t tmp[8];
    _mm512_store_epi64(tmp, qul_tot_vec);
    for (int ii = 0; ii < 8; ii++)
        qul_tot += tmp[ii];
    for (; i < slen; i++) {
        char b = bases[i] & 0x07;
        if (b == 3 || b == 7) gc_cnt++;
        int q = max(0, quals[i] - phredSub);
        qul_tot += q;
        if (q >= 30) {
            q20bases_++;
            q30bases_++;
        } else if (q >= 20) {
            q20bases_++;
        }
        pos_cnt_[i * 8 + b]++;
        pos_qul_[i * 8 + b] += q;
        if (bases[i] == 'N') flag = 5;
        int val = valAGCT[bases[i] & 0x07];
        kmer = ((kmer << 2) & 0x3FC) | val;
        if (flag <= 0) kmer_[kmer]++;
        flag--;
    }


#elif Vec256
    int i = 0;
    const long long int *p1, *p2, *p3, *p4, *p5, *p6;
    __m256i ad0, ad1, ad2, ad3, ad4, v1, v2, v3, v4, v5, v6, sub33, quamm;
    __m128i bse, and7, add8, idx;
    __m128i ide;
    bse = _mm_set_epi32(3 * 8, 2 * 8, 1 * 8, 0 * 8);
    ad1 = _mm256_set1_epi64x(1);
    and7 = _mm_set1_epi32(0x07);
    add8 = _mm_set1_epi32(32);
    sub33 = _mm256_set1_epi64x(phredSub);
    __m128i ones_128 = _mm_set1_epi32(-1);
    __m256i ones_256 = _mm256_set1_epi64x(-1);
    __m256i q20_vec = _mm256_set1_epi64x(20);
    __m256i q30_vec = _mm256_set1_epi64x(30);
    __m256i q0_vec = _mm256_set1_epi64x(0);
    __m128i con3 = _mm_set1_epi32(3);
    __m128i con7 = _mm_set1_epi32(7);
    __m256i qul_tot_vec = _mm256_set1_epi64x(0);


    int tag = 0;
    for (; i + 4 <= slen; i += 4) {
        //int q = max(0, quals[i] - phredSub);
        ide = _mm_loadu_si128((__m128i const*)(quals + i));
        quamm = _mm256_cvtepi8_epi64(ide);
        ad2 = _mm256_sub_epi64(quamm, sub33);



        //if (q >= 30) {
        //    q20bases_++;
        //    q30bases_++;       
        //} else if (q >= 20) {
        //    q20bases_++;        
        //}
        __m256i res30 = _mm256_xor_si256(_mm256_cmpgt_epi64(q30_vec, ad2), ones_256);
        int q30_mask = _mm256_movemask_epi8(res30);
        __m256i res20 = _mm256_xor_si256(_mm256_cmpgt_epi64(q20_vec, ad2), ones_256);
        int q20_mask = _mm256_movemask_epi8(res20);
        q20bases_ += _mm_popcnt_u32(q20_mask) >> 3;
        q30bases_ += _mm_popcnt_u32(q30_mask) >> 3;


        //qul_tot += q;
        qul_tot_vec = _mm256_add_epi64(qul_tot_vec, ad2);


        //char b = bases[i] & 0x07;
        ide = _mm_loadu_si128((__m128i const*)(bases + i));
        idx = _mm_cvtepi8_epi32(ide);
        idx = _mm_and_si128(idx, and7);


        //if (b == 3 || b == 7) gc_cnt++;
        __m128i res3 = _mm_cmpeq_epi32(idx, con3);
        __m128i res7 = _mm_cmpeq_epi32(idx, con7);
        int gc3_mask = _mm_movemask_epi8(res3);
        int gc7_mask = _mm_movemask_epi8(res7);
        gc_cnt += _mm_popcnt_u32(gc3_mask | gc7_mask) >> 2;



        //i * 8 + b
        idx = _mm_add_epi32(bse, idx);

        //i++
        bse = _mm_add_epi32(bse, add8);

        //pos_cnt_[i * 8 + b]++
        p3 = (const long long int*)pos_cnt_;
        v3 = _mm256_i32gather_epi64(p3, idx, 8);
        v3 = _mm256_add_epi64(v3, ad1);

        //pos_qul_[i * 8 + b] += q
        p4 = (const long long int*)pos_qul_;
        v4 = _mm256_i32gather_epi64(p4, idx, 8);
        v4 = _mm256_add_epi64(v4, ad2);

        int idx_tmp[4];
        memcpy(idx_tmp, &idx, 4 * sizeof(int));
        int64_t data_tmp1[4];
        memcpy(data_tmp1, &v3, 4 * sizeof(int64_t));
        int64_t data_tmp2[4];
        memcpy(data_tmp2, &v4, 4 * sizeof(int64_t));


        for (int j = i, k = 0; j < i + 4; j++, k++) {
            pos_cnt_[idx_tmp[k]] = data_tmp1[k];
            pos_qul_[idx_tmp[k]] = data_tmp2[k];

            if (bases[j] == 'N') flag = 5;
            int val = valAGCT[bases[j] & 0x07];
            kmer = ((kmer << 2) & 0x3FC) | val;
            if (flag <= 0) kmer_[kmer]++;
            flag--;
        }
    }
    long long int tmp[4];
    _mm256_maskstore_epi64(tmp, ones_256, qul_tot_vec);
    for (int ii = 0; ii < 4; ii++)
        qul_tot += tmp[ii];
    for (; i < slen; i++) {
        char b = bases[i] & 0x07;
        if (b == 3 || b == 7) gc_cnt++;
        int q = max(0, quals[i] - phredSub);
        qul_tot += q;
        if (q >= 30) {
            q20bases_++;
            q30bases_++;
        } else if (q >= 20) {
            q20bases_++;
        }
        pos_cnt_[i * 8 + b]++;
        pos_qul_[i * 8 + b] += q;
        if (bases[i] == 'N') flag = 5;
        int val = valAGCT[bases[i] & 0x07];
        kmer = ((kmer << 2) & 0x3FC) | val;
        if (flag <= 0) kmer_[kmer]++;
        flag--;
    }
#else
    for (int i = 0; i < slen; i++) {
        char b = bases[i] & 0x07;
        if (b == 3 || b == 7) gc_cnt++;
        int q = max(0, quals[i] - phredSub);

        qul_tot += q;
        if (q >= 30) {
            q20bases_++;
            q30bases_++;
        } else if (q >= 20) {
            q20bases_++;
        }
        pos_cnt_[i * 8 + b]++;
        pos_qul_[i * 8 + b] += q;
        if (bases[i] == 'N') flag = 5;
        int val = valAGCT[b];
        kmer = ((kmer << 2) & 0x3FC) | val;
        if (flag <= 0) kmer_[kmer]++;
        flag--;
    }
#endif

    tot_bases_ += slen;
    gc_bases_ += gc_cnt;
    gc_cnt_[int(1.0 * gcMax * gc_cnt / slen)]++;
    qul_cnt_[int(1.0 * qul_tot / slen)]++;
    // do overrepresentation analysis for 1 of every 20 reads
    
#ifdef Verbose
    double t0 = GetTime();
#endif
    if (do_over_represent_analyze_) {
        if (lines_ % over_representation_sampling_ == 0) {
            StateORP(ref);
        }
    }
#ifdef Verbose
    orpCost += GetTime() - t0;
#endif
    lines_++;
}

/*
//this function is seems to be moew accurate.
//however change to the next function to ensure that the output is exactly the same as fastp.
void State::StateORP(neoReference &ref){
int slen = ref.lseq;
int eva_len = 0;
if (is_read2_)eva_len = cmd_info_->eva_len2_;
else eva_len = cmd_info_->eva_len_;
int steps[5] = {10, 20, 40, 100, min(150, eva_len - 2)};
char *seq = reinterpret_cast<char *>(ref.base + ref.pseq);
string seqs=string(seq,slen);
for(int s=0;s<5;s++){
if(steps[s]>slen)continue;
unsigned k=steps[s];
uint64_t hVal = NT64(seq, k);
HashQueryAndAdd(hVal, 0, steps[s], eva_len);
for (size_t i = 0; i < slen - k; i++){
hVal = NTF64(hVal,k, seq[i], seq[i+k]);
HashQueryAndAdd(hVal, i+1, steps[s], eva_len);

}

}

}
*/

void State::StateORP(neoReference &ref) {
    int slen = ref.lseq;
    int eva_len = 0;
    if (is_read2_) eva_len = cmd_info_->eva_len2_;
    else
        eva_len = cmd_info_->eva_len_;
    int steps[5] = {10, 20, 40, 100, min(150, eva_len - 2)};
    char *seq = reinterpret_cast<char *>(ref.base + ref.pseq);
    string seqs = string(seq, slen);
    for (int s = 0; s < 5; s++) {
        if (steps[s] > slen) continue;
        int k = steps[s];
        uint64_t hVal = NT64(seq, k);
        int lim = 0;
        bool res = HashQueryAndAdd(hVal, 0, steps[s], eva_len);
        if (res) lim = k;
        for (int i = 0; i < slen - k - 1; i++) {
            hVal = NTF64(hVal, k, seq[i], seq[i + k]);
            if (lim == 0) {
                res = HashQueryAndAdd(hVal, i + 1, steps[s], eva_len);
                if (res) lim = k;
            } else
                lim--;
        }
    }
}

void State::Summarize() {
    if (has_summarize_) return;

    kmer_min_ = kmer_[0];
    kmer_max_ = kmer_[0];
    for (int i = 0; i < kmer_buf_len_; i++) {
        //        fprintf(stderr, "%d ", kmer_[i]);
        if (kmer_[i] > kmer_max_)
            kmer_max_ = kmer_[i];
        if (kmer_[i] < kmer_min_)
            kmer_min_ = kmer_[i];
    }
    //    fprintf(stderr, "\n");

    has_summarize_ = true;
}

/**
 * @brief Merge some states to one state
 * @param states
 * @return
 */
State *State::MergeStates(const vector<State *> &states) {
    int now_seq_len = 0;
    for (auto item: states) {
        item->Summarize();
        now_seq_len = max(now_seq_len, item->real_seq_len_);
    }
    auto *res_state = new State(states[0]->cmd_info_, now_seq_len, states[0]->qul_range_, states[0]->is_read2_);
    res_state->real_seq_len_ = now_seq_len;
    int eva_len = 0;
    if (states[0]->is_read2_) eva_len = states[0]->cmd_info_->eva_len2_;
    else
        eva_len = states[0]->cmd_info_->eva_len_;

    for (auto item: states) {
        res_state->orpCost = max(res_state->orpCost, item->orpCost);
        res_state->q20bases_ += item->q20bases_;
        res_state->q30bases_ += item->q30bases_;
        res_state->lines_ += item->lines_;
        res_state->tot_bases_ += item->tot_bases_;
        res_state->gc_bases_ += item->gc_bases_;
        res_state->pass_reads_ += item->pass_reads_;
        res_state->fail_short_ += item->fail_short_;
        res_state->fail_long_ += item->fail_long_;
        res_state->fail_N_ += item->fail_N_;
        res_state->fail_lowq_ += item->fail_lowq_;
        res_state->trim_adapter_ += item->trim_adapter_;
        res_state->trim_adapter_bases_ += item->trim_adapter_bases_;

        for (auto itmap: item->adapter_map_) {
            string seq = itmap.first;
            int cnt = itmap.second;
            if (res_state->adapter_map_.count(seq)) {
                int cnt_now = res_state->adapter_map_[seq];
                res_state->adapter_map_[seq] = cnt + cnt_now;
            } else {
                res_state->adapter_map_[seq] = cnt;
            }
        }


        for (int i = 0; i < item->real_seq_len_; i++) {
            for (int j = 0; j < 8; j++) {
                res_state->pos_cnt_[i * 8 + j] += item->pos_cnt_[i * 8 + j];
                res_state->pos_qul_[i * 8 + j] += item->pos_qul_[i * 8 + j];
            }
            res_state->len_cnt_[i] += item->len_cnt_[i];
        }
        for (int i = 0; i < res_state->qul_range_; i++) {
            res_state->qul_cnt_[i] += item->qul_cnt_[i];
        }
        //TODO find a better method to get smoother line
        int gcPart = gcMax / 100;
        for (int i = 0; i < 100; i++) {
            double sum = 0;
            for (int j = i * gcPart; j < (i + 1) * gcPart; j++) {
                sum = max(sum, 1.0 * item->gc_cnt_[j]);
            }
            //res_state->gc_cnt_[i] += sum / gcPart;
            res_state->gc_cnt_[i] += sum;
        }
        for (int i = 0; i < res_state->kmer_buf_len_; i++) {
            res_state->kmer_[i] += item->kmer_[i];
        }

        if (res_state->do_over_represent_analyze_) {
            // merge over rep seq
            int hash_num = res_state->hash_num_;
            res_state->over_representation_qcnt_ += item->over_representation_qcnt_;
            res_state->over_representation_pcnt_ += item->over_representation_pcnt_;
            for (int i = 0; i < hash_num; i++) {
                res_state->hash_graph_[i].cnt += item->hash_graph_[i].cnt;
                for (int j = 0; j < eva_len; j++) {
                    res_state->hash_graph_[i].dist[j] += item->hash_graph_[i].dist[j];
                }
            }
        }
    }

    res_state->avg_len = 1.0 * res_state->tot_bases_ / res_state->lines_;
    res_state->Summarize();
    return res_state;
}

void State::PrintAdapterToFile(const State *state) {
    string srr_name;
    bool is_pe = state->cmd_info_->in_file_name2_.length() > 0;
    int64_t tot_trimmed_base = 0;
    ofstream ooff;
    static int cntt = 1;
    if (cntt == 1) {
        srr_name = state->cmd_info_->in_file_name1_;
    } else if (cntt == 2) {
        srr_name = state->cmd_info_->in_file_name2_;
    }
    srr_name = PaseFileName(srr_name);
    string file_name = (is_pe ? "pe_" : "se_") + srr_name + "_trimmed_adapters.txt";
    cntt++;
    ooff.open(file_name);
    ooff << "adapter count\n";
    for (auto item: state->adapter_map_) {
        tot_trimmed_base += item.second;
        //tot_trimmed_base += (item.first.length()) * (item.second);
    }
    for (auto item: state->adapter_map_) {
        if(1.0 * item.second / tot_trimmed_base > 0.01)
            ooff << item.first << " " << item.second << endl;
    }
    ooff.close();
}

void State::PrintFilterResults(const State *state) {

    fprintf(stderr, "filtering result:\n");
    fprintf(stderr, "pass filter read number: %lld\n", state->pass_reads_);
    fprintf(stderr, "not pass filter due to too short: %lld\n", state->fail_short_);
    fprintf(stderr, "not pass filter due to too long: %lld\n", state->fail_long_);
    fprintf(stderr, "not pass filter due to too many N: %lld\n", state->fail_N_);
    fprintf(stderr, "not pass filter due to low quality: %lld\n", state->fail_lowq_);
    if (state->cmd_info_->print_what_trimmed_)
        fprintf(stderr, "trimmed adapter read number: %lld (all trimmed adapter (len >= %d) can be find in *_trimmed_adapters.txt)\n",
                state->trim_adapter_, state->cmd_info_->adapter_len_lim_);
    else
        fprintf(stderr, "trimmed adapter read number: %lld \n", state->trim_adapter_);
    fprintf(stderr, "trimmed base number due to adapter: %lld\n", state->trim_adapter_bases_);
}

/**
 * @brief Print state information to screen
 * @param state
 */
void State::PrintStates(const State *state) {
    fprintf(stderr, "total bases %lld\n", state->tot_bases_);
    fprintf(stderr, "q20 bases %lld\n", state->q20bases_);
    fprintf(stderr, "q30 bases %lld\n", state->q30bases_);
    fprintf(stderr, "read number %lld\n", state->lines_);
    fprintf(stderr, "average read length %.3f\n", state->avg_len);
    //fprintf(stderr, "kmer max is %lld\n", state->kmer_max_);
    //fprintf(stderr, "kmer min is %lld\n", state->kmer_min_);
    //if(state->do_over_represent_analyze_){
    //	fprintf(stderr, "orp qcnt %lld\n",state->over_representation_qcnt_);
    //	fprintf(stderr, "orp pcnt %lld\n",state->over_representation_pcnt_);
    //}
    //    int now_seq_len = state->real_seq_len_;
    //    fprintf(stderr, "position--quality :\n");
    //    for (int i = 0; i < now_seq_len; i++) {
    //        int64_t tot_cnt = 0;
    //        int64_t tot_qul = 0;
    //        for (int j = 0; j < 8; j++) {
    //            tot_cnt += state->pos_cnt_[i * 8 + j];
    //            tot_qul += state->pos_qul_[i * 8 + j];
    //        }
    //        fprintf(stderr, "pos %d, quality %.5f\n", i, 1.0 * tot_qul / tot_cnt);
    //    }
    //    fprintf(stderr, "mean_quality--ref_number :\n");
    //    for (int i = 0; i < state->qul_range_; i++) {
    //        fprintf(stderr, "quality %d, ref_number %lld\n", i, state->qul_cnt_[i]);
    //    }
    //    fprintf(stderr, "gc%%--ref_number :\n");
    //    for (int i = 0; i <= 100; i++) {
    //        fprintf(stderr, "gc%% %d, ref_number %lld\n", i, state->gc_cnt_[i]);
    //    }
    //    fprintf(stderr, "seq_len--ref_number :\n");
    //    for (int i = 0; i < state->real_seq_len_; i++) {
    //        fprintf(stderr, "seq_len %d, ref_number %lld\n", i + 1, state->len_cnt_[i]);
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

//const unordered_map<string, int64_t> &State::GetHotSeqsInfo() const {
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

double State::GetOrpCost() const {
    return orpCost;
}

double State::GetAvgLen() const {
    return avg_len;
}

void State::AddPassReads() {
    pass_reads_++;
}

void State::AddFailShort() {
    fail_short_++;
}


void State::AddFailLong() {
    fail_long_++;
}

void State::AddFailN() {
    fail_N_++;
}

void State::AddFailLowq() {
    fail_lowq_++;
}

void State::AddTrimAdapter() {
    trim_adapter_++;
}

void State::AddTrimAdapterBase(int cnt) {
    trim_adapter_bases_ += cnt;
}

unordered_map<string, int> State::GetAdapterMap() {
    return adapter_map_;
}

string State::list2string(int64_t *list, int size) {
    stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

string State::list2string(double *list, int size) {
    stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}
