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
#define gcMax 100
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
    pos_cnt_ = new int[malloc_seq_len_ * 4];
    memset(pos_cnt_, 0, malloc_seq_len_ * 4 * sizeof(int));
    pos_qul_ = new int[malloc_seq_len_];
    memset(pos_qul_, 0, malloc_seq_len_ * sizeof(int));
    len_cnt_ = new int[malloc_seq_len_];
    memset(len_cnt_, 0, malloc_seq_len_ * sizeof(int));
    gc_cnt_ = new int[gcMax + 100];
    memset(gc_cnt_, 0, (gcMax + 1) * sizeof(int));
    qul_cnt_ = new int[qul_range_];
    memset(qul_cnt_, 0, qul_range_ * sizeof(int));

    //kmer_buf_len_ = 2 << (5 * 2);
    //kmer_ = new int64_t[kmer_buf_len_];
    //memset(kmer_, 0, sizeof(int64_t) * kmer_buf_len_);
    //TODO

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
    //delete[] kmer_;
    //TODO
    if (do_over_represent_analyze_) {
        for (int i = 0; i < hash_num_; i++)
            delete[] hash_graph_[i].dist;
        delete[] hash_graph_;
        delete[] head_hash_graph_;
    }
}


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
    printf("print hash table state info ===========================\n");
    for (int i = 1; i < 100; i++) {
        if (cntTotal[i]) printf("%d %d\n", i, cntTotal[i]);
    }
    printf("=======================================================\n");
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
    int *newBuf = NULL;

    newBuf = new int[new_len * 4];
    memset(newBuf, 0, sizeof(int) * new_len * 4);
    memcpy(newBuf, pos_cnt_, sizeof(int) * old_len * 4);
    delete pos_cnt_;
    pos_cnt_ = newBuf;

    newBuf = new int[new_len];
    memset(newBuf, 0, sizeof(int) * new_len);
    memcpy(newBuf, pos_qul_, sizeof(int) * old_len);
    delete pos_qul_;
    pos_qul_ = newBuf;


    newBuf = new int[new_len];
    memset(newBuf, 0, sizeof(int) * new_len);
    memcpy(newBuf, len_cnt_, sizeof(int) * old_len);
    delete len_cnt_;
    len_cnt_ = newBuf;
}
/*
void print128_8(__m128i var){
    char val[16];
    memcpy(val, &var, sizeof(val));
    printf("Numerical: %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c %c\n", 
            val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7], val[8], val[9], val[10], val[11], val[12], val[13], val[14], val[15]);
}
void print128_32(__m128i var){
    int val[4];
    memcpy(val, &var, sizeof(val));
    printf("Numerical: %d %d %d %d\n", 
            val[0], val[1], val[2], val[3]);
}
void print128_64(__m128i var) {
    int64_t v64val[2];
    memcpy(v64val, &var, sizeof(v64val));
    printf("Numerical: %lld %lld\n", v64val[0], v64val[1]);
}


void print256_32(__m256i var){
    int val[8];
    memcpy(val, &var, sizeof(val));
    printf("Numerical: %d %d %d %d %d %d %d %d\n", 
            val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
}
void print256_64(__m256i var) {
    int64_t v64val[4];
    memcpy(v64val, &var, sizeof(v64val));
    printf("Numerical: %lld %lld %lld %lld\n", v64val[0], v64val[1], v64val[2], v64val[3]);
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
        pos_cnt_[i * 4 + valAGCT[b]]++;
        pos_qul_[i] += q;
        if (bases[i] == 'N') flag = 5;
        int val = valAGCT[b];
        kmer = ((kmer << 2) & 0x3FC) | val;
        if (flag <= 0) kmer_[kmer]++;
        flag--;
    }

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

    //kmer_min_ = kmer_[0];
    //kmer_max_ = kmer_[0];
    //for (int i = 0; i < kmer_buf_len_; i++) {
    //    //        printf("%d ", kmer_[i]);
    //    if (kmer_[i] > kmer_max_)
    //        kmer_max_ = kmer_[i];
    //    if (kmer_[i] < kmer_min_)
    //        kmer_min_ = kmer_[i];
    //}
    //    printf("\n");

    has_summarize_ = true;
    //TODO
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
            for (int j = 0; j < 4; j++) {
                res_state->pos_cnt_[i * 4 + j] += item->pos_cnt_[i * 4 + j];
            }
            res_state->pos_qul_[i] += item->pos_qul_[i];
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
        //for (int i = 0; i < res_state->kmer_buf_len_; i++) {
        //    res_state->kmer_[i] += item->kmer_[i];
        //}

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

    printf("filtering result:\n");
    printf("pass filter read number: %lld\n", state->pass_reads_);
    printf("not pass filter due to too short: %lld\n", state->fail_short_);
    printf("not pass filter due to too long: %lld\n", state->fail_long_);
    printf("not pass filter due to too many N: %lld\n", state->fail_N_);
    printf("not pass filter due to low quality: %lld\n", state->fail_lowq_);
    if (state->cmd_info_->print_what_trimmed_)
        printf("trimmed adapter read number: %lld (all trimmed adapter (len >= %d) can be find in *_trimmed_adapters.txt)\n",
                state->trim_adapter_, state->cmd_info_->adapter_len_lim_);
    else
        printf("trimmed adapter read number: %lld \n", state->trim_adapter_);
    printf("trimmed base number due to adapter: %lld\n", state->trim_adapter_bases_);
}

/**
 * @brief Print state information to screen
 * @param state
 */
void State::PrintStates(const State *state) {
    printf("total bases %lld\n", state->tot_bases_);
    printf("q20 bases %lld\n", state->q20bases_);
    printf("q30 bases %lld\n", state->q30bases_);
    printf("read number %lld\n", state->lines_);
    printf("average read length %.3f\n", state->avg_len);
    //printf("kmer max is %lld\n", state->kmer_max_);
    //printf("kmer min is %lld\n", state->kmer_min_);
    //if(state->do_over_represent_analyze_){
    //	printf("orp qcnt %lld\n",state->over_representation_qcnt_);
    //	printf("orp pcnt %lld\n",state->over_representation_pcnt_);
    //}
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

int *State::GetPosQul() const {
    return pos_qul_;
}

int *State::GetPosCnt() const {
    return pos_cnt_;
}

int *State::GetLenCnt() const {
    return len_cnt_;
}

int *State::GetGcCnt() const {
    return gc_cnt_;
}

int *State::GetQulCnt() const {
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
