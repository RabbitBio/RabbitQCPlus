//
// Created by ylf9811 on 2021/7/8.
//

#ifndef RERABBITQC_STATE_H
#define RERABBITQC_STATE_H

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "Globals.h"
#include "Reference.h"
#include "cmdinfo.h"
#include "nthash.hpp"

struct node {
    int pre, cnt;
    uint64_t v;
    int slen;
    std::string seq;
    int64_t *dist;
};

class State {
public:
    State(CmdInfo *cmd_info, int seq_len, int qul_range, bool is_reed2);

    ~State();

    void StateInfo(neoReference &ref);

    void StateORP(neoReference &ref);

    void ExtendBuffer(int old_len, int new_len);

    void Summarize();

    static State *MergeStates(const std::vector<State *> &states);

    static void PrintStates(const State *state);

    static void PrintFilterResults(const State *state);

    static void PrintAdapterToFile(const State *state);

    static std::string list2string(int64_t *list, int size);

    static std::string list2string(double *list, int size);

    void HashInsert(const char *seq, int len, int eva_len);

    void HashState();

    //void HashQueryAndAdd(uint64_t now, int offset, int len, int eva_len);

    bool HashQueryAndAdd(uint64_t now, int offset, int len, int eva_len);

    int *GetHeadHashGraph() const;

    node *GetHashGraph() const;

    int GetHashNum() const;

    int64_t GetQ20Bases() const;

    int64_t GetQ30Bases() const;

    int64_t GetLines() const;

    int GetMallocSeqLen() const;

    int GetQulRange() const;

    int GetRealSeqLen() const;

    int GetKmerBufLen() const;

    int64_t *GetPosQul() const;

    int64_t *GetPosCnt() const;

    int64_t *GetLenCnt() const;

    int64_t *GetGcCnt() const;

    int64_t *GetQulCnt() const;

    int64_t *GetKmer() const;

    int64_t GetKmerMin() const;

    int64_t GetKmerMax() const;

    bool IsHasSummarize() const;

    int64_t GetTotBases() const;

    int64_t GetGcBases() const;

    double GetOrpCost() const;

    double GetAvgLen() const;

    void AddPassReads();

    void AddFailShort();

    void AddFailLong();

    void AddFailN();

    void AddFailLowq();

    void AddTrimAdapter();

    void AddTrimAdapterBase(int cnt);


    unordered_map<string, int> GetAdapterMap();

    CmdInfo *GetCmdInfo() const;

private:
    CmdInfo *cmd_info_;

    int64_t q20bases_;
    int64_t q30bases_;
    int64_t lines_;
    int malloc_seq_len_;
    int qul_range_;
    int real_seq_len_;
    int kmer_buf_len_;
    int64_t *pos_qul_;
    int64_t *pos_cnt_;
    int64_t *len_cnt_;
    int64_t *gc_cnt_;
    int64_t *qul_cnt_;
    int64_t *kmer_;
    int64_t kmer_min_;
    int64_t kmer_max_;
    int64_t tot_bases_;
    int64_t gc_bases_;
    double avg_len;
    bool has_summarize_;
    int *head_hash_graph_;
    node *hash_graph_;
    int hash_num_;
    bool is_read2_;
    bool do_over_represent_analyze_;
    int over_representation_sampling_;
    int64_t over_representation_qcnt_;
    int64_t over_representation_pcnt_;
    int64_t *bf_zone_;
    double orpCost;
public:
    int64_t pass_reads_;
    int64_t fail_short_;
    int64_t fail_long_;
    int64_t fail_N_;
    int64_t fail_lowq_;
    int64_t trim_adapter_;
    int64_t trim_adapter_bases_;
    unordered_map<string, int> adapter_map_;
};


#endif //RERABBITQC_STATE_H
