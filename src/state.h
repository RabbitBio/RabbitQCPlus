//
// Created by ylf9811 on 2021/7/8.
//

#ifndef RERABBITQC_STATE_H
#define RERABBITQC_STATE_H

#include <cstring>
#include <iostream>
#include "Globals.h"
#include "Reference.h"
#include "cmdinfo.h"

class State {
public:
    State(CmdInfo *cmd_info, int seq_len, int qul_range);

    ~State();

    void StateInfo(neoReference &ref);

    void ExtendBuffer(int old_len, int new_len);

    void Summarize();

    static State *MergeStates(const std::vector<State *> &states);

    static void PrintStates(const State *state);

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


private:
    CmdInfo *cmd_info_;
public:
    const std::map<std::string, int64_t *> &GetHotSeqsDist() const;

public:
    CmdInfo *GetCmdInfo() const;

private:
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
    bool has_summarize_;
    std::map<std::string, int64_t> hot_seqs_info_;
    std::map<std::string, int64_t *> hot_seqs_dist_;
    bool do_over_represent_analyze_;
public:
    const std::map<std::string, int64_t> &GetHotSeqsInfo() const;

private:
    int over_representation_sampling_;


};


#endif //RERABBITQC_STATE_H
