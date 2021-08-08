//
// Created by ylf9811 on 2021/7/8.
//

#ifndef RABBITQCPLUS_STATE_H
#define RABBITQCPLUS_STATE_H

#include <cstring>

#include "Globals.h"
#include "Reference.h"
#include "cmdinfo.h"

class State {
public:
    State(int seq_len, int qul_range);

    ~State();

    void StateInfo(neoReference &ref);

    void Summarize();

    int get_seq_len();

    static State *MergeStates(const std::vector<State *> &states);

    static void PrintStates(const State *state);

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

    bool has_summarize_;


};


#endif //RABBITQCPLUS_STATE_H
