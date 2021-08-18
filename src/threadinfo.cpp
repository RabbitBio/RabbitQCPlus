//
// Created by ylf9811 on 2021/7/7.
//

#include "threadinfo.h"

ThreadInfo::ThreadInfo(CmdInfo *cmd_info, bool is_pe) {
    is_pe_ = is_pe;
    cmd_info_ = cmd_info;
    if (cmd_info->is_TGS_) {
        TGS_state_ = new TGSStats(cmd_info->TGS_min_len_);
    } else {
        if (is_pe) {
            pre_state1_ = new State(cmd_info, cmd_info->seq_len_, cmd_info->qul_range_, false);
            pre_state2_ = new State(cmd_info, cmd_info->seq_len_, cmd_info->qul_range_, true);
            aft_state1_ = new State(cmd_info, cmd_info->seq_len_, cmd_info->qul_range_, false);
            aft_state2_ = new State(cmd_info, cmd_info->seq_len_, cmd_info->qul_range_, true);
        } else {
            pre_state1_ = new State(cmd_info, cmd_info->seq_len_, cmd_info->qul_range_, false);
            aft_state1_ = new State(cmd_info, cmd_info->seq_len_, cmd_info->qul_range_, false);
        }
    }
    insert_size_dist_ = NULL;
    if (cmd_info->no_insert_size_ == 0) {
        insert_size_dist_ = new int64_t[cmd_info->max_insert_size_ + 1];
        memset(insert_size_dist_, 0, sizeof(int64_t) * (cmd_info->max_insert_size_ + 1));
    }
}

ThreadInfo::~ThreadInfo() {
    if (cmd_info_->is_TGS_) {
        delete TGS_state_;
    } else {
        if (is_pe_) {
            delete pre_state1_;
            delete pre_state2_;
            delete aft_state1_;
            delete aft_state2_;
        } else {
            delete pre_state1_;
            delete aft_state1_;
        }

    }
    if (cmd_info_->no_insert_size_ == 0) {
        delete[]insert_size_dist_;
    }

}
