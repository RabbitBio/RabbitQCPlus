//
// Created by ylf9811 on 2021/7/7.
//

#include "threadinfo.h"

ThreadInfo::ThreadInfo(CmdInfo *cmd_info) {
    cmd_info_ = cmd_info;
    if (cmd_info->is_TGS_) {
        TGS_state_ = new TGSStats(cmd_info->TGS_min_len_);
    } else {
        pre_state1_ = new State(cmd_info->seq_len_, cmd_info->qul_range_);
        pre_state2_ = new State(cmd_info->seq_len_, cmd_info->qul_range_);
        aft_state1_ = new State(cmd_info->seq_len_, cmd_info->qul_range_);
        aft_state2_ = new State(cmd_info->seq_len_, cmd_info->qul_range_);
    }

}

ThreadInfo::~ThreadInfo() {
    if (cmd_info_->is_TGS_) {
        delete TGS_state_;
    } else {
        delete pre_state1_;
        delete pre_state2_;
        delete aft_state1_;
        delete aft_state2_;
    }

}
