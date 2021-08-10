//
// Created by ylf9811 on 2021/7/7.
//

#include "threadinfo.h"

ThreadInfo::ThreadInfo(CmdInfo *cmd_info) {
    cmd_info_ = cmd_info;
    pre_state1_ = new State(cmd_info->seq_len_, cmd_info->qul_range_);
    pre_state2_ = new State(cmd_info->seq_len_, cmd_info->qul_range_);
    aft_state1_ = new State(cmd_info->seq_len_, cmd_info->qul_range_);
    aft_state2_ = new State(cmd_info->seq_len_, cmd_info->qul_range_);
}

ThreadInfo::~ThreadInfo() {
    delete pre_state1_;
    delete pre_state2_;
    delete aft_state1_;
    delete aft_state2_;
}
