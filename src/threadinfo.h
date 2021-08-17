//
// Created by ylf9811 on 2021/7/7.
//

#ifndef RERABBITQC_THREADINFO_H
#define RERABBITQC_THREADINFO_H

#include "Globals.h"
#include "state.h"
#include "cmdinfo.h"
#include "tgsstate.h"

class ThreadInfo {
public:

    ThreadInfo(CmdInfo *cmd_info,bool is_pre);

    ~ThreadInfo();

public:
    bool is_pe_;
    CmdInfo *cmd_info_;
    TGSStats *TGS_state_;
    State *pre_state1_;
    State *pre_state2_;
    State *aft_state1_;
    State *aft_state2_;


};


#endif //RERABBITQC_THREADINFO_H
