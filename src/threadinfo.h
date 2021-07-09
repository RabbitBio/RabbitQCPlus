//
// Created by ylf9811 on 2021/7/7.
//

#ifndef RABBITQCPLUS_THREADINFO_H
#define RABBITQCPLUS_THREADINFO_H

#include "Globals.h"
#include "state.h"
#include "cmdinfo.h"

class ThreadInfo {
public:

    ThreadInfo(CmdInfo *cmd_info);

public:
    CmdInfo *cmd_info_;
    State *pre_state1_;
    State *pre_state2_;
    State *aft_state1_;
    State *aft_state2_;


};


#endif //RABBITQCPLUS_THREADINFO_H
