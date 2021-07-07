//
// Created by ylf9811 on 2021/7/7.
//

#ifndef RABBITQCPLUS_FILTER_H
#define RABBITQCPLUS_FILTER_H

#include "Globals.h"
#include "Reference.h"
#include "cmdinfo.h"

class Filter {
public:
    Filter();

    Filter(CmdInfo *cmd_info);

    int ReadFiltering(neoReference &ref);

private:
    CmdInfo *cmd_info;
};


#endif //RABBITQCPLUS_FILTER_H
