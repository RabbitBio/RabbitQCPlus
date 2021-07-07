//
// Created by ylf9811 on 2021/7/7.
//

#ifndef RABBITQCPLUS_THREADINFO_H
#define RABBITQCPLUS_THREADINFO_H

#include "Globals.h"

class ThreadInfo {
public:
    ThreadInfo();

public:
    int64_t q20bases_;
    int64_t q30bases_;
    int64_t lines_;
    int64_t pass_lines_;

};


#endif //RABBITQCPLUS_THREADINFO_H
