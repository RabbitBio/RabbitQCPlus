//
// Created by ylf9811 on 2021/7/20.
//

#ifndef RABBITQCPLUS_UMIER_H
#define RABBITQCPLUS_UMIER_H

#include <cstring>
#include "cmdinfo.h"
#include "Reference.h"

class Umier {
public:
    Umier(CmdInfo *cmd_info);

    ~Umier();

    std::string firstIndex(std::string name);

    std::string lastIndex(std::string name);

    void ProcessSe(neoReference &ref);

    void ProcessPe(neoReference &r1, neoReference &r2);

    void addUmiToName(neoReference &ref, std::string umi);


private:
    CmdInfo *cmd_info_;
};


#endif //RABBITQCPLUS_UMIER_H
