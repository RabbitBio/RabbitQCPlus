//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RABBITQCPLUS_CMDINFO_H
#define RABBITQCPLUS_CMDINFO_H

#include "Globals.h"

class CmdInfo {
public:
    CmdInfo();

public:
    std::string in_file_name1_;
    std::string in_file_name2_;
    std::string out_file_name1_;
    std::string out_file_name2_;
    bool no_adapter_detect_;
    int thread_number_;
    int n_number_limit_;
    int low_qual_perc_limit_;
    int base_len_limit_;
    bool write_data_;
    int64_t out_block_size_;
    int64_t in_file_size1_;
    int64_t in_file_size2_;
    int seq_len_;
    int qul_range_;


};


#endif //RABBITQCPLUS_CMDINFO_H
