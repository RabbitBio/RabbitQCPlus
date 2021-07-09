//
// Created by ylf9811 on 2021/7/6.
//


#include "cmdinfo.h"


CmdInfo::CmdInfo() {
    in_file_name1_ = "";
    in_file_name2_ = "";
    out_file_name1_ = "";
    out_file_name2_ = "";
    no_adapter_detect_ = true;
    thread_number_ = 1;
    n_number_limit_ = 5;
    low_qual_perc_limit_ = 40;
    base_len_limit_ = 15;
    write_data_ = false;
    in_file_size1_ = 0;
    in_file_size2_ = 0;
    //TODO change
    out_block_size_ = 1 << 21;
    //TODO change
    seq_len_ = 200;
    qul_range_ = 50;
}