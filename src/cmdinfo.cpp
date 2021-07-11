//
// Created by ylf9811 on 2021/7/6.
//


#include "cmdinfo.h"


CmdInfo::CmdInfo() {
    in_file_name1_ = "";
    in_file_name2_ = "";
    out_file_name1_ = "";
    out_file_name2_ = "";

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



    //TODO add this function and set value true
    trim_adapter_ = true;
    auto_detect_adapter_ = true;
    correct_data_ = false;
    overlap_diff_limit_ = 5;
    overlap_require_ = 30;
    detect_adapter1_ = false;
    detect_adapter2_ = false;
    adapter_seq1_ = "";
    adapter_seq2_ = "";
}