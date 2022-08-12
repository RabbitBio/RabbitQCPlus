//
// Created by ylf9811 on 2021/7/6.
//


#include "cmdinfo.h"


CmdInfo::CmdInfo() {
    in_file_name1_ = "";
    in_file_name2_ = "";
    out_file_name1_ = "";
    out_file_name2_ = "";
    command_ = "";
    isPhred64_ = false;
    isStdin_ = false;
    isStdout_ = false;


    thread_number_ = 4;
    n_number_limit_ = 5;
    low_qual_perc_limit_ = 40;
    length_required_ = 15;
    length_limit_ = 500;
    trim_5end_ = false;
    trim_3end_ = false;
    cut_window_size_ = 4;
    cut_mean_quality_ = 20;
    trim_front1_ = 0;
    trim_tail1_ = 0;
    trim_front2_ = 0;
    trim_tail2_ = 0;

    write_data_ = false;
    in_file_size1_ = 0;
    in_file_size2_ = 0;
    out_block_size_ = 1 << 21;
    seq_len_ = 200;
    qul_range_ = 50;


    trim_adapter_ = true;
    no_trim_adapter_ = false;
    se_auto_detect_adapter_ = true;
    pe_auto_detect_adapter_ = false;
    correct_data_ = false;
    overlap_diff_limit_ = 5;
    overlap_require_ = 30;
    detect_adapter1_ = false;
    detect_adapter2_ = false;
    adapter_seq1_ = "";
    adapter_seq2_ = "";
    analyze_overlap_ = false;
    adapter_len_lim_ = 0;
    print_what_trimmed_ = false;

    //duplicate
    state_duplicate_ = true;

    //polyx
    trim_polyg_ = false;
    trim_polyx_ = false;
    trim_poly_len_ = 10;

    //umi
    umi_loc_ = UMI_LOC_NONE;
    umi_len_ = 0;
    umi_prefix_ = "";
    umi_skip_ = 0;
    add_umi_ = false;

    //TGS
    is_TGS_ = false;
    TGS_min_len_ = 200;

    //overrepresentation
    do_overrepresentation_ = false;
    overrepresentation_sampling_ = 20;
    eva_len_ = 150;
    eva_len2_ = 150;
    print_ORP_seqs_ = false;

    //insert size
    no_insert_size_ = false;
    max_insert_size_ = 512;

    //gz
    compression_level_ = 2;

    //interleaved
    interleaved_in_ = false;
    interleaved_out_ = false;

    //parallel gz
    use_pugz_ = false;
    use_pigz_ = false;
    pugz_threads_ = 2;
    pigz_threads_ = 2;
}
