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

    overWrite_ = false;
    thread_number_ = 1;
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
    adapter_fasta_file_ = "";

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
    compression_level_ = 4;

    //interleaved
    interleaved_in_ = false;
    interleaved_out_ = false;

    //parallel gz
    use_pugz_ = false;
    use_pigz_ = false;
    pugz_threads_ = 2;
    pigz_threads_ = 2;
}

void CmdInfo::CmdInfoCopy(CmdInfo *cmd_info) {
    in_file_name1_ = cmd_info->in_file_name1_;
    in_file_name2_ = cmd_info->in_file_name2_;
    out_file_name1_ = cmd_info->out_file_name1_;
    out_file_name2_ = cmd_info->out_file_name2_;
    command_ = cmd_info->command_;
    isPhred64_ = cmd_info->isPhred64_;
    isStdin_ = cmd_info->isStdin_;
    isStdout_ = cmd_info->isStdout_;

    overWrite_ = cmd_info->overWrite_;
    thread_number_ = cmd_info->thread_number_;
    n_number_limit_ = cmd_info->n_number_limit_;
    low_qual_perc_limit_ = cmd_info->low_qual_perc_limit_;
    length_required_ = cmd_info->length_required_;
    length_limit_ = cmd_info->length_limit_;
    trim_5end_ = cmd_info->trim_5end_;
    trim_3end_ = cmd_info->trim_3end_;
    cut_window_size_ = cmd_info->cut_window_size_;
    cut_mean_quality_ = cmd_info->cut_mean_quality_;
    trim_front1_ = cmd_info->trim_front1_;
    trim_tail1_ = cmd_info->trim_tail1_;
    trim_front2_ = cmd_info->trim_front2_;
    trim_tail2_ = cmd_info->trim_tail2_;

    write_data_ = cmd_info->write_data_;
    in_file_size1_ = cmd_info->in_file_size1_;
    in_file_size2_ = cmd_info->in_file_size2_;
    out_block_size_ = cmd_info->out_block_size_;
    seq_len_ = cmd_info->seq_len_;
    qul_range_ = cmd_info->qul_range_;


    trim_adapter_ = cmd_info->trim_adapter_;
    no_trim_adapter_ = cmd_info->no_trim_adapter_;
    se_auto_detect_adapter_ = cmd_info->se_auto_detect_adapter_;
    pe_auto_detect_adapter_ = cmd_info->pe_auto_detect_adapter_;
    correct_data_ = cmd_info->correct_data_;
    overlap_diff_limit_ = cmd_info->overlap_diff_limit_;
    overlap_require_ = cmd_info->overlap_require_;
    detect_adapter1_ = cmd_info->detect_adapter1_;
    detect_adapter2_ = cmd_info->detect_adapter2_;
    adapter_seq1_ = cmd_info->adapter_seq1_;
    adapter_seq2_ = cmd_info->adapter_seq2_;
    analyze_overlap_ = cmd_info->analyze_overlap_;
    adapter_len_lim_ = cmd_info->adapter_len_lim_;
    print_what_trimmed_ = cmd_info->print_what_trimmed_;
    adapter_fasta_file_ = cmd_info->adapter_fasta_file_;

    //duplicate
    state_duplicate_ = cmd_info->state_duplicate_;

    //polyx
    trim_polyg_ = cmd_info->trim_polyg_;
    trim_polyx_ =cmd_info->trim_polyx_;
    trim_poly_len_ = cmd_info->trim_poly_len_;

    //umi
    umi_loc_ = cmd_info->umi_loc_;
    umi_len_ = cmd_info->umi_len_;
    umi_prefix_ = cmd_info->umi_prefix_;
    umi_skip_ = cmd_info->umi_skip_;
    add_umi_ = cmd_info->add_umi_;

    //TGS
    is_TGS_ = cmd_info->is_TGS_;
    TGS_min_len_ = cmd_info->TGS_min_len_;

    //overrepresentation
    do_overrepresentation_ = cmd_info->do_overrepresentation_;
    overrepresentation_sampling_ = cmd_info->overrepresentation_sampling_;
    eva_len_ = cmd_info->eva_len_;
    eva_len2_ = cmd_info->eva_len2_;
    print_ORP_seqs_ = cmd_info->print_ORP_seqs_;

    //insert size
    no_insert_size_ = cmd_info->no_insert_size_;
    max_insert_size_ = cmd_info->max_insert_size_;

    //gz
    compression_level_ = cmd_info->compression_level_;

    //interleaved
    interleaved_in_ = cmd_info->interleaved_in_;
    interleaved_out_ = cmd_info->interleaved_out_;

    //parallel gz
    use_pugz_ = cmd_info->use_pugz_;
    use_pigz_ = cmd_info->use_pigz_;
    pugz_threads_ = cmd_info->pugz_threads_;
    pigz_threads_ = cmd_info->pigz_threads_;
}
