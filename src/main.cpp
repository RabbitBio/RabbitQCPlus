#include <iostream>
#include "CLI11.hpp"
#include "seqc.h"
#include "peqc.h"
#include "cmdinfo.h"
#include "Globals.h"

void error_exit(string log_info) {
    cout << log_info << endl;
    exit(0);
}

int main(int argc, char **argv) {

    CmdInfo cmd_info;
    CLI::App app("RabbitQCPlus");
    auto opt = app.add_option("-i,--inFile1", cmd_info.in_file_name1_, "input fastq name 1, can not be ''");
    opt->required();
    app.add_option("-I,--inFile2", cmd_info.in_file_name2_, "input fastq name 2, can be '' when single data");
    app.add_option("-o,--outFile1", cmd_info.out_file_name1_, "output fastq name 1");
    app.add_option("-O,--outFile2", cmd_info.out_file_name2_, "output fastq name 2");

    app.add_flag("-a,--noTrimAdapter", cmd_info.no_trim_adapter_, "no trim adapter");
    app.add_flag("--decAdaForSe", cmd_info.se_auto_detect_adapter_, "detect adapter for se data");
    app.add_flag("--decAdaForPe", cmd_info.pe_auto_detect_adapter_, "detect adapter for pe data");
    app.add_option("--adapterSeq1", cmd_info.adapter_seq1_, "input adapter sequence1");
    app.add_option("--adapterSeq2", cmd_info.adapter_seq2_, "input adapter sequence2");

    app.add_flag("-c,--correctData", cmd_info.correct_data_, "correct data");

    app.add_option("-w,--threadNum", cmd_info.thread_number_, "number thread used to solve fastq data");

    //filter
    app.add_flag("-5,--trim5End", cmd_info.trim_5end_, "do sliding window 5end trim");
    app.add_flag("-3,--trim3End", cmd_info.trim_3end_, "do sliding window 3end trim");
    app.add_option("--trimFront1", cmd_info.trim_front1_, "ref1 trim front size");
    app.add_option("--trimFront2", cmd_info.trim_front2_, "ref2 trim front size");
    app.add_option("--trimTail1", cmd_info.trim_tail1_, "ref1 trim tail size");
    app.add_option("--trimTail2", cmd_info.trim_tail2_, "ref2 trim tail size");


    app.add_flag("-g,--trimPolyg", cmd_info.trim_polyg_, "do polyg trim");
    app.add_flag("-x,--trimPolyx", cmd_info.trim_polyx_, "do polyx trim");


    app.add_flag("-u,--addUmi", cmd_info.add_umi_, "do umi add");
    app.add_option("--umiLen", cmd_info.umi_len_, "");
    string umiLoc = "";
    app.add_option("--umiLoc", umiLoc, "");
    app.add_option("--umiPrefix", cmd_info.umi_prefix_, "");
    app.add_option("--umiSkip", cmd_info.umi_skip_, "");

    app.add_option("--seqLen", cmd_info.seq_len_, "");


    //TGS
    //TODO now just se
    app.add_flag("--TGS", cmd_info.is_TGS_, "process TGS");


    //overrepresentation
    app.add_flag("-p,--doOverrepresentation", cmd_info.do_overrepresentation_, "do overrepresentation");
    app.add_option("-P,--overrepresentationSampling", cmd_info.overrepresentation_sampling_,
                   "do overrepresentation every [] reads");

    //insert size analyze
    app.add_flag("--noInsertSize", cmd_info.no_insert_size_, "no insert size analyze");


    CLI11_PARSE(app, argc, argv);
    printf("in1 is %s\n", cmd_info.in_file_name1_.c_str());
    if (cmd_info.in_file_name2_.length())printf("in2 is %s\n", cmd_info.in_file_name2_.c_str());
    if (cmd_info.out_file_name1_.length())printf("out1 is %s\n", cmd_info.out_file_name1_.c_str());
    if (cmd_info.out_file_name2_.length())printf("out2 is %s\n", cmd_info.out_file_name2_.c_str());

    if (cmd_info.no_trim_adapter_)cmd_info.trim_adapter_ = false;
    ASSERT(cmd_info.no_trim_adapter_ != cmd_info.trim_adapter_);
    if (!cmd_info.trim_adapter_) {
        cmd_info.se_auto_detect_adapter_ = false;
        cmd_info.pe_auto_detect_adapter_ = false;
        cmd_info.adapter_seq1_ = "";
        cmd_info.adapter_seq2_ = "";
        printf("no adapter trim!\n");
    }

    if (cmd_info.trim_5end_) {
        printf("now do 5end trim\n");
    }
    if (cmd_info.trim_3end_) {
        printf("now do 3end trim\n");
    }
    if (cmd_info.trim_polyg_) {
        printf("now do polyg trim\n");
    }
    if (cmd_info.trim_polyx_) {
        printf("now do polyx trim\n");
    }


    if (cmd_info.add_umi_) {
        printf("now doing umi add\n");
        printf("umi location is %s\n", umiLoc.c_str());
        printf("umi len is %d\n", cmd_info.umi_len_);
        printf("umi skip is %d\n", cmd_info.umi_skip_);
        printf("umi prefix is %s\n", cmd_info.umi_prefix_.c_str());

        if (umiLoc.empty())
            error_exit("You've enabled UMI by (--addUmi), you should specify the UMI location by (--umiLoc)");
        if (umiLoc != "index1" && umiLoc != "index2" && umiLoc != "read1" && umiLoc != "read2" &&
            umiLoc != "per_index" && umiLoc != "per_read") {
            error_exit("UMI location can only be index1/index2/read1/read2/per_index/per_read");
        }
        if (cmd_info.in_file_name2_.length() == 0 && (umiLoc == "index2" || umiLoc == "read2"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the input data is not paired end.");
        if (cmd_info.umi_len_ == 0 && (umiLoc == "read1" || umiLoc == "read2" || umiLoc == "per_read"))
            error_exit(
                    "You specified the UMI location as " + umiLoc + ", but the length is not specified (--umiLen).");
        if (umiLoc == "index1") {
            cmd_info.umi_loc_ = UMI_LOC_INDEX1;
        } else if (umiLoc == "index2") {
            cmd_info.umi_loc_ = UMI_LOC_INDEX2;
        } else if (umiLoc == "read1") {
            cmd_info.umi_loc_ = UMI_LOC_READ1;
        } else if (umiLoc == "read2") {
            cmd_info.umi_loc_ = UMI_LOC_READ2;
        } else if (umiLoc == "per_index") {
            cmd_info.umi_loc_ = UMI_LOC_PER_INDEX;
        } else if (umiLoc == "per_read") {
            cmd_info.umi_loc_ = UMI_LOC_PER_READ;
        }
    }

    if (cmd_info.do_overrepresentation_) {
        printf("now overrepresentation\n");
        printf("overrepresentation sampling is %d\n", cmd_info.overrepresentation_sampling_);
    }

    printf("now use %d thread\n", cmd_info.thread_number_);

    double t1 = GetTime();
    if (cmd_info.in_file_name2_.length()) {
        if (cmd_info.out_file_name1_.length() > 0 && cmd_info.out_file_name2_.length() > 0) {
            cmd_info.write_data_ = true;
            printf("auto set write_data_ 1\n");
        }
        //calculate file size and estimate reads number
        FILE *p_file;
        p_file = fopen(cmd_info.in_file_name1_.c_str(), "r");
        fseek(p_file, 0, SEEK_END);
        int64_t total_size = ftell(p_file);
        cmd_info.in_file_size1_ = total_size;
        printf("in file total size is %lld\n", total_size);
        printf("my evaluate readNum is %lld\n", int64_t(total_size / 200.0));

        //adapter
        if (cmd_info.adapter_seq1_.length()) {
            printf("input adapter1 is %s\n", cmd_info.adapter_seq1_.c_str());
            if (cmd_info.adapter_seq2_.length() == 0) {
                cmd_info.adapter_seq2_ = cmd_info.adapter_seq1_;
            }
            printf("input adapter2 is %s\n", cmd_info.adapter_seq2_.c_str());
            cmd_info.pe_auto_detect_adapter_ = false;
            cmd_info.detect_adapter1_ = true;
            cmd_info.detect_adapter2_ = true;
        }
        if (cmd_info.pe_auto_detect_adapter_) {
            double t2 = GetTime();
            printf("now auto detect adapter\n");
            cmd_info.adapter_seq1_ = Adapter::AutoDetect(cmd_info.in_file_name1_, cmd_info.trim_tail1_);
            cmd_info.adapter_seq2_ = Adapter::AutoDetect(cmd_info.in_file_name2_, cmd_info.trim_tail1_);
            if (cmd_info.adapter_seq1_.length()) {
                printf("find adapter %s\n", cmd_info.adapter_seq1_.c_str());
                cmd_info.detect_adapter1_ = true;
            } else {
                printf("not find adapter\n");
            }
            if (cmd_info.adapter_seq2_.length()) {
                printf("find adapter %s\n", cmd_info.adapter_seq2_.c_str());
                cmd_info.detect_adapter2_ = true;
            } else {
                printf("not find adapter\n");
            }
            printf("detect adapter cost %.5f\n", GetTime() - t2);

        }
        if (cmd_info.correct_data_) {
            printf("now correct data\n");
        }
        if (cmd_info.trim_adapter_ || cmd_info.correct_data_ || !cmd_info.no_insert_size_) {
            cmd_info.analyze_overlap_ = true;
            printf("now do overlap analyze\n");
        }

        if (cmd_info.trim_front1_) {
            printf("ref1 trim front %d bases\n", cmd_info.trim_front1_);
            cmd_info.trim_front2_ = cmd_info.trim_front1_;
            printf("ref2 trim front %d bases\n", cmd_info.trim_front2_);

        }
        if (cmd_info.trim_tail1_) {
            printf("ref1 trim tail %d bases\n", cmd_info.trim_tail1_);
            cmd_info.trim_tail2_ = cmd_info.trim_tail1_;
            printf("ref2 trim tail %d bases\n", cmd_info.trim_tail2_);
        }

        if (cmd_info.do_overrepresentation_) {
            double t2 = GetTime();
            printf("now pre over represent\n");
            Adapter::PreOverAnalyze(cmd_info.in_file_name1_, cmd_info.hot_seqs_, cmd_info.eva_len_);
            Adapter::PreOverAnalyze(cmd_info.in_file_name2_, cmd_info.hot_seqs2_, cmd_info.eva_len2_);
            printf("pre over represent done\n");
            printf("total %d hot sqes1\n", cmd_info.hot_seqs_.size());
            printf("total %d hot sqes2\n", cmd_info.hot_seqs2_.size());
            printf("pre over representation cost %.5f\n", GetTime() - t2);
        }

        double tp = GetTime();
        PeQc *pe_qc = new PeQc(&cmd_info);
        pe_qc->ProcessPeFastq();
        delete pe_qc;
        printf("part %.5f\n", GetTime() - tp);

    } else {
        cmd_info.no_insert_size_ = 1;
        if (cmd_info.out_file_name1_.length() > 0) {
            cmd_info.write_data_ = true;
            printf("auto set write_data_ 1\n");
        }
        //calculate file size and estimate reads number
        FILE *p_file;
        p_file = fopen(cmd_info.in_file_name1_.c_str(), "r");
        fseek(p_file, 0, SEEK_END);
        int64_t total_size = ftell(p_file);
        cmd_info.in_file_size1_ = total_size;
        printf("in file total size is %lld\n", total_size);
        printf("my evaluate readNum is %lld\n", int64_t(total_size / 200.0));

        //adapter
        if (cmd_info.adapter_seq1_.length()) {
            printf("input adapter is %s\n", cmd_info.adapter_seq1_.c_str());
            cmd_info.se_auto_detect_adapter_ = false;
            cmd_info.detect_adapter1_ = true;
        }
        if (cmd_info.se_auto_detect_adapter_) {
            double t2 = GetTime();
            printf("now auto detect adapter\n");
            cmd_info.adapter_seq1_ = Adapter::AutoDetect(cmd_info.in_file_name1_, cmd_info.trim_tail1_);
            if (cmd_info.adapter_seq1_.length()) {
                printf("find adapter %s\n", cmd_info.adapter_seq1_.c_str());
                cmd_info.detect_adapter1_ = true;
            } else {
                printf("not find adapter\n");
            }
            printf("detect adapter cost %.5f\n", GetTime() - t2);

        }
        if (cmd_info.trim_front1_) {
            printf("trim front %d bases\n", cmd_info.trim_front1_);
        }
        if (cmd_info.trim_tail1_) {
            printf("trim tail %d bases\n", cmd_info.trim_tail1_);
        }

        if (cmd_info.do_overrepresentation_) {
            double t2 = GetTime();
            printf("now pre over represent\n");
            Adapter::PreOverAnalyze(cmd_info.in_file_name1_, cmd_info.hot_seqs_, cmd_info.eva_len_);
            printf("pre over represent done\n");
            printf("total %d hot sqes\n", cmd_info.hot_seqs_.size());
            printf("pre over representation cost %.5f\n", GetTime() - t2);
        }

        double tp = GetTime();
        SeQc *se_qc = new SeQc(&cmd_info);
        if (cmd_info.is_TGS_) {
            se_qc->ProcessSeTGS();
        } else {
            se_qc->ProcessSeFastq();
        }
        delete se_qc;
        printf("part %.5f\n", GetTime() - tp);

    }
    printf("total cost %.5f\n", GetTime() - t1);

    return 0;


}
/*
@SRR2496709.17 17 length=100
CCTTCCCCTCAAGCTCAGGGCCAAGCTGTCCGCCAACCTCGGCTCCTCCGGGCAGCCCTCGCCCGGGGTGCGCCCCGGGGCAGGACCCCCAGCCCACGCC
+SRR2496709.17 17 length=100
18+8?=>?>?==>?>?:>?1>,<<?)>><@@8>@>=>@?@8/?>@=?>=;;@A?BA>+>?8>,=;=?=01>;>>=(8->-8=9-2=>==;19>>A?.7>=
@SRR2496709.17 17 length=100
GGAGGGCAGGGGCGGGCCCTGGGCGTGGGCTGGGGGTCCTGCCCCGGGGCGCACCCCGGGCGAGGGCAGCCCGGAGGAGCCGAGGCTTGCGGACAGCTGG
+SRR2496709.17 17 length=100
819190-<>/;/:<//:===219-*<09.>0?>....891:8=@@;?><-6:1-=>############################################
===================================================
ov 21 79 4 2
reference before correct :
@SRR2496709.17 17 length=100
CCTTCCCCTCAAGCTCAGGGCCAAGCTGTCCGCCAACCTCGGCTCCTCCGGGCAGCCCTCGCCCGGGGTGCGCCCCGGGGCAGGACCCCCAGCCCACGCC
+SRR2496709.17 17 length=100
18+8?=>?>?==>?>?:>?1>,<<?)>><@@8>@>=>@?@8/?>@=?>=;;@A?BA>+>?8>,=;=?=01>;>>=(8->-8=9-2=>==;19>>A?.7>=
@SRR2496709.17 17 length=100
GGAGGGCAGGGGCGGGCCCTGGGCGTGGGCTGGGGGTCCTGCCCCGGGGCGCACCCCGGGCGAGGGCTGCCCGGAGGAGCCGAGGCTGGCGGACAGCTGG
+SRR2496709.17 17 length=100
819190-<>/;/:<//:===219-*<09.>0?>....891:8=@@;?><-6:1-=>###########?###################@############
===================================================




@SRR2496709.17 17 length=100
CCTTCCCCTCAAGCTCAGGGCCAAGCTGTCCGCCAACCTCGGCTCCTCCGGGCAGCCCTCGCCCGGGGTGCGCCCCGGGGCAGGACCCCCAGCCCACGCC
+SRR2496709.17 17 length=100
18+8?=>?>?==>?>?:>?1>,<<?)>><@@8>@>=>@?@8/?>@=?>=;;@A?BA>+>?8>,=;=?=01>;>>=(8->-8=9-2=>==;19>>A?.7>=
@SRR2496709.17 17 length=100
GGAGGGCAGGGGCGGGCCCTGGGCGTGGGCTGGGGGTCCTGCCCCGGGGCGCACCCCGGGCGAGGGCAGCCCGGAGGAGCCGAGGCTTGCGGACAGCTGG
+SRR2496709.17 17 length=100
819190-<>/;/:<//:===219-*<09.>0?>....891:8=@@;?><-6:1-=>############################################
===================================================
ov 21 79 4 2
reference after correct :
@SRR2496709.17 17 length=100
CCTTCCCCTCAAGCTCAGGGCCAAGCTGTCCGCCAACCTCGGCTCCTCCGGGCAGCCCTCGCCCGGGGTGCGCCCCGGGGCAGGACCCCCAGCCCACGCC
+SRR2496709.17 17 length=100
18+8?=>?>?==>?>?:>?1>,<<?)>><@@8>@>=>@?@8/?>@=?>=;;@A?BA>+>?8>,=;=?=01>;>>=(8->-8=9-2=>==;19>>A?.7>=
@SRR2496709.17 17 length=100
GGAGGGCAGGGGCGGGCCCTGGGCGTGGGCTGGGGGTCCTGCCCCGGGGCGCACCCCGGGCGAGGGCTGCCCGGAGGAGCCGAGGCTGGCGGACAGCTGG
+SRR2496709.17 17 length=100
819190-<>/;/:<//:===219-*<09.>0?>@...891:8=@@;?><-6:1?=>############################################
===================================================

 */