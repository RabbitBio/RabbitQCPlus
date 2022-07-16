#include "CLI11.hpp"
#include "Globals.h"
#include "cmdinfo.h"
#include "peqc.h"
#include "seqc.h"
#include <iostream>


int main(int argc, char **argv) {

    CmdInfo cmd_info;
    CLI::App app("RabbitQCPlus");
    app.add_option("-i,--inFile1", cmd_info.in_file_name1_, "input fastq name 1");
    app.add_option("-I,--inFile2", cmd_info.in_file_name2_, "input fastq name 2");
    app.add_option("-o,--outFile1", cmd_info.out_file_name1_, "output fastq name 1");
    app.add_option("-O,--outFile2", cmd_info.out_file_name2_, "output fastq name 2");

    app.add_flag("--phred64", cmd_info.isPhred64_, "input is using phred64 scoring, default is phred33");
    app.add_flag("--stdin", cmd_info.isStdin_,
                 "input from stdin, or -i /dev/stdin, only for se data or interleaved pe data(which means use --interleavedIn)");
    app.add_flag("--stdout", cmd_info.isStdout_,
                 "output to stdout, or -o /dev/stdout, only for se data or interleaved pe data(which means use --interleavedOut)");


    app.add_flag("-a,--noTrimAdapter", cmd_info.no_trim_adapter_, "don't trim adapter, default is off");
    app.add_flag("--decAdaForSe", cmd_info.se_auto_detect_adapter_, "detect adapter for se data, default is on");
    app.add_flag("--decAdaForPe", cmd_info.pe_auto_detect_adapter_,
                 "detect adapter for pe data, default is off, tool prefers to use overlap to find adapter");
    app.add_flag("--printWhatTrimmed", cmd_info.print_what_trimmed_,
                 "if print what trimmed to *_trimmed_adapters.txt or not, default is not");
    app.add_option("--adapterSeq1", cmd_info.adapter_seq1_, "specify adapter sequence for read1");
    app.add_option("--adapterSeq2", cmd_info.adapter_seq2_, "specify adapter sequence for read2");
    //app.add_option("--adapterLengthLimit", cmd_info.adapter_len_lim_, "minimum adapter length when trimming, default is 0");

    app.add_flag("-c,--correctData", cmd_info.correct_data_, "correcting low quality bases using information from overlap, default is off");

    app.add_option("-w,--threadNum", cmd_info.thread_number_, "number thread used to do QC, default is 4");

    //filter
    app.add_flag("-5,--trim5End", cmd_info.trim_5end_, "do sliding window from 5end to 3end to trim low quality bases, default is off");
    app.add_flag("-3,--trim3End", cmd_info.trim_3end_, "do sliding window from 5end to 3end to trim low quality bases, default is off");
    app.add_option("--trimFront1", cmd_info.trim_front1_, "number of bases to trim from front in read1, default is 0");
    app.add_option("--trimFront2", cmd_info.trim_front2_, "number of bases to trim from front in read2, default is 0");
    app.add_option("--trimTail1", cmd_info.trim_tail1_, "number of bases to trim from tail in read1, default is 0");
    app.add_option("--trimTail2", cmd_info.trim_tail2_, "number of bases to trim from tail in read2, default is 0");


    app.add_flag("-g,--trimPolyg", cmd_info.trim_polyg_, "do polyg tail trim, default is off");
    app.add_flag("-x,--trimPolyx", cmd_info.trim_polyx_, "do polyx tail trim, default is off");


    app.add_flag("-u,--addUmi", cmd_info.add_umi_, "do unique molecular identifier (umi) processing, default is off");
    app.add_option("--umiLen", cmd_info.umi_len_, "umi length if it is in read1/read2, default is 0");
    string umiLoc = "";
    app.add_option("--umiLoc", umiLoc, "specify the location of umi, can be (index1/index2/read1/read2/per_index/per_read), default is 0");
    app.add_option("--umiPrefix", cmd_info.umi_prefix_, "identification to be added in front of umi, default is no prefix");
    app.add_option("--umiSkip", cmd_info.umi_skip_, "the number bases to skip if umi exists, default is 0");

    //app.add_option("--seqLen", cmd_info.seq_len_, "max sequence length, default is 200");


    //TGS
    //TODO now just se
    app.add_flag("--TGS", cmd_info.is_TGS_, "process third generation sequencing (TGS) data (only for se data), default is off");


    //overrepresentation
    app.add_flag("-p,--doOverrepresentation", cmd_info.do_overrepresentation_,
                 "do over-representation sequence analysis, default is off");
    app.add_option("-P,--overrepresentationSampling", cmd_info.overrepresentation_sampling_,
                   "do overrepresentation every [] reads, default is 1");

    app.add_flag("--printORPSeqs", cmd_info.print_ORP_seqs_,
                 "if print overrepresentation sequences to *ORP_sequences.txt or not, default is not");
    //insert size analyze
    app.add_flag("--noInsertSize", cmd_info.no_insert_size_, "no insert size analysis (only for pe data), default is to do insert size analysis");


    //interleaved
    app.add_flag("--interleavedIn", cmd_info.interleaved_in_,
                 "use interleaved input (only for pe data), default is off");
    app.add_flag("--interleavedOut", cmd_info.interleaved_out_,
                 "use interleaved output (only for pe data), default is off");


    //parallel gz
    app.add_flag("--usePugz", cmd_info.use_pugz_, "use pugz to decompress data, default is off");
    app.add_flag("--usePigz", cmd_info.use_pigz_, "use pigz to compress data, default is off");
    app.add_option("--pugzThread", cmd_info.pugz_threads_, "pugz thread number, default is 2");
    app.add_option("--pigzThread", cmd_info.pigz_threads_, "pigz thread number, default is 2");
    stringstream ss;
    for (int i = 0; i < argc; i++) {
        ss << argv[i] << " ";
    }
    string command = ss.str();

    cmd_info.command_ = command;


    bool quVersion = false;
    app.add_flag("-V,--version", quVersion, "application version");

    CLI11_PARSE(app, argc, argv);

    if (quVersion) {
        printf("0.0.2\n");
        return 0;
    }

    if (ends_with(cmd_info.in_file_name1_, ".gz") == 0) {
        cmd_info.use_pugz_ = false;
    }
    if (ends_with(cmd_info.out_file_name1_, ".gz") == 0 || cmd_info.pigz_threads_ == 1) {
        cmd_info.use_pigz_ = false;
    }
    //    if (cmd_info.use_pigz_) {
    //        string out_name1 = cmd_info.out_file_name1_;
    //        cmd_info.out_file_name1_ = out_name1.substr(0, out_name1.find(".gz"));
    //        string out_name2 = cmd_info.out_file_name2_;
    //        cmd_info.out_file_name2_ = out_name2.substr(0, out_name2.find(".gz"));
    //    }
    if (cmd_info.isStdin_) {
        cmd_info.in_file_name1_ = "/dev/stdin";
    }
    if (cmd_info.in_file_name1_ == "/dev/stdin") {
        cmd_info.se_auto_detect_adapter_ = false;
        cmd_info.pe_auto_detect_adapter_ = false;
        cmd_info.do_overrepresentation_ = false;
    }
    if (!quVersion && cmd_info.in_file_name1_.length() == 0) {
        error_exit("-i/--inFile1 can't be null");
    }
    printf("inFile1 is %s\n", cmd_info.in_file_name1_.c_str());
    if (cmd_info.isStdout_) cmd_info.out_file_name1_ = "/dev/stdout";
    if (cmd_info.in_file_name2_.length()) printf("inFile2 is %s\n", cmd_info.in_file_name2_.c_str());
    if (cmd_info.out_file_name1_.length()) {
        remove(cmd_info.out_file_name1_.c_str());
        printf("outFile1 is %s\n", cmd_info.out_file_name1_.c_str());
    }
    if (cmd_info.out_file_name2_.length()) {
        remove(cmd_info.out_file_name2_.c_str());
        printf("outFile2 is %s\n", cmd_info.out_file_name2_.c_str());
    }


    if (cmd_info.no_trim_adapter_) cmd_info.trim_adapter_ = false;
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
        printf("now doing overrepresentation\n");
        printf("overrepresentation sampling is %d\n", cmd_info.overrepresentation_sampling_);
    }

    if (cmd_info.isPhred64_) {
        printf("now use phred64 input\n");
    }


    if (cmd_info.use_pugz_) {
        printf("now use pugz, pugz thread number is %d\n", cmd_info.pugz_threads_);
    }
    if (cmd_info.use_pigz_) {
        printf("now use pigz, pigz thread number is %d\n", cmd_info.pigz_threads_);
    }
    if (cmd_info.thread_number_ == 1)
        printf("now use %d thread to do QC operations\n", cmd_info.thread_number_);
    else
        printf("now use %d threads to do QC operations\n", cmd_info.thread_number_);
    int mx_len = Adapter::EvalMaxLen(cmd_info.in_file_name1_);
    //printf("auto detect max seqs len is %d\n", mx_len);
    cmd_info.seq_len_ = mx_len;
    double t1 = GetTime();

    if (cmd_info.in_file_name2_.length() || cmd_info.interleaved_in_) {
        //PE
        if ((cmd_info.out_file_name1_.length() > 0 && cmd_info.out_file_name2_.length() > 0) ||
            cmd_info.interleaved_out_) {
            cmd_info.write_data_ = true;
        }

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
                printf("find adapter %s in read1\n", cmd_info.adapter_seq1_.c_str());
                cmd_info.detect_adapter1_ = true;
            } else {
                printf("not find adapter in read1\n");
            }
            if (cmd_info.adapter_seq2_.length()) {
                printf("find adapter %s in read2\n", cmd_info.adapter_seq2_.c_str());
                cmd_info.detect_adapter2_ = true;
            } else {
                printf("not find adapter in read2\n");
            }
#ifdef Verbose
            printf("detect adapter cost %.5f\n", GetTime() - t2);
#endif
        }
#ifdef Verbose
        if (cmd_info.correct_data_) {
            printf("now correct data\n");
        }
#endif
        if (cmd_info.trim_adapter_ || cmd_info.correct_data_ || !cmd_info.no_insert_size_) {
            cmd_info.analyze_overlap_ = true;
            printf("for PE data, overlap analysis is used to find adapter by default\n");
        }

        if (cmd_info.trim_front1_) {
            printf("read1 trim front %d bases\n", cmd_info.trim_front1_);
            cmd_info.trim_front2_ = cmd_info.trim_front1_;
            printf("read2 trim front %d bases\n", cmd_info.trim_front2_);
        }
        if (cmd_info.trim_tail1_) {
            printf("read1 trim tail %d bases\n", cmd_info.trim_tail1_);
            cmd_info.trim_tail2_ = cmd_info.trim_tail1_;
            printf("read2 trim tail %d bases\n", cmd_info.trim_tail2_);
        }

        if (cmd_info.do_overrepresentation_) {
            double t2 = GetTime();
            printf("now doing overrepresent preprocessing part\n");
            Adapter::PreOverAnalyze(cmd_info.in_file_name1_, cmd_info.hot_seqs_, cmd_info.eva_len_);
            Adapter::PreOverAnalyze(cmd_info.in_file_name2_, cmd_info.hot_seqs2_, cmd_info.eva_len2_);
            printf("overrepresent preprocessing part done\n");
            printf("read1 has %d hot sequence\n", cmd_info.hot_seqs_.size());
            printf("read2 has %d hot sequence\n", cmd_info.hot_seqs2_.size());
#ifdef Verbose
            printf("pre over representation cost %.5f\n", GetTime() - t2);
#endif
        }

        if (cmd_info.interleaved_in_) {
            printf("now input use interleaved pe data\n");
        }
        if (cmd_info.interleaved_out_) {
            printf("now output use interleaved pe data\n");
        }
        if (cmd_info.adapter_seq1_.length() > 0)
            cmd_info.adapter_len_lim_ = min(cmd_info.adapter_len_lim_, int(cmd_info.adapter_seq1_.length()));
        if (cmd_info.adapter_seq2_.length() > 0)
            cmd_info.adapter_len_lim_ = min(cmd_info.adapter_len_lim_, int(cmd_info.adapter_seq2_.length()));

        PeQc *pe_qc = new PeQc(&cmd_info);
        pe_qc->ProcessPeFastq();
        delete pe_qc;

    } else {
        //SE
        cmd_info.no_insert_size_ = 1;
        if (cmd_info.out_file_name1_.length() > 0) {
            cmd_info.write_data_ = true;
        }
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
#ifdef Verbose
            printf("detect adapter cost %.5f\n", GetTime() - t2);
#endif
        }
        if (cmd_info.trim_front1_) {
            printf("trim front %d bases\n", cmd_info.trim_front1_);
        }
        if (cmd_info.trim_tail1_) {
            printf("trim tail %d bases\n", cmd_info.trim_tail1_);
        }

        if (cmd_info.do_overrepresentation_) {
            double t2 = GetTime();
            printf("now doing overrepresent preprocessing part\n");
            Adapter::PreOverAnalyze(cmd_info.in_file_name1_, cmd_info.hot_seqs_, cmd_info.eva_len_);
            printf("overrepresent preprocessing part done\n");
            printf("total %d hot sqes\n", cmd_info.hot_seqs_.size());
#ifdef Verbose
            printf("pre over representation cost %.5f\n", GetTime() - t2);
#endif
        }
        if (cmd_info.adapter_seq1_.length() > 0)
            cmd_info.adapter_len_lim_ = min(cmd_info.adapter_len_lim_, int(cmd_info.adapter_seq1_.length()));
        SeQc *se_qc = new SeQc(&cmd_info);
        if (cmd_info.is_TGS_) {
            se_qc->ProcessSeTGS();
        } else {
            se_qc->ProcessSeFastq();
        }
        delete se_qc;
    }
    printf("cmd is %s\n", command.c_str());
    printf("total cost %.5fs\n", GetTime() - t1);
    return 0;
}
