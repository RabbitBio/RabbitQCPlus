#include "CLI11.hpp"
#include "Globals.h"
#include "cmdinfo.h"
#include "peqc.h"
#include "seqc.h"
#include "th_ass.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

inline bool exists_file (const std::string& name) {
    if(name == "/dev/stdout") return 0;
    ifstream f(name.c_str());
    return f.good();
}

int main(int argc, char **argv) {
    srand(time(0));
    CmdInfo cmd_info;
    CLI::App app("RabbitQCPlus");
    app.add_option("-i,--inFile1", cmd_info.in_file_name1_, "input fastq name 1");
    app.add_option("-I,--inFile2", cmd_info.in_file_name2_, "input fastq name 2");
    app.add_option("-o,--outFile1", cmd_info.out_file_name1_, "output fastq name 1");
    app.add_option("-O,--outFile2", cmd_info.out_file_name2_, "output fastq name 2");

    app.add_flag("--notKeepOrder", cmd_info.notKeepOrder_, "do not keep order as input when outputting reads (may slightly improve performance if opened), default is off");

    app.add_option("--compressLevel", cmd_info.compression_level_, "output file compression level (1 - 9), default is 4");
    app.add_flag("--useIgzip", cmd_info.use_igzip_, "use igzip instead of pugz/zlib, default is off");

    app.add_flag("--overWrite", cmd_info.overWrite_, "overwrite out file if already exists.");
    app.add_flag("--phred64", cmd_info.isPhred64_, "input is using phred64 scoring, default is phred33");
    app.add_flag("--stdin", cmd_info.isStdin_,
            "input from stdin, or -i /dev/stdin, only for se data or interleaved pe data(which means use --interleavedIn)");
    app.add_flag("--stdout", cmd_info.isStdout_,
            "output to stdout, or -o /dev/stdout, only for se data or interleaved pe data(which means use --interleavedOut)");

    app.add_option("-w,--threadNum", cmd_info.thread_number_, "number of thread used to do QC, including (de)compression for compressed data, default is 8");

    //TGS
    app.add_flag("--TGS", cmd_info.is_TGS_, "process third generation sequencing (TGS) data (only for se data, does not support trimming and will not produce output files), default is off");
    //insert size analyze
    app.add_flag("--noInsertSize", cmd_info.no_insert_size_, "no insert size analysis (only for pe data), default is to do insert size analysis");


    //interleaved
    app.add_flag("--interleavedIn", cmd_info.interleaved_in_,
            "use interleaved input (only for pe data), default is off");
    app.add_flag("--interleavedOut", cmd_info.interleaved_out_,
            "use interleaved output (only for pe data), default is off");




    //trim adapter
    app.add_flag("-a,--noTrimAdapter", cmd_info.no_trim_adapter_, "don't trim adapter, default is off")->group("Trim Adapter");
    app.add_flag("--decAdaForSe", cmd_info.se_auto_detect_adapter_, "detect adapter for se data, default is on")->group("Trim Adapter");
    app.add_flag("--decAdaForPe", cmd_info.pe_auto_detect_adapter_,
            "detect adapter for pe data, default is off, tool prefers to use overlap to find adapter")->group("Trim Adapter");
    app.add_flag("--printWhatTrimmed", cmd_info.print_what_trimmed_,
            "if print what trimmed to *_trimmed_adapters.txt or not, default is not")->group("Trim Adapter");
    app.add_option("--adapterSeq1", cmd_info.adapter_seq1_, "specify adapter sequence for read1")->group("Trim Adapter");
    app.add_option("--adapterSeq2", cmd_info.adapter_seq2_, "specify adapter sequence for read2")->group("Trim Adapter");
    app.add_option("--adapterFastaFile", cmd_info.adapter_fasta_file_, "specify adapter sequences use fasta file")->group("Trim Adapter");
    //app.add_option("--adapterLengthLimit", cmd_info.adapter_len_lim_, "minimum adapter length when trimming, default is 0");




    //filter
    app.add_flag("-5,--trim5End", cmd_info.trim_5end_, "do sliding window from 5end to 3end to trim low quality bases, default is off")->group("Filter");
    app.add_flag("-3,--trim3End", cmd_info.trim_3end_, "do sliding window from 5end to 3end to trim low quality bases, default is off")->group("Filter");
    app.add_option("--trimFront1", cmd_info.trim_front1_, "number of bases to trim from front in read1, default is 0")->group("Filter");
    app.add_option("--trimFront2", cmd_info.trim_front2_, "number of bases to trim from front in read2, default is 0")->group("Filter");
    app.add_option("--trimTail1", cmd_info.trim_tail1_, "number of bases to trim from tail in read1, default is 0")->group("Filter");
    app.add_option("--trimTail2", cmd_info.trim_tail2_, "number of bases to trim from tail in read2, default is 0")->group("Filter");
    app.add_flag("-g,--trimPolyg", cmd_info.trim_polyg_, "do polyg tail trim, default is off")->group("Filter");
    app.add_flag("-x,--trimPolyx", cmd_info.trim_polyx_, "do polyx tail trim, default is off")->group("Filter");


    //umi
    app.add_flag("-u,--addUmi", cmd_info.add_umi_, "do unique molecular identifier (umi) processing, default is off")->group("UMI Operation");
    app.add_option("--umiLen", cmd_info.umi_len_, "umi length if it is in read1/read2, default is 0")->group("UMI Operation");
    string umiLoc = "";
    app.add_option("--umiLoc", umiLoc, "specify the location of umi, can be (index1/index2/read1/read2/per_index/per_read), default is 0")->group("UMI Operation");
    app.add_option("--umiPrefix", cmd_info.umi_prefix_, "identification to be added in front of umi, default is no prefix")->group("UMI Operation");
    app.add_option("--umiSkip", cmd_info.umi_skip_, "the number bases to skip if umi exists, default is 0")->group("UMI Operation");

    //app.add_option("--seqLen", cmd_info.seq_len_, "max sequence length, default is 200");


    //overrepresentation
    app.add_flag("-p,--doOverrepresentation", cmd_info.do_overrepresentation_,
            "do over-representation sequence analysis, default is off")->group("Over-representation Operation");
    app.add_option("-P,--overrepresentationSampling", cmd_info.overrepresentation_sampling_,
            "do overrepresentation every [] reads, default is 20")->group("Over-representation Operation");
    app.add_flag("--printORPSeqs", cmd_info.print_ORP_seqs_,
            "if print overrepresentation sequences to *ORP_sequences.txt or not, default is not")->group("Over-representation Operation");



    app.add_flag("--correctData", cmd_info.correct_data_, "sample correcting low quality bases using information from overlap (faster but less accurate), default is off")->group("Error Correction Operation");

    //care correction engine para
    //
    app.add_flag("--correctWithCare", cmd_info.do_correction_with_care_, "correct data use care engine (slower but much more accurate), default is off")->group("Error Correction Operation");
    app.add_option("--pairmode", cmd_info.pairmode_, "SE (single-end) or PE (paired-end), default is SE")->group("Error Correction Operation");
    app.add_option("--coverage", cmd_info.coverage_, "Estimated coverage of input file (i.e. number_of_reads * read_length / genome_size)")->group("Error Correction Operation");
    app.add_option("--hashmaps", cmd_info.hashmaps_, "The requested number of hash maps. Must be greater than 0. The actual number of used hash maps may be lower to respect the set memory limit, default is 48")->group("Error Correction Operation");
    app.add_option("--kmerlength", cmd_info.kmerlength_, "The kmer length for minhashing. If 0 or missing, it is automatically determined, default is 0")->group("Error Correction Operation");
    app.add_flag("--enforceHashmapCount", cmd_info.enforceHashmapCount_, "If the requested number of hash maps cannot be fulfilled, the program terminates without error correction, default is off")->group("Error Correction Operation");
    app.add_flag("--useQualityScores", cmd_info.useQualityScores_, "If set, quality scores (if any) are considered during read correction, default is off")->group("Error Correction Operation");
    app.add_option("--qualityScoreBits", cmd_info.qualityScoreBits_, "How many bits should be used to store a single quality score. Allowed values: 1,2,8. If not 8, a lossy compression via binning is used, default is 8")->group("Error Correction Operation");
    app.add_flag("--excludeAmbiguous", cmd_info.excludeAmbiguous_, "If set, reads which contain at least one ambiguous nucleotide will not be corrected, default is off")->group("Error Correction Operation");
    app.add_option("--maxmismatchratio", cmd_info.maxmismatchratio_, "Overlap between anchor and candidate must contain at most (maxmismatchratio * overlapsize) mismatches, default is 0.2")->group("Error Correction Operation");
    app.add_option("--minalignmentoverlap", cmd_info.minalignmentoverlap_, "Overlap between anchor and candidate must be at least this long, default is 30")->group("Error Correction Operation");
    app.add_option("--minalignmentoverlapratio", cmd_info.minalignmentoverlapratio_, "Overlap between anchor and candidate must be at least as long as (minalignmentoverlapratio * candidatelength), default is 0.3")->group("Error Correction Operation");
    app.add_option("--errorfactortuning", cmd_info.errorfactortuning_, "errorfactortuning, default is 0.06")->group("Error Correction Operation");
    app.add_option("--coveragefactortuning", cmd_info.coveragefactortuning_, "coveragefactortuning, default is 0.6")->group("Error Correction Operation");
    app.add_flag("--showProgress", cmd_info.showProgress_, "If set, progress bar is shown during correction, default is off")->group("Error Correction Operation");
    app.add_option("--tempdir", cmd_info.tempdir_, "Directory to store temporary files, default is the output directory")->group("Error Correction Operation");
    app.add_option("--save-preprocessedreads-to", cmd_info.save_preprocessedreads_to_, "Save binary dump of data structure which stores input reads to disk")->group("Error Correction Operation");
    app.add_option("--load-preprocessedreads-from", cmd_info.load_preprocessedreads_from_, "Load binary dump of read data structure from disk")->group("Error Correction Operation");
    app.add_option("--save-hashtables-to", cmd_info.save_hashtables_to_, "Save binary dump of hash tables to disk")->group("Error Correction Operation");
    app.add_option("--load-hashtables-from", cmd_info.load_hashtables_from_, "Load binary dump of hash tables from disk")->group("Error Correction Operation");
    app.add_option("--memHashtables", cmd_info.memHashtables_, "Memory limit in bytes for hash tables and hash table construction. Can use suffix K,M,G, e.g. 20G means 20 gigabyte. This option is not a hard limit, default is a bit less than memory")->group("Error Correction Operation");
    app.add_option("--memTotal", cmd_info.memTotal_, "Total memory limit in bytes. Can use suffix K,M,G, e.g. 20G means 20 gigabyte. This option is not a hard limit, default is all free memory")->group("Error Correction Operation");
    app.add_option("--hashloadfactor", cmd_info.hashloadfactor_, "Load factor of hashtables. 0.0 < hashloadfactor < 1.0. Smaller values can improve the runtime at the expense of greater memory usage, default is 0.8")->group("Error Correction Operation");
    app.add_option("--fixedNumberOfReads", cmd_info.fixedNumberOfReads_, "Process only the first n reads, default is 0")->group("Error Correction Operation");
    app.add_flag("--singlehash", cmd_info.singlehash_, "Use 1 hashtables with h smallest unique hashes, default is off")->group("Error Correction Operation");
    //app.add_flag("--gzoutput", cmd_info.gzoutput_, "gz compressed output (very slow), default is off");
    app.add_flag("--correctionQualityLabels", cmd_info.correctionQualityLabels_, "If set, correction quality label will be appended to output read headers, default is off")->group("Error Correction Operation");
    app.add_flag("--candidateCorrection", cmd_info.candidateCorrection_, "If set, candidate reads will be corrected,too, default is off")->group("Error Correction Operation");
    app.add_option("--candidateCorrectionNewColumns", cmd_info.candidateCorrectionNewColumns_, "If candidateCorrection is set, a candidates with an absolute shift of candidateCorrectionNewColumns compared to anchor are corrected, default is 15")->group("Error Correction Operation");
    app.add_option("--correctionType", cmd_info.correctionType_, "0: Classic, 2: Print . Print is only supported in the cpu version, default is 0")->group("Error Correction Operation");
    app.add_option("--correctionTypeCands", cmd_info.correctionTypeCands_, "0: Classic, 2: Print. Print is only supported in the cpu version, default is 0")->group("Error Correction Operation");



    //parallel gz
    //app.add_flag("--usePugz", cmd_info.use_pugz_, "use pugz to decompress data, default is off");
    //app.add_flag("--usePigz", cmd_info.use_pigz_, "use pigz to compress data, default is off");
    app.add_option("--pugzThread", cmd_info.pugz_threads_, "pugz thread number for each input file, automatic assignment according to the number of total threads (-w) by default. Note: must >=2 threads when specified manually\n");
    app.add_option("--pigzThread", cmd_info.pigz_threads_, "pigz thread number for each output file, automatic assignment according to the number of total threads (-w) by default. Note: must >=2 threads when specified manually\n");

    stringstream ss;
    for (int i = 0; i < argc; i++) {
        ss << argv[i] << " ";
    }
    string command = ss.str();

    cmd_info.command_ = command;


    bool quVersion = false;
    app.add_flag("-V,--version", quVersion, "application version");

    CLI11_PARSE(app, argc, argv);


    cmd_info.correct_threadnum_ = cmd_info.thread_number_;

    if(cmd_info.compression_level_ < 1 || cmd_info.compression_level_ > 9){
        fprintf(stderr, "error : compression level should in [1, 9]!\n");
        return 0;
    }

    if (quVersion) {
        fprintf(stderr, "2.3.0\n");
        return 0;
    }

    if(cmd_info.notKeepOrder_) {
        //fprintf(stderr, "now do not print data keeping order.\n");
    }

    if(cmd_info.is_TGS_){
        if(cmd_info.in_file_name2_.length() > 0){
            fprintf(stderr, "WARNING : the TGS module does not support pe inputs, so ignore the -I parameter.\n");
            cmd_info.in_file_name2_ = "";
        }
        if(cmd_info.out_file_name1_.length() != 0 || cmd_info.out_file_name2_.length() != 0){
            fprintf(stderr, "WARNING : the TGS module does not support trimming and will not produce output files, so ignore the -o and -O parameter.\n");
            cmd_info.out_file_name1_ = "";
            cmd_info.out_file_name2_ = "";
        }
    }

    if(cmd_info.in_file_name2_.length() > 0 && cmd_info.out_file_name1_.length() > 0 && cmd_info.out_file_name2_.length() == 0 && cmd_info.interleaved_out_ == 0){
        error_exit("when interleaved_out module is off, two input files can not correspond to one output.\n");
    }


    if(cmd_info.out_file_name2_.length() > 0 && cmd_info.in_file_name2_.length() == 0 && cmd_info.interleaved_in_ == 0){
        error_exit("when interleaved_in module is off, one input file can not correspond to two outputs.\n");
    }


    int is_pe = cmd_info.in_file_name2_.length() != 0;

    //check if files is ok
    int files_ok = 1;
    if (is_pe) {
        if (!cmd_info.interleaved_in_){
            if(ends_with(cmd_info.in_file_name1_, ".gz") != ends_with(cmd_info.in_file_name2_, ".gz"))
                files_ok = 0;
        }
        if (!cmd_info.interleaved_out_){
            if(ends_with(cmd_info.out_file_name1_, ".gz") != ends_with(cmd_info.out_file_name2_, ".gz"))
                files_ok = 0;
        }
    }
    if(files_ok == 0){
        fprintf(stderr, "error : for PE data, both files must be of the same type, i.e. cannot be one compressed and one uncompressed!\n");
        return 0;
    }

    int manual_pugz = 0;
    int manual_pigz = 0;
    if(cmd_info.pugz_threads_ != -1) {
        cmd_info.use_pugz_ = 1;
        manual_pugz = 1;
    }
    if(cmd_info.pigz_threads_ != -1) {
        cmd_info.use_pigz_ = 1;
        manual_pigz = 1;
    }
    int base_mul = 1;
    if(cmd_info.in_file_name2_.length() != 0) base_mul = 2;

    int except_threas = 0;
    if(manual_pugz) except_threas += cmd_info.pugz_threads_ * base_mul;
    if(manual_pigz) except_threas += cmd_info.pigz_threads_ * base_mul;

    if(except_threas >= cmd_info.thread_number_) {
        fprintf(stderr, "error : number of threads used for parallel (de)compression >= total number of threads (-w)!  note: processing PE files may require twice the number of threads to (de)compress the file\n");
        return 0;
    }
    if(manual_pugz > 0 && !ends_with(cmd_info.in_file_name1_, ".gz")) {
        fprintf(stderr, "error : cannot specify pugzThread for uncompressed input file!\n");
        return 0;
    }
    if(manual_pigz > 0 && !ends_with(cmd_info.out_file_name1_, ".gz")) {
        fprintf(stderr, "error : cannot specify pigzThread for uncompressed output file!\n");
        return 0;
    }
    if(manual_pugz == 1 && cmd_info.pugz_threads_ == 1) {
        fprintf(stderr, "error : must be >=2 when using the usePugz parameter.\n");
        return 0;
    }
    if(manual_pigz == 1 && cmd_info.pigz_threads_ == 1) {
        fprintf(stderr, "error : must be >=2 when using the usePigz parameter.\n");
        return 0;
    }
    cmd_info.thread_number_ -= except_threas;


    int ori_threads = cmd_info.thread_number_;

    //TODO not good when use care
    int in_gz = ends_with(cmd_info.in_file_name1_, ".gz");
    int out_gz = ends_with(cmd_info.out_file_name1_, ".gz");
    int in_module = 0;
    if (in_gz && cmd_info.use_igzip_ == 0 && manual_pugz == 0) in_module |= 1;
    if (out_gz && manual_pigz == 0) in_module |= 2;

    //if(cmd_info.in_file_name2_.length() == 0)cmd_info.thread_number_ = min(32, cmd_info.thread_number_);
    //else cmd_info.thread_number_ = min(cmd_info.thread_number_, 64);
    if (is_pe) {
        //PE
        //pugz 16 * 2, w 32, pigz 48 * 2
        if (in_module == 3){
            cmd_info.thread_number_ = min(cmd_info.thread_number_, 160);
        } else if (in_module == 2){
            cmd_info.thread_number_ = min(cmd_info.thread_number_, 128);
        } else if (in_module == 1){
            cmd_info.thread_number_ = min(cmd_info.thread_number_, 64);
        } else {
            cmd_info.thread_number_ = min(cmd_info.thread_number_, 32);
        }
    }else {
        //SE
        //pugz 16, w 20, pigz 48
        if (in_module == 3){
            cmd_info.thread_number_ = min(cmd_info.thread_number_, 84); 
        } else if (in_module == 2){
            cmd_info.thread_number_ = min(cmd_info.thread_number_, 68); 
        } else if (in_module == 1){
            cmd_info.thread_number_ = min(cmd_info.thread_number_, 36); 
        } else {
            cmd_info.thread_number_ = min(cmd_info.thread_number_, 20); 
        }
    }



    int now_threads = cmd_info.thread_number_;

    int t1, t2, t3;
    if (is_pe) {
        //PE
        int tot_threads = cmd_info.thread_number_;
        if (in_module == 3){
            t1 = thread_assignment_pe3[tot_threads].pugz_t;
            t2 = thread_assignment_pe3[tot_threads].w_t;
            t3 = thread_assignment_pe3[tot_threads].pigz_t;
        } else if (in_module == 2){
            t1 = thread_assignment_pe2[tot_threads].pugz_t;
            t2 = thread_assignment_pe2[tot_threads].w_t;
            t3 = thread_assignment_pe2[tot_threads].pigz_t;
        } else if (in_module == 1){
            t1 = thread_assignment_pe1[tot_threads].pugz_t;
            t2 = thread_assignment_pe1[tot_threads].w_t;
            t3 = thread_assignment_pe1[tot_threads].pigz_t;
        } else {
            t1 = thread_assignment_pe0[tot_threads].pugz_t;
            t2 = thread_assignment_pe0[tot_threads].w_t;
            t3 = thread_assignment_pe0[tot_threads].pigz_t;
        }
    }else {
        //SE
        int tot_threads = cmd_info.thread_number_;
        if (in_module == 3){
            t1 = thread_assignment_se3[tot_threads].pugz_t;
            t2 = thread_assignment_se3[tot_threads].w_t;
            t3 = thread_assignment_se3[tot_threads].pigz_t;
        } else if (in_module == 2){
            t1 = thread_assignment_se2[tot_threads].pugz_t;
            t2 = thread_assignment_se2[tot_threads].w_t;
            t3 = thread_assignment_se2[tot_threads].pigz_t;
        } else if (in_module == 1){
            t1 = thread_assignment_se1[tot_threads].pugz_t;
            t2 = thread_assignment_se1[tot_threads].w_t;
            t3 = thread_assignment_se1[tot_threads].pigz_t;
        } else {
            t1 = thread_assignment_se0[tot_threads].pugz_t;
            t2 = thread_assignment_se0[tot_threads].w_t;
            t3 = thread_assignment_se0[tot_threads].pigz_t;
        }
    }


    cmd_info.thread_number_ = t2;
    if (t1 <= 1){
        if(manual_pugz == 0) cmd_info.use_pugz_ = 0;
    } else {
        cmd_info.use_pugz_ = 1;
        cmd_info.pugz_threads_ = t1;
    }
    if (t3 <= 1){
        if(manual_pigz == 0) cmd_info.use_pigz_ = 0;
    } else {
        cmd_info.use_pigz_ = 1;
        cmd_info.pigz_threads_ = t3;
    }

    if(cmd_info.interleaved_in_ || cmd_info.interleaved_out_){
        cmd_info.use_pugz_ = 0;
        cmd_info.use_pigz_ = 0;
    }


    if(ori_threads > now_threads) {
         double threads_rate = 1.0 * ori_threads / now_threads;
         if(cmd_info.use_pugz_ && manual_pugz == 0) cmd_info.pugz_threads_ = ceil(threads_rate * cmd_info.pugz_threads_);
         if(cmd_info.use_pigz_ && manual_pigz == 0) cmd_info.pigz_threads_ = ceil(threads_rate * cmd_info.pigz_threads_);
         cmd_info.thread_number_ = ceil(threads_rate * cmd_info.thread_number_);
    }

#ifdef Verbose
    fprintf(stderr, "t1, t2, t3 is %d %d %d\n", t1, t2, t3);

    fprintf(stderr, " pugz is %d\n", cmd_info.use_pugz_);
    fprintf(stderr, " pigz is %d\n", cmd_info.use_pigz_);
    fprintf(stderr, " pugz threads are %d\n", cmd_info.pugz_threads_);
    fprintf(stderr, " pigz threads are %d\n", cmd_info.pigz_threads_);
    fprintf(stderr, " qc threads are %d\n", cmd_info.thread_number_);
#endif

    if (cmd_info.isStdin_) {
        cmd_info.in_file_name1_ = "/dev/stdin";
    }
    bool do_max_len_eval = true;
    if (cmd_info.in_file_name1_ == "/dev/stdin") {
        cmd_info.se_auto_detect_adapter_ = false;
        cmd_info.pe_auto_detect_adapter_ = false;
        cmd_info.do_overrepresentation_ = false;
        do_max_len_eval = false;
    }
    if (!quVersion && cmd_info.in_file_name1_.length() == 0) {
        error_exit("-i/--inFile1 can't be null");
    }
    fprintf(stderr, "inFile1 is %s\n", cmd_info.in_file_name1_.c_str());
    if (cmd_info.isStdout_) cmd_info.out_file_name1_ = "/dev/stdout";
    if (cmd_info.in_file_name2_.length()) fprintf(stderr, "inFile2 is %s\n", cmd_info.in_file_name2_.c_str());
    if (cmd_info.out_file_name1_.length()) {
        bool res = exists_file(cmd_info.out_file_name1_);
        if(res){
            string tmps;
            if(cmd_info.overWrite_){
                tmps = "y";
            }else{
                char tmp[100];
                fprintf(stderr, "\n");
                fprintf(stderr, "%s already exists, overwrite it or not ? (y/n)\n", cmd_info.out_file_name1_.c_str());
                scanf("%s", tmp);
                tmps = string(tmp);
                while(tmps != "y" && tmps != "n"){
                    fprintf(stderr, "please input y or n\n");
                    scanf("%s", tmp);
                    tmps = string(tmp);
                }
            }
            if(tmps == "y"){
                remove(cmd_info.out_file_name1_.c_str());
            } else if(tmps == "n"){
                return 0;
            } else{
                assert(0);
            }
        }
        fprintf(stderr, "outFile1 is %s\n", cmd_info.out_file_name1_.c_str());
    }
    if (cmd_info.out_file_name2_.length()) {
        bool res = exists_file(cmd_info.out_file_name2_);
        if(res){
            string tmps;
            if(cmd_info.overWrite_){
                tmps = "y";
            }else{
                char tmp[100];
                fprintf(stderr, "\n");
                fprintf(stderr, "%s already exists, overwrite it or not ? (y/n)\n", cmd_info.out_file_name2_.c_str());
                scanf("%s", tmp);
                tmps = string(tmp);
                while(tmps != "y" && tmps != "n"){
                    fprintf(stderr, "please input y or n\n");
                    scanf("%s", tmp);
                    tmps = string(tmp);
                }
            }
            if(tmps == "y"){
                remove(cmd_info.out_file_name2_.c_str());
            } else if(tmps == "n"){
                return 0;
            } else{
                assert(0);
            }
        }
        fprintf(stderr, "outFile2 is %s\n", cmd_info.out_file_name2_.c_str());
    }

    string prefix_name1;
    string prefix_name2;
    string suffix_name1;
    string suffix_name2;

    if(out_gz){
        //fprintf(stderr, "===out1 is %s\n", cmd_info.out_file_name1_.c_str());
        if(cmd_info.out_file_name2_.length() > 0){
            //fprintf(stderr, "===out2 is %s\n", cmd_info.out_file_name2_.c_str());
        }
        //fprintf(stderr, "now change out file name...\n");
        int real_begin1 = 0;
        for(int i = cmd_info.out_file_name1_.length() - 1; i >= 0; i--){
            if(cmd_info.out_file_name1_[i] == '/'){
                real_begin1 = i + 1;
                break;
            }
        }
        prefix_name1 = cmd_info.out_file_name1_.substr(0, real_begin1);
        suffix_name1 = cmd_info.out_file_name1_.substr(real_begin1, cmd_info.out_file_name1_.length() - real_begin1);
        //fprintf(stderr, "prefix_name1 %s\n", prefix_name1.c_str());
        //fprintf(stderr, "suffix_name1 %s\n", suffix_name1.c_str());
        cmd_info.out_file_name1_ = prefix_name1 + "tmp_" + to_string(rand()) + suffix_name1;
        if(cmd_info.out_file_name2_.length() > 0){
            int real_begin2 = 0;
            for(int i = cmd_info.out_file_name2_.length() - 1; i >= 0; i--){
                if(cmd_info.out_file_name2_[i] == '/'){
                    real_begin2 = i + 1;
                    break;
                }
            }
            prefix_name2 = cmd_info.out_file_name2_.substr(0, real_begin2);
            suffix_name2 = cmd_info.out_file_name2_.substr(real_begin2, cmd_info.out_file_name2_.length() - real_begin2);
            //fprintf(stderr, "prefix_name2 %s\n", prefix_name2.c_str());
            //fprintf(stderr, "suffix_name2 %s\n", suffix_name2.c_str());
            cmd_info.out_file_name2_ = prefix_name2 + "tmp_" + to_string(rand()) + suffix_name2;
        }
        //fprintf(stderr, "===out1 is %s\n", cmd_info.out_file_name1_.c_str());
        if(cmd_info.out_file_name2_.length() > 0){
            //fprintf(stderr, "===out2 is %s\n", cmd_info.out_file_name2_.c_str());
        }
    }

    if (cmd_info.no_trim_adapter_) cmd_info.trim_adapter_ = false;
    ASSERT(cmd_info.no_trim_adapter_ != cmd_info.trim_adapter_);
    if (!cmd_info.trim_adapter_) {
        cmd_info.se_auto_detect_adapter_ = false;
        cmd_info.pe_auto_detect_adapter_ = false;
        cmd_info.adapter_seq1_ = "";
        cmd_info.adapter_seq2_ = "";
        cmd_info.adapter_fasta_file_ = "";
        fprintf(stderr, "no adapter trim (ignore '--adapterSeq*' and '--adapterFastaFile' options) because using the '-a (--noTrimAdapter)' option!\n");
    }

    if (cmd_info.trim_5end_) {
        fprintf(stderr, "now do 5end trim\n");
    }
    if (cmd_info.trim_3end_) {
        fprintf(stderr, "now do 3end trim\n");
    }
    if (cmd_info.trim_polyg_) {
        fprintf(stderr, "now do polyg trim\n");
    }
    if (cmd_info.trim_polyx_) {
        fprintf(stderr, "now do polyx trim\n");
    }


    if (cmd_info.add_umi_) {
        fprintf(stderr, "now doing umi add\n");
        fprintf(stderr, "umi location is %s\n", umiLoc.c_str());
        fprintf(stderr, "umi len is %d\n", cmd_info.umi_len_);
        fprintf(stderr, "umi skip is %d\n", cmd_info.umi_skip_);
        fprintf(stderr, "umi prefix is %s\n", cmd_info.umi_prefix_.c_str());

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
        fprintf(stderr, "now doing overrepresentation\n");
        fprintf(stderr, "overrepresentation sampling is %d\n", cmd_info.overrepresentation_sampling_);
    }

    //if (cmd_info.isPhred64_) {
    //    fprintf(stderr, "now use phred64 input\n");
    //}


    if (cmd_info.use_pugz_) {
        if(cmd_info.pugz_threads_ > 8){
            //fprintf(stderr, "pugz thread number must <= 8, now set pugz thread number == 8.\n");
            //cmd_info.pugz_threads_ = 8;
        }
        //fprintf(stderr, "now use pugz, pugz thread number is %d\n", cmd_info.pugz_threads_);
    }
    if (cmd_info.use_pigz_) {
        //fprintf(stderr, "now use pigz, pigz thread number is %d\n", cmd_info.pigz_threads_);
    }
    //if (cmd_info.thread_number_ == 1)
    //    fprintf(stderr, "now use %d thread to do QC operations\n", cmd_info.thread_number_);
    //else
    //    fprintf(stderr, "now use %d threads to do QC operations\n", cmd_info.thread_number_);
    int mx_len = 150;
    if(do_max_len_eval) mx_len = Adapter::EvalMaxLen(cmd_info.in_file_name1_);
    //fprintf(stderr, "auto detect max seqs len is %d\n", mx_len);
    cmd_info.seq_len_ = mx_len;
    if(cmd_info.adapter_fasta_file_.length() > 0){
        fprintf(stderr, "loading adatper from %s\n", cmd_info.adapter_fasta_file_.c_str());
        cmd_info.adapter_from_fasta_ = Adapter::LoadAdaptersFromFasta(cmd_info.adapter_fasta_file_);
        sort(cmd_info.adapter_from_fasta_.begin(), cmd_info.adapter_from_fasta_.end());
        for(auto item : cmd_info.adapter_from_fasta_)fprintf(stderr, " --- %s ---\n", item.c_str());
    }
    double t_start = GetTime();

    if (cmd_info.in_file_name2_.length() || cmd_info.interleaved_in_) {
        //PE
        if ((cmd_info.out_file_name1_.length() > 0 && cmd_info.out_file_name2_.length() > 0) ||
                cmd_info.interleaved_out_) {
            cmd_info.write_data_ = true;
        }

        //adapter
        if (cmd_info.adapter_seq1_.length()) {
            fprintf(stderr, "input adapter1 is %s\n", cmd_info.adapter_seq1_.c_str());
            if (cmd_info.adapter_seq2_.length() == 0) {
                cmd_info.adapter_seq2_ = cmd_info.adapter_seq1_;
            }
            fprintf(stderr, "input adapter2 is %s\n", cmd_info.adapter_seq2_.c_str());
            cmd_info.pe_auto_detect_adapter_ = false;
            cmd_info.detect_adapter1_ = true;
            cmd_info.detect_adapter2_ = true;
        }
        if (cmd_info.pe_auto_detect_adapter_) {
            double t2 = GetTime();
            fprintf(stderr, "now auto detect adapter\n");
            cmd_info.adapter_seq1_ = Adapter::AutoDetect(cmd_info.in_file_name1_, cmd_info.trim_tail1_);
            cmd_info.adapter_seq2_ = Adapter::AutoDetect(cmd_info.in_file_name2_, cmd_info.trim_tail1_);
            if (cmd_info.adapter_seq1_.length()) {
                fprintf(stderr, "find adapter %s in read1\n", cmd_info.adapter_seq1_.c_str());
                cmd_info.detect_adapter1_ = true;
            } else {
                fprintf(stderr, "not find adapter in read1\n");
            }
            if (cmd_info.adapter_seq2_.length()) {
                fprintf(stderr, "find adapter %s in read2\n", cmd_info.adapter_seq2_.c_str());
                cmd_info.detect_adapter2_ = true;
            } else {
                fprintf(stderr, "not find adapter in read2\n");
            }
#ifdef Verbose
            fprintf(stderr, "detect adapter cost %.5f\n", GetTime() - t2);
#endif
        }
#ifdef Verbose
        if (cmd_info.correct_data_) {
            fprintf(stderr, "now correct data\n");
        }
#endif
        if (cmd_info.trim_adapter_ || cmd_info.correct_data_ || !cmd_info.no_insert_size_) {
            cmd_info.analyze_overlap_ = true;
            fprintf(stderr, "for PE data, overlap analysis is used to find adapter by default\n");
        }

        if (cmd_info.trim_front1_) {
            fprintf(stderr, "read1 trim front %d bases\n", cmd_info.trim_front1_);
            cmd_info.trim_front2_ = cmd_info.trim_front1_;
            fprintf(stderr, "read2 trim front %d bases\n", cmd_info.trim_front2_);
        }
        if (cmd_info.trim_tail1_) {
            fprintf(stderr, "read1 trim tail %d bases\n", cmd_info.trim_tail1_);
            cmd_info.trim_tail2_ = cmd_info.trim_tail1_;
            fprintf(stderr, "read2 trim tail %d bases\n", cmd_info.trim_tail2_);
        }

        if (cmd_info.do_overrepresentation_) {
            double t2 = GetTime();
            fprintf(stderr, "now doing overrepresent preprocessing part\n");
            Adapter::PreOverAnalyze(cmd_info.in_file_name1_, cmd_info.hot_seqs_, cmd_info.eva_len_);
            Adapter::PreOverAnalyze(cmd_info.in_file_name2_, cmd_info.hot_seqs2_, cmd_info.eva_len2_);
            fprintf(stderr, "overrepresent preprocessing part done\n");
            fprintf(stderr, "read1 has %d hot sequence\n", cmd_info.hot_seqs_.size());
            fprintf(stderr, "read2 has %d hot sequence\n", cmd_info.hot_seqs2_.size());
#ifdef Verbose
            fprintf(stderr, "pre over representation cost %.5f\n", GetTime() - t2);
#endif
        }

        if (cmd_info.interleaved_in_) {
            fprintf(stderr, "now input use interleaved pe data\n");
        }
        if (cmd_info.interleaved_out_) {
            fprintf(stderr, "now output use interleaved pe data\n");
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
            fprintf(stderr, "input adapter is %s\n", cmd_info.adapter_seq1_.c_str());
            cmd_info.se_auto_detect_adapter_ = false;
            cmd_info.detect_adapter1_ = true;
        }
        if (cmd_info.se_auto_detect_adapter_) {
            double t2 = GetTime();
            fprintf(stderr, "now auto detect adapter\n");
            cmd_info.adapter_seq1_ = Adapter::AutoDetect(cmd_info.in_file_name1_, cmd_info.trim_tail1_);
            if (cmd_info.adapter_seq1_.length()) {
                fprintf(stderr, "find adapter %s\n", cmd_info.adapter_seq1_.c_str());
                cmd_info.detect_adapter1_ = true;
            } else {
                fprintf(stderr, "not find adapter\n");
            }
#ifdef Verbose
            fprintf(stderr, "detect adapter cost %.5f\n", GetTime() - t2);
#endif
        }
        if (cmd_info.trim_front1_) {
            fprintf(stderr, "trim front %d bases\n", cmd_info.trim_front1_);
        }
        if (cmd_info.trim_tail1_) {
            fprintf(stderr, "trim tail %d bases\n", cmd_info.trim_tail1_);
        }

        if (cmd_info.do_overrepresentation_) {
            double t2 = GetTime();
            fprintf(stderr, "now doing overrepresent preprocessing part\n");
            Adapter::PreOverAnalyze(cmd_info.in_file_name1_, cmd_info.hot_seqs_, cmd_info.eva_len_);
            fprintf(stderr, "overrepresent preprocessing part done\n");
            fprintf(stderr, "total %d hot sqes\n", cmd_info.hot_seqs_.size());
#ifdef Verbose
            fprintf(stderr, "pre over representation cost %.5f\n", GetTime() - t2);
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
    if(out_gz){
        string out_name1 = cmd_info.out_file_name1_;
        out_name1 = out_name1.substr(0, out_name1.find(".gz"));
        remove(out_name1.c_str());
        if(cmd_info.out_file_name2_.length() > 0){
            string out_name2 = cmd_info.out_file_name2_;
            out_name2 = out_name2.substr(0, out_name2.find(".gz"));
            remove(out_name2.c_str());
        }
        string init_name1 = prefix_name1 + suffix_name1;
        rename(cmd_info.out_file_name1_.c_str(), init_name1.c_str());
        if(cmd_info.out_file_name2_.length() > 0){
            string init_name2 = prefix_name2 + suffix_name2;
            rename(cmd_info.out_file_name2_.c_str(), init_name2.c_str());
        }
    }
    fprintf(stderr, "cmd is %s\n", command.c_str());
    fprintf(stderr, "total cost %.5fs\n", GetTime() - t_start);
    return 0;
}
