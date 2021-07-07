#include <sys/time.h>
#include <string>
#include <iostream>
#include "FastxStream.h"
#include "CLI11.hpp"
#include "seqc.h"
#include "peqc.h"
#include "cmdinfo.h"

double GetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}


int main(int argc, char **argv) {

    CmdInfo cmd_info;
    CLI::App app("RabbitQCPlus");
    auto opt = app.add_option("-i,--inFile1", cmd_info.in_file_name1_, "input fastq name 1, can not be ''");
    opt->required();
    app.add_option("-I,--inFile2", cmd_info.in_file_name2_, "input fastq name 2, can be '' when single data");
    app.add_option("-o,--outFile1", cmd_info.out_file_name1_, "output fastq name 1");
    app.add_option("-O,--outFile2", cmd_info.out_file_name2_, "output fastq name 2");
    app.add_flag("-a,--noAda", cmd_info.no_adapter_detect_, "do not detect adapter");
    app.add_option("-w,--threadNum", cmd_info.thread_number_, "number thread used to solve fastq data");

    CLI11_PARSE(app, argc, argv);

    std::cout << "in1 is " << cmd_info.in_file_name1_ << std::endl;
    if (cmd_info.in_file_name2_.length())std::cout << "in2 is " << cmd_info.in_file_name2_ << std::endl;
    if (cmd_info.out_file_name1_.length())std::cout << "out1 is " << cmd_info.out_file_name1_ << std::endl;
    if (cmd_info.out_file_name2_.length())std::cout << "out2 is " << cmd_info.out_file_name2_ << std::endl;
    if (cmd_info.no_adapter_detect_)std::cout << "now don't detect adapter" << std::endl;
    else std::cout << "now detect adapter" << std::endl;
    std::cout << "now use " << cmd_info.thread_number_ << " thread" << (cmd_info.thread_number_ > 1 ? "s" : "")
              << std::endl;
    double t1 = GetTime();
    if (cmd_info.in_file_name2_.length()) {
        PeQc pe_qc(&cmd_info);
        pe_qc.ProcessPeFastq();
    } else {
        SeQc se_qc(&cmd_info);
        se_qc.ProcessSeFastq();
    }
    printf("cost %.5f\n", GetTime() - t1);

    return 0;


}