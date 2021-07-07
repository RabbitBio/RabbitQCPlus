//
// Created by ylf9811 on 2021/7/6.
//

#ifndef RABBITQCPLUS_CMDINFO_H
#define RABBITQCPLUS_CMDINFO_H

#include <iostream>

using namespace std;

class CmdInfo {
public:
    CmdInfo();

public:
    string in_file_name1_;
    string in_file_name2_;
    string out_file_name1_;
    string out_file_name2_;
    bool no_adapter_detect_;
    int thread_number_;
    int n_number_limit_;
    int low_qual_perc_limit_;
    int base_len_limit_;
    bool write_data_;

};


#endif //RABBITQCPLUS_CMDINFO_H
