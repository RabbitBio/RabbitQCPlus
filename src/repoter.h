//
// Created by ylf9811 on 2021/7/12.
//

#ifndef RERABBITQC_REPOTER_H
#define RERABBITQC_REPOTER_H

#include "state.h"
#include "Reference.h"
#include "state.h"
#include "tgsstate.h"

class Repoter {
public:
    Repoter();


    static void PrintRef(neoReference &ref);

    static void PrintRef(Reference &ref);

    static void ReportHtmlTGS(std::string html_name, std::string command, TGSStats *tgs_stats, std::string file_name);

    static void ReportHtmlSe(std::string html_name, State *state1, State *state2, std::string file_name, double dup);

    static void
    ReportHtmlPe(std::string html_name, State *pre_state1, State *pre_state2, State *aft_state1, State *aft_state2,
                 std::string file_name1,
                 std::string file_name2, double dup, int64_t *size_info);
//    static void ReportHtmlPe(State *state1, State *state2);

    static Reference GetRevRef(neoReference &ref);

};


#endif //RERABBITQC_REPOTER_H
