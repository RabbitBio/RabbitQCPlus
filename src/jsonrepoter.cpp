#include "jsonrepoter.h"
#include <fstream>

using json = nlohmann::json;

namespace {
    json GenSummary(State *state) {
        int64_t total_reads = state->GetLines();
        int64_t total_bases = state->GetTotBases();
        int64_t q20 = state->GetQ20Bases();
        int64_t q30 = state->GetQ30Bases();
        double gc = total_bases > 0 ? (double)state->GetGcBases() / total_bases : 0.0;

        json j;
        j["total_reads"] = total_reads;
        j["total_bases"] = total_bases;
        j["q20_bases"] = q20;
        j["q30_bases"] = q30;
        j["q20_rate"] = total_bases > 0 ? (double)q20 / total_bases : 0.0;
        j["q30_rate"] = total_bases > 0 ? (double)q30 / total_bases : 0.0;
        j["read1_mean_length"] = state->GetAvgLen();
        j["gc_content"] = gc;
        return j;
    }
}

void JsonReporter::ReportSe(const std::string &json_path, State *before, State *after,
                            const std::string &file_name, double dup_rate) {
    json j;
    j["file_name"] = file_name;
    j["summary"]["before_filtering"] = GenSummary(before);
    j["summary"]["after_filtering"] = GenSummary(after);
    j["filtering_result"]["passed_filter_reads"] = after->pass_reads_;
    j["filtering_result"]["low_quality_reads"] = after->fail_lowq_;
    j["filtering_result"]["too_many_N_reads"] = after->fail_N_;
    j["filtering_result"]["too_short_reads"] = after->fail_short_;
    j["filtering_result"]["too_long_reads"] = after->fail_long_;
    j["duplication"]["rate"] = dup_rate;

    std::ofstream out(json_path);
    out << j.dump(4);
}

void JsonReporter::ReportPe(const std::string &json_path, State *pre1, State *pre2,
                            State *aft1, State *aft2,
                            const std::string &file1, const std::string &file2,
                            double dup_rate, int64_t *size_info) {
    json j;
    j["file_name_1"] = file1;
    j["file_name_2"] = file2;

    j["summary"]["before_filtering"]["read1"] = GenSummary(pre1);
    j["summary"]["before_filtering"]["read2"] = GenSummary(pre2);

    j["summary"]["after_filtering"]["read1"] = GenSummary(aft1);
    j["summary"]["after_filtering"]["read2"] = GenSummary(aft2);

    j["filtering_result"]["passed_filter_reads"] = aft1->pass_reads_;
    j["filtering_result"]["low_quality_reads"] = aft1->fail_lowq_;
    j["filtering_result"]["too_many_N_reads"] = aft1->fail_N_;
    j["filtering_result"]["too_short_reads"] = aft1->fail_short_;
    j["filtering_result"]["too_long_reads"] = aft1->fail_long_;

    j["duplication"]["rate"] = dup_rate;

    j["size_info"]["size1"] = size_info[0];
    j["size_info"]["size2"] = size_info[1];

    std::ofstream out(json_path);
    out << j.dump(4);
}
