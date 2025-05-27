//
// Created by ylf9811 on 2025/5/27.
//

#ifndef JSON_RQCP_JSONREPOTER_H
#define JSON_RQCP_JSONREPOTER_H

#include <string>
#include "state.h"
#include <nlohmann/json.hpp>

class JsonReporter {
public:
    static void ReportSe(const std::string &json_path, State *before, State *after,
                         const std::string &file_name, double dup_rate);

    static void ReportPe(const std::string &json_path, State *pre1, State *pre2,
                         State *aft1, State *aft2,
                         const std::string &file1, const std::string &file2,
                         double dup_rate, int64_t *size_info);
};


#endif //JSON_RQCP_JSONREPOTER_H
