//
// Created by ylf9811 on 2021/7/7.
//

#include "filter.h"

Filter::Filter(CmdInfo *cmd_info1) {
    cmd_info = cmd_info1;
}

/**
	@brief filtering read
	@return filtering result
    0:pass
    1:fail because too short
    2:fail because too long
    3:fail because too many bases < q15
    4:fail because too many N
 */
int Filter::ReadFiltering(neoReference &ref) {
    int seq_len = ref.lseq;
    int qul_len = ref.lqual;
    ASSERT(seq_len == qul_len);
    if (seq_len == 0 || seq_len < cmd_info->base_len_limit_) {
        return 1;
    }
    int n_number = 0;
    int low_qual_number = 0;// < q15
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);
    const char q15 = '0';

    for (int i = 0; i < seq_len; i++) {
        if (quals[i] < q15) {
            low_qual_number++;
        }
        if (bases[i] == 'N') {
            n_number++;
        }
    }
    if (low_qual_number > cmd_info->low_qual_perc_limit_ * seq_len / 100.0)return 3;
    if (n_number > cmd_info->n_number_limit_)return 4;
    return 0;

}

