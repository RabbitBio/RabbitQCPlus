//
// Created by ylf9811 on 2021/7/7.
//

#include "filter.h"

Filter::Filter(CmdInfo *cmd_info1) {
    cmd_info = cmd_info1;
    memset(filter_res_cnt_, 0, 16 * sizeof(int64_t));
}

/**
	@brief filtering read
	@return filtering result
    0:pass
    1:fail because too many N
    2:fail because too short
    3:fail because too long
    4:fail because too many bases < q15
 */
int Filter::ReadFiltering(neoReference &ref, bool trim_res, bool isPhred64) {
    if (!trim_res) return 2;
    int seq_len = ref.lseq;
    int qul_len = ref.lqual;
    ASSERT(seq_len == qul_len);
    if (seq_len >= cmd_info->length_limit_) {
        return 3;
    }
    if (seq_len == 0 || seq_len < cmd_info->length_required_) {
        return 2;
    }
    int n_number = 0;
    int low_qual_number = 0;// < q15
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);

    int phredSub = 33;
    if (isPhred64) phredSub = 64;
    for (int i = 0; i < seq_len; i++) {
        int q = std::max(0, quals[i] - phredSub);
        if (q < 15) {
            low_qual_number++;
        }
        if (bases[i] == 'N') {
            n_number++;
        }
    }
    if (low_qual_number > cmd_info->low_qual_perc_limit_ * seq_len / 100.0) return 4;
    if (n_number > cmd_info->n_number_limit_) return 1;
    return 0;
}


bool Filter::TrimSeq(neoReference &ref, int front, int tail) {

    int new_seq_len = ref.lseq - front - tail;
    if (new_seq_len <= 0) return false;

    int w = cmd_info->cut_window_size_;
    int l = ref.lseq;
    const char *seq = reinterpret_cast<const char *>(ref.base + ref.pseq);
    const char *qualstr = reinterpret_cast<const char *>(ref.base + ref.pqual);
    int phredSub = 33;
    if (cmd_info->isPhred64_) phredSub = 64;
    // quality cutting forward
    if (cmd_info->trim_5end_) {
        int s = front;
        if (l - front - tail - w <= 0) {
            return 0;
        }

        int totalQual = 0;

        // preparing rolling
        for (int i = 0; i < w - 1; i++)
            totalQual += qualstr[s + i];

        for (s = front; s + w < l - tail; s++) {
            totalQual += qualstr[s + w - 1];
            // rolling
            if (s > front) {
                totalQual -= qualstr[s - 1];
            }
            // add 33 for phred33 transforming
            if ((double) totalQual / (double) w >= phredSub + cmd_info->cut_mean_quality_)
                break;
        }

        // the trimming in front is forwarded and rlen is recalculated
        if (s > 0)
            s = s + w - 1;
        while (s < l && seq[s] == 'N')
            s++;
        front = s;
        new_seq_len = l - front - tail;
    }

    // quality cutting backward
    if (cmd_info->trim_3end_) {
        if (l - front - tail - w <= 0) {
            return 0;
        }

        int totalQual = 0;
        int t = l - tail - 1;

        // preparing rolling
        for (int i = 0; i < w - 1; i++)
            totalQual += qualstr[t - i];

        for (t = l - tail - 1; t - w >= front; t--) {
            totalQual += qualstr[t - w + 1];
            // rolling
            if (t < l - tail - 1) {
                totalQual -= qualstr[t + 1];
            }
            // add 33 for phred33 transforming
            if ((double) totalQual / (double) w >= phredSub + cmd_info->cut_mean_quality_)
                break;
        }

        if (t < l - 1)
            t = t - w + 1;
        while (t >= 0 && seq[t] == 'N')
            t--;
        new_seq_len = t - front + 1;
    }

    if (new_seq_len <= 0 || front >= l - 1) {
        return 0;
    }


    ref.pseq += front;
    ref.pqual += front;
    ref.lseq = new_seq_len;
    ref.lqual = new_seq_len;

    return true;
}


/**
 * @brief Print filtering result, for example xx refs because too short......
 */
void Filter::PrintResult() {
}
