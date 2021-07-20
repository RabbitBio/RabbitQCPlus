//
// Created by ylf9811 on 2021/7/11.
//

#include "adapter.h"


void PrintRef(neoReference &ref) {
    printf("%s\n", std::string((char *) ref.base + ref.pname, ref.lname).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pseq, ref.lseq).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pstrand, ref.lstrand).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pqual, ref.lqual).c_str());
}

void PrintRef(Reference &ref) {
    printf("%s\n", ref.name.c_str());
    printf("%s\n", ref.seq.c_str());
    printf("%s\n", ref.strand.c_str());
    printf("%s\n", ref.quality.c_str());
}

Reference GetRevRef(neoReference &ref) {
    Reference res;
    res.name = std::string((char *) ref.base + ref.pname, ref.lname);
    res.seq = std::string((char *) ref.base + ref.pseq, ref.lseq);
    res.strand = std::string((char *) ref.base + ref.pstrand, ref.lstrand);
    res.quality = std::string((char *) ref.base + ref.pqual, ref.lqual);
    reverse(res.quality.begin(), res.quality.end());
    std::string tmp = std::string(ref.lseq, 'A');
    for (int i = 0; i < ref.lseq; i++) {
        tmp[i] = rabbit::reMap[res.seq[ref.lseq - i - 1]];
    }
    res.seq = tmp;
    return res;
}


/**
 * @brief Analyze overlap information for pe data
 * @param r1
 * @param r2
 * @param overlap_diff_limit : the max number of different bases in overlap part
 * @param overlap_require : the min length of overlap part
 * @return bool overlaped;int offset;int overlap_len;int diff_num;
 */
OverlapRes Adapter::AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require) {

    int len1 = r1.lseq;
    int len2 = r2.lseq;
    const char *str1 = reinterpret_cast<const char *>(r1.base + r1.pseq);
    const char *str2 = reinterpret_cast<const char *>(r2.base + r2.pseq);
    int complete_compare_require = 50;
    int overlap_len = 0;
    int offset = 0;
    int diff_num = 0;

    // forward
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1 - overlap_require) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = std::min(len1 - offset, len2);

        diff_num = 0;
        int i = 0;
        for (i = 0; i < overlap_len; i++) {
            if (str1[offset + i] != rabbit::reMap[str2[len2 - 1 - i]]) {
                diff_num += 1;
                if (diff_num > overlap_diff_limit && i < complete_compare_require)
                    break;
            }
        }

        if (diff_num <= overlap_diff_limit || (diff_num > overlap_diff_limit && i > complete_compare_require)) {
            return {true, offset, overlap_len, diff_num};
        }

        offset += 1;
    }


    // reverse
    // in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    // check if distance can get smaller if offset goes negative
    // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    // we go reversely
    offset = 0;
    while (offset > -(len2 - overlap_require)) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = std::min(len1, len2 - abs(offset));

        diff_num = 0;
        int i = 0;
        for (i = 0; i < overlap_len; i++) {
            if (str1[i] != rabbit::reMap[str2[len2 - 1 + offset - i]]) {
                diff_num += 1;
                if (diff_num > overlap_diff_limit && i < complete_compare_require)
                    break;
            }
        }

        if (diff_num <= overlap_diff_limit || (diff_num > overlap_diff_limit && i > complete_compare_require)) {
            return {true, offset, overlap_len, diff_num};
        }

        offset -= 1;
    }
    return {false, 0, 0, 0};

}
/**
 * @brief Trim adapter by giving adapter sequence
 * @param ref
 * @param adapter_seq
 * @param isR2
 * @return
 */

bool Adapter::TrimAdapter(neoReference &ref, std::string &adapter_seq, bool isR2) {
    const int matchReq = 4;
    const int allowOneMismatchForEach = 8;

    int rlen = ref.lseq;
    int alen = adapter_seq.length();

    const char *adata = adapter_seq.c_str();
    const char *rdata = reinterpret_cast<const char *>(ref.base + ref.pseq);

    if (alen < matchReq)
        return false;

    int pos = 0;
    bool found = false;
    int start = 0;
    if (alen >= 16)
        start = -4;
    else if (alen >= 12)
        start = -3;
    else if (alen >= 8)
        start = -2;
    // we start from negative numbers since the Illumina adapter dimer usually have the first A skipped as A-tailing
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = std::min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for (int i = std::max(0, -pos); i < cmplen; i++) {
            if (adata[i] != rdata[i + pos]) {
                mismatch++;
                if (mismatch > allowedMismatch) {
                    matched = false;
                    break;
                }
            }
        }
        if (matched) {
            found = true;
            break;
        }

    }

    if (found) {
        if (pos < 0) {
            std::string adapter = adapter_seq.substr(0, alen + pos);
            ref.lseq = 0;
            ref.lqual = 0;
//            if (fr) {
//                fr->addAdapterTrimmed(adapter, isR2);
//            }

        } else {
            std::string new_adapter_seq = std::string(reinterpret_cast<const char *>(ref.base + ref.pseq + pos),
                                                      rlen - pos);
            ref.lseq = pos;
            ref.lqual = pos;
//            if (fr) {
//                fr->addAdapterTrimmed(adapter, isR2);
//            }
        }
        return true;
    }

    return false;
}

/**
 * @brief Trim adapter by overlap information
 * @param r1
 * @param r2
 * @param offset
 * @param overlap_len
 * @return true if successfully trim adapter, false else
 */
bool Adapter::TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len) {
//    if(ov.diff<=5 && ov.overlapped && ov.offset < 0 && ol > r1->length()/3)

    if (overlap_len > 0 && offset < 0) {
        std::string adapter1 = std::string(reinterpret_cast<const char *>(r1.base + r1.pseq + overlap_len),
                                           r1.lseq - overlap_len);
        std::string adapter2 = std::string(reinterpret_cast<const char *>(r2.base + r1.pseq + overlap_len),
                                           r2.lseq - overlap_len);

//        printf("adapter1 %s\n", adapter1.c_str());
//        printf("adapter2 %s\n", adapter2.c_str());
//        printf("overlap : %d , %d\n", offset, overlap_len);
//        PrintRef(r1);
//        Reference tmp = GetRevRef(r2);
//        PrintRef(tmp);
//        printf("\n");
        r1.lseq = overlap_len;
        r1.lqual = overlap_len;
        r2.lseq = overlap_len;
        r2.lqual = overlap_len;
//        printf("after trim adapter:\n");
//        PrintRef(r1);
//        tmp = GetRevRef(r2);
//        PrintRef(tmp);
//        printf("\n");

//TODO add adapters to some ....
//        fr->addAdapterTrimmed(adapter1, adapter2);
        return true;
    }
    return false;
}

int Adapter::CorrectData(neoReference &r1, neoReference &r2, OverlapRes &overlap_res) {
    if (overlap_res.diff_num == 0 || overlap_res.diff_num > 5)
        return 0;

    int ol = overlap_res.overlap_len;
    //TODO check right ?

    int start1 = std::max(0, overlap_res.offset);
    int start2 = r2.lseq - std::max(0, -overlap_res.offset) - 1;

    char *seq1 = reinterpret_cast< char *>(r1.base + r1.pseq);
    char *seq2 = reinterpret_cast< char *>(r2.base + r2.pseq);
    char *qual1 = reinterpret_cast< char *>(r1.base + r1.pqual);
    char *qual2 = reinterpret_cast< char *>(r2.base + r2.pqual);

    const char GOOD_QUAL = '?';//30 + 33
    const char BAD_QUAL = '/';//14 + 33

//    printf("GOOD_QUAL %d\n", GOOD_QUAL);
//    printf("BAD_QUAL %d\n", BAD_QUAL);


    int corrected = 0;
    int uncorrected = 0;
    bool r1Corrected = false;
    bool r2Corrected = false;
    for (int i = 0; i < ol; i++) {
        int p1 = start1 + i;
        int p2 = start2 - i;

        if (seq1[p1] != rabbit::reMap[seq2[p2]]) {
            if (qual1[p1] >= GOOD_QUAL && qual2[p2] <= BAD_QUAL) {
                // use R1
                seq2[p2] = rabbit::reMap[seq1[p1]];
                qual2[p2] = qual1[p1];
                corrected++;
                r2Corrected = true;
                //TODO add to filter result
//                if (fr) {
//                    fr->addCorrection(seq2[p2], complement(seq1[p1]));
//                }
            } else if (qual2[p2] >= GOOD_QUAL && qual1[p1] <= BAD_QUAL) {
                // use R2

                seq1[p1] = rabbit::reMap[seq2[p2]];
                qual1[p1] = qual2[p2];
                corrected++;
                r1Corrected = true;
//                if (fr) {
//                    fr->addCorrection(seq1[p1], complement(seq2[p2]));
//                }
            } else {
                uncorrected++;
            }
        }
    }

    // should never happen
    if (uncorrected + corrected != overlap_res.diff_num) {
        static bool warned = false;
        if (!warned) {
            printf("WARNING: the algorithm is wrong! uncorrected + corrected != overlap_res.diff_num");
            warned = true;
        }
    }

//    if (corrected > 0 && fr) {
//        if (r1Corrected && r2Corrected)
//            fr->incCorrectedReads(2);
//        else
//            fr->incCorrectedReads(1);
//    }

    return corrected;
}