//
// Created by ylf9811 on 2021/7/14.
//

#include "polyx.h"
#include <fstream>

PolyX::PolyX() {
}


PolyX::~PolyX() {
}

void PolyX::trimPolyG(neoReference &r1, neoReference &r2, int compareReq) {
    trimPolyG(r1, compareReq);
    trimPolyG(r2, compareReq);
}

void PolyX::trimPolyG(neoReference &ref, int compareReq) {
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;

    const char *data = reinterpret_cast<const char *>(ref.base + ref.pseq);

    int rlen = ref.lseq;

    int mismatch = 0;
    int i = 0;
    int firstGPos = rlen - 1;
    for (i = 0; i < rlen; i++) {
        if (data[rlen - i - 1] != 'G') {
            mismatch++;
        } else {
            firstGPos = rlen - i - 1;
        }

        int allowedMismatch = (i + 1) / allowOneMismatchForEach;
        if (mismatch > maxMismatch || (mismatch > allowedMismatch && i >= compareReq - 1))
            break;
    }

    if (i >= compareReq) {
        ref.lseq = firstGPos;
        ref.lqual = firstGPos;
    }
}

void PolyX::trimPolyX(neoReference &r1, neoReference &r2, int compareReq) {
    trimPolyX(r1, compareReq);
    trimPolyX(r2, compareReq);
}

void PolyX::trimPolyX(neoReference &ref, int compareReq) {
    static bool geok = 0;
    if (geok == 0) {
        geok = 1;
        fprintf(stderr, "compareReq %d\n", compareReq);
    }
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;

    const char *data = reinterpret_cast<const char *>(ref.base + ref.pseq);

    int rlen = ref.lseq;

    int atcgNumbers[4] = {0, 0, 0, 0};
    char atcgBases[4] = {'A', 'T', 'C', 'G'};
    int pos = 0;
    for (pos = 0; pos < rlen; pos++) {
        switch (data[rlen - pos - 1]) {
            case 'A':
                atcgNumbers[0]++;
                break;
            case 'T':
                atcgNumbers[1]++;
                break;
            case 'C':
                atcgNumbers[2]++;
                break;
            case 'G':
                atcgNumbers[3]++;
                break;
            case 'N':
                atcgNumbers[0]++;
                atcgNumbers[1]++;
                atcgNumbers[2]++;
                atcgNumbers[3]++;
                break;
            default:
                break;
        }

        int cmp = (pos + 1);
        int allowedMismatch = std::min(maxMismatch, cmp / allowOneMismatchForEach);

        bool needToBreak = true;
        for (int b = 0; b < 4; b++) {
            if (cmp - atcgNumbers[b] <= allowedMismatch)
                needToBreak = false;
        }
        if (needToBreak && (pos >= allowOneMismatchForEach || pos + 1 >= compareReq - 1)) {
            break;
        }
    }


    // has polyX
    if (pos + 1 >= compareReq && pos < rlen) {
        //        std::string namer = std::string(reinterpret_cast<const char *>(ref.base + ref.pname), ref.lname);
        //
        //        static int cnt = 0;
        //        static bool okk = 0;
        //        static std::ofstream out_stream_;
        //        if (okk == 0) {
        //            out_stream_.open("info.txt");
        //            okk = 1;
        //        }
        //        out_stream_ << "name " << namer << " " << cnt << " " << pos << " ";
        //        if (cnt == 45691) {
        //            out_stream_ << std::string(reinterpret_cast<const char *>(ref.base + ref.pseq), ref.lseq) << std::endl;
        //            out_stream_ << atcgNumbers[0] << " " << atcgNumbers[1] << " " << atcgNumbers[2] << " " << atcgNumbers[3]
        //                        << std::endl;
        //            out_stream_ << ref.lseq << std::endl;
        //        }
        //        cnt++;

        // find the poly
        int poly;
        int maxCount = -1;
        for (int b = 0; b < 4; b++) {
            if (atcgNumbers[b] > maxCount) {
                maxCount = atcgNumbers[b];
                poly = b;
            }
        }
        char polyBase = atcgBases[poly];
        while (data[rlen - pos - 1] != polyBase && pos >= 0)
            pos--;
        //        out_stream_ << pos << std::endl;
        ref.lseq = rlen - pos - 1;
        ref.lqual = rlen - pos - 1;
        //          r->resize(rlen - pos - 1);
    }
}
