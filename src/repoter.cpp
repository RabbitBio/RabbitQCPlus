//
// Created by ylf9811 on 2021/7/12.
//

#include "repoter.h"

void Repoter::PrintRef(neoReference &ref) {
    printf("%s\n", std::string((char *) ref.base + ref.pname, ref.lname).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pseq, ref.lseq).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pstrand, ref.lstrand).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pqual, ref.lqual).c_str());
}

void Repoter::PrintRef(Reference &ref) {
    printf("%s\n", ref.name.c_str());
    printf("%s\n", ref.seq.c_str());
    printf("%s\n", ref.strand.c_str());
    printf("%s\n", ref.quality.c_str());
}

Reference Repoter::GetRevRef(neoReference &ref) {
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
