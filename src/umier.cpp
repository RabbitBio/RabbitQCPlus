//
// Created by ylf9811 on 2021/7/20.
//

#include "umier.h"

Umier::Umier(CmdInfo *cmd_info) {
    cmd_info_ = cmd_info;
}


Umier::~Umier() {
}

std::string Umier::firstIndex(std::string name) {
    int len = name.length();
    int end = len;
    if (len < 5)
        return "";
    for (int i = len - 3; i >= 0; i--) {
        if (name[i] == '+')
            end = i - 1;
        if (name[i] == ':') {
            return name.substr(i + 1, end - i);
        }
    }
    return "";
}

std::string Umier::lastIndex(std::string name) {
    int len = name.length();
    if (len < 5)
        return "";
    for (int i = len - 3; i >= 0; i--) {
        if (name[i] == ':' || name[i] == '+') {
            return name.substr(i + 1, len - i);
        }
    }
    return "";
}

void Umier::ProcessSe(neoReference &ref) {
    if (!cmd_info_->add_umi_)
        return;
    int seq_len = ref.lseq;
    std::string umi;
    std::string name = std::string(reinterpret_cast<const char *>(ref.base + ref.pname), ref.lname);
    std::string seq = std::string(reinterpret_cast<const char *>(ref.base + ref.pseq), ref.lseq);
    if (cmd_info_->umi_loc_ == UMI_LOC_INDEX1) {
        umi = firstIndex(name);
    } else if (cmd_info_->umi_loc_ == UMI_LOC_INDEX2) {
    } else if (cmd_info_->umi_loc_ == UMI_LOC_READ1) {
        umi = seq.substr(0, std::min(seq_len, cmd_info_->umi_len_));
        int trim_len = std::min(seq_len - 1, int(umi.length() + cmd_info_->umi_skip_));
        ref.pseq += trim_len;
        ref.pqual += trim_len;
        ref.lseq -= trim_len;
        ref.lqual -= trim_len;
    } else if (cmd_info_->umi_loc_ == UMI_LOC_READ2) {
    } else if (cmd_info_->umi_loc_ == UMI_LOC_PER_INDEX) {
        std::string umiMerged = firstIndex(name);
        addUmiToName(ref, umiMerged);
    } else if (cmd_info_->umi_loc_ == UMI_LOC_PER_READ) {
        umi = seq.substr(0, std::min(seq_len, cmd_info_->umi_len_));
        int trim_len = std::min(seq_len, int(umi.length() + cmd_info_->umi_skip_));
        ref.pseq += trim_len;
        ref.pqual += trim_len;
        ref.lseq -= trim_len;
        ref.lqual -= trim_len;
        addUmiToName(ref, umi);
    }

    if (cmd_info_->umi_loc_ != UMI_LOC_PER_INDEX && cmd_info_->umi_loc_ != UMI_LOC_PER_READ) {
        if (!umi.empty())
            addUmiToName(ref, umi);
    }
}

void Umier::ProcessPe(neoReference &r1, neoReference &r2) {
    if (!cmd_info_->add_umi_)
        return;
    int seq_len1 = r1.lseq;
    int seq_len2 = r2.lseq;
    std::string umi;
    std::string name1 = std::string(reinterpret_cast<const char *>(r1.base + r1.pname), r1.lname);
    std::string name2 = std::string(reinterpret_cast<const char *>(r2.base + r2.pname), r2.lname);
    std::string seq1 = std::string(reinterpret_cast<const char *>(r1.base + r1.pseq), r1.lseq);
    std::string seq2 = std::string(reinterpret_cast<const char *>(r2.base + r2.pseq), r2.lseq);

    if (cmd_info_->umi_loc_ == UMI_LOC_INDEX1) {
        umi = firstIndex(name1);
    } else if (cmd_info_->umi_loc_ == UMI_LOC_INDEX2) {
        umi = lastIndex(name2);
    } else if (cmd_info_->umi_loc_ == UMI_LOC_READ1) {
        umi = seq1.substr(0, std::min(seq_len1, cmd_info_->umi_len_));
        int trim_len = std::min(seq_len1 - 1, int(umi.length() + cmd_info_->umi_skip_));
        r1.pseq += trim_len;
        r1.pqual += trim_len;
        r1.lseq -= trim_len;
        r1.lqual -= trim_len;
    } else if (cmd_info_->umi_loc_ == UMI_LOC_READ2) {
        umi = seq2.substr(0, std::min(seq_len2, cmd_info_->umi_len_));
        int trim_len = std::min(seq_len1 - 1, int(umi.length() + cmd_info_->umi_skip_));
        r2.pseq += trim_len;
        r2.pqual += trim_len;
        r2.lseq -= trim_len;
        r2.lqual -= trim_len;
    } else if (cmd_info_->umi_loc_ == UMI_LOC_PER_INDEX) {
        std::string umiMerged = firstIndex(name1);
        umiMerged = umiMerged + "_" + lastIndex(name2);
        addUmiToName(r1, umiMerged);
        addUmiToName(r2, umiMerged);

    } else if (cmd_info_->umi_loc_ == UMI_LOC_PER_READ) {
        std::string umi1 = seq1.substr(0, std::min(seq_len1, cmd_info_->umi_len_));
        std::string umiMerged = umi1;
        int trim_len1 = std::min(seq_len1 - 1, int(umi1.length() + cmd_info_->umi_skip_));
        r1.pseq += trim_len1;
        r1.pqual += trim_len1;
        r1.lseq -= trim_len1;
        r1.lqual -= trim_len1;
        std::string umi2 = seq2.substr(0, std::min(seq_len2, cmd_info_->umi_len_));
        umiMerged = umiMerged + "_" + umi2;
        int trim_len2 = std::min(seq_len2 - 1, int(umi2.length() + cmd_info_->umi_skip_));
        r2.pseq += trim_len2;
        r2.pqual += trim_len2;
        r2.lseq -= trim_len2;
        r2.lqual -= trim_len2;
        addUmiToName(r1, umiMerged);
        addUmiToName(r2, umiMerged);
    }

    if (cmd_info_->umi_loc_ != UMI_LOC_PER_INDEX && cmd_info_->umi_loc_ != UMI_LOC_PER_READ) {
        if (!umi.empty())
            addUmiToName(r1, umi);
        if (!umi.empty())
            addUmiToName(r2, umi);
    }
}

void Umier::addUmiToName(neoReference &ref, std::string umi) {
    std::string tag;
    if (cmd_info_->umi_prefix_.empty())
        tag = ":" + umi;
    else
        tag = ":" + cmd_info_->umi_prefix_ + "_" + umi;
    int spacePos = -1;
    const char *name = reinterpret_cast<const char *>(ref.base + ref.pname);
    for (int i = 0; i < ref.lname; i++) {
        if (name[i] == ' ') {
            spacePos = i;
            break;
        }
    }
    if (spacePos == -1) {
        char *new_name = new char[ref.lname + tag.length()];
        int pos = 0;
        memcpy(new_name + pos, ref.base + ref.pname, ref.lname);
        pos += ref.lname;
        memcpy(new_name + pos, tag.c_str(), tag.length());
        ref.lname = ref.lname + tag.length();
        ref.pname = reinterpret_cast<uint64_t>(new_name) - reinterpret_cast<uint64_t>(ref.base);
    } else {
        char *new_name = new char[ref.lname + tag.length()];
        int pos = 0;
        memcpy(new_name + pos, ref.base + ref.pname, spacePos);
        pos += spacePos;
        memcpy(new_name + pos, tag.c_str(), tag.length());
        pos += tag.length();
        memcpy(new_name + pos, ref.base + ref.pname + spacePos, ref.lname - spacePos);
        ref.lname = ref.lname + tag.length();
        ref.pname = reinterpret_cast<uint64_t>(new_name) - reinterpret_cast<uint64_t>(ref.base);
    }
}
