#include "adapter.h"

#define mod 9999991

#ifdef Vec512

#include <immintrin.h>

#endif

#ifdef Vec256

#include <immintrin.h>

#endif

using namespace std;

struct AdapterSeedInfo {
    int recordsID;
    int pos;
};

inline map<string, string> getKnownAdapter() {
    map<string, string> knownAdapters;

    knownAdapters["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"] = "Illumina TruSeq Adapter Read 1";
    knownAdapters["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"] = "Illumina TruSeq Adapter Read 2";
    knownAdapters["GATCGTCGGACTGTAGAACTCTGAACGTGTAGA"] = "Illumina Small RNA Adapter Read 2";

    return knownAdapters;
}

void PrintRef(neoReference &ref) {
    printf("%s\n", string((char *) ref.base + ref.pname, ref.lname).c_str());
    printf("%s\n", string((char *) ref.base + ref.pseq, ref.lseq).c_str());
    printf("%s\n", string((char *) ref.base + ref.pstrand, ref.lstrand).c_str());
    printf("%s\n", string((char *) ref.base + ref.pqual, ref.lqual).c_str());
}

void PrintRef(Reference &ref) {
    printf("%s\n", ref.name.c_str());
    printf("%s\n", ref.seq.c_str());
    printf("%s\n", ref.strand.c_str());
    printf("%s\n", ref.quality.c_str());
}

string Reverse(string seq) {
    string res(seq);
    reverse(res.begin(), res.end());
    return res;
}

Reference GetRevRef(neoReference &ref) {
    Reference res;
    res.name = string((char *) ref.base + ref.pname, ref.lname);
    res.seq = string((char *) ref.base + ref.pseq, ref.lseq);
    res.strand = string((char *) ref.base + ref.pstrand, ref.lstrand);
    res.quality = string((char *) ref.base + ref.pqual, ref.lqual);
    reverse(res.quality.begin(), res.quality.end());
    string tmp = string(ref.lseq, 'A');
    for (int i = 0; i < ref.lseq; i++) {
        tmp[i] = rabbit::reMap[res.seq[ref.lseq - i - 1]];
    }
    res.seq = tmp;
    return res;
}


string Adapter::int2seq(unsigned int val, int seqlen) {
    char bases[4] = {'A', 'T', 'C', 'G'};
    string ret(seqlen, 'N');
    int done = 0;
    while (done < seqlen) {
        ret[seqlen - done - 1] = bases[val & 0x03];
        val = (val >> 2);
        done++;
    }
    return ret;
}


int Adapter::seq2int(const char *seq, int pos, int keylen, int lastVal) {
    if (lastVal >= 0) {
        const int mask = (1 << (keylen * 2)) - 1;
        int key = (lastVal << 2) & mask;
        char base = seq[pos + keylen - 1];
        int val = valAGCT[base & 0x07];
        if (val == -1) return -1;
        key += val;
        return key;
    } else {
        int key = 0;
        for (int i = pos; i < keylen + pos; i++) {
            key = (key << 2);
            char base = seq[i];
            int val = valAGCT[base & 0x07];
            if (val == -1) return -1;
            key += val;
        }
        return key;
    }
}


string Adapter::matchKnownAdapter(string seq) {
    map<string, string> knownAdapters = getKnownAdapter();
    map<string, string>::iterator iter;
    for (iter = knownAdapters.begin(); iter != knownAdapters.end(); iter++) {
        string adapter = iter->first;
        string desc = iter->second;
        if (seq.length() < adapter.length()) {
            continue;
        }
        int diff = 0;
        for (int i = 0; i < adapter.length() && i < seq.length(); i++) {
            if (adapter[i] != seq[i])
                diff++;
        }
        if (diff == 0)
            return adapter;
    }
    return "";
}

void Adapter::PreOverAnalyze(string file_name, vector<string> &hot_seqs, int &eva_len) {

#ifdef Verbose
    double t0 = GetTime();
#endif
    auto *fastq_data_pool = new rabbit::fq::FastqDataPool(128, 1 << 22);
    auto fqFileReader = new rabbit::fq::FastqFileReader(file_name, *fastq_data_pool, "",
            file_name.find(".gz") != string::npos);
    int64_t n_chunks = 0;
    const long BASE_LIMIT = 151 * 10000;
    long records = 0;
    long bases = 0;
    int maxn = 1e7;
    int maxm = 1e7;
    int *head = new int[maxn];
    int *pre = new int[maxm];
    int *cnt = new int[maxm];
    int64_t *v = new int64_t[maxm];
    string *seqs = new string[maxm];
    int num = 0;
    for (int i = 0; i < maxn; i++) head[i] = -1;


    vector<neoReference> loadedReads;
    vector<rabbit::fq::FastqDataChunk *> chunks;

    while (bases < BASE_LIMIT) {
        rabbit::fq::FastqDataChunk *fqdatachunk;
        fqdatachunk = fqFileReader->readNextChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        vector<neoReference> data;
        rabbit::fq::chunkFormat(fqdatachunk, data, true);
        chunks.push_back(fqdatachunk);
        for (auto item: data) {
            bases += item.lseq;
            loadedReads.push_back(item);
            records++;
            if (bases >= BASE_LIMIT)
                break;
        }
    }
#ifdef Verbose
    t0 = GetTime();
#endif

    int seqlen = 0;
    for (int i = 0; i < loadedReads.size(); i++) {
        auto item = loadedReads[i];
        seqlen = max(seqlen, int(item.lseq));
    }
    eva_len = seqlen;


    for (int it = 0; it < loadedReads.size(); it++) {
        auto item = loadedReads[it];
        int rlen = item.lseq;

        int steps[5] = {10, 20, 40, 100, min(150, seqlen - 2)};

        for (int s = 0; s < 5; s++) {
            int step = steps[s];
            for (int offset = 0; offset < rlen - step; offset++) {
                string seq = string(reinterpret_cast<const char *>(item.base + item.pseq + offset), step);
                int64_t now = 0;
                for (int j = 0; j < step; j++) {
                    now = now * 5;
                    now += valAGCT2[seq[j] & 0x07];
                }
                int ha = (now % mod + mod) % mod;
                int ok_find = 0;
                for (int i = head[ha]; i != -1; i = pre[i]) {
                    if (now == v[i]) {
                        ok_find = 1;
                        //#pragma omp critical
                        {
                            cnt[i]++;
                        }
                        break;
                    }
                }
                if (ok_find == 0) {
                    //#pragma omp critical
                    {
                        v[num] = now;
                        pre[num] = head[ha];
                        cnt[num] = 1;
                        seqs[num] = seq;
                        head[ha] = num;
                        num++;
                    }
                    if (num >= maxm) {
                        printf("hot sequence is too many\n");
                        exit(0);
                    }
                }
            }
        }
    }
#ifdef Verbose
    t0 = GetTime();
#endif

    vector<pair<string, int>> hot_s;

    for (int i = 0; i < num; i++) {
        if (seqs[i].length() >= seqlen - 1) {
            if (cnt[i] >= 3) {
                hot_s.push_back({seqs[i], cnt[i]});
            }
        } else if (seqs[i].length() >= 100) {
            if (cnt[i] >= 5) {
                hot_s.push_back({seqs[i], cnt[i]});
            }
        } else if (seqs[i].length() >= 40) {
            if (cnt[i] >= 20) {
                hot_s.push_back({seqs[i], cnt[i]});
            }
        } else if (seqs[i].length() >= 20) {
            if (cnt[i] >= 100) {
                hot_s.push_back({seqs[i], cnt[i]});
            }
        } else if (seqs[i].length() >= 10) {
            if (cnt[i] >= 500) {
                hot_s.push_back({seqs[i], cnt[i]});
            }
        }
    }

#ifdef Verbose
    t0 = GetTime();
#endif
    // remove substrings
    for (int i = 0; i < hot_s.size(); i++) {
        auto item = hot_s[i];
        auto seq = item.first;
        auto count = item.second;
        bool isSubString = false;
        for (auto item: hot_s) {
            auto seq2 = item.first;
            auto count2 = item.second;
            if (seq != seq2 && seq2.find(seq) != string::npos && count / count2 < 10) {
                isSubString = true;
                break;
            }
        }
        if (!isSubString) {
            {
                hot_seqs.push_back(seq);
            }
        }
    }
#ifdef Verbose
    t0 = GetTime();
#endif
#ifdef Verbose
    t0 = GetTime();
#endif
    for (auto item: chunks)
        fastq_data_pool->Release(item);
    delete fastq_data_pool;
    delete fqFileReader;
    delete[] head;
    delete[] pre;
    delete[] v;
    delete[] cnt;
    delete[] seqs;
}


vector<string> Adapter::LoadAdaptersFromFasta(string file_name) {
    ifstream iff;
    iff.open(file_name);
    string name, seq;
    vector<string> res;
    while(iff >> name >> seq) {
        res.push_back(seq);
    }
    iff.close();
    return res;
}

int Adapter::EvalMaxLen(string file_name) {
    auto *fastq_data_pool = new rabbit::fq::FastqDataPool(128, 1 << 22);
    auto fqFileReader = new rabbit::fq::FastqFileReader(file_name, *fastq_data_pool, "",
            file_name.find(".gz") != string::npos);
    int64_t n_chunks = 0;
    // stat up to 256K reads
    const long READ_LIMIT = 256 * 1024;
    const long BASE_LIMIT = 151 * READ_LIMIT;
    long records = 0;
    long bases = 0;
    int mx_len = 0;
    vector<rabbit::fq::FastqDataChunk *> chunks;

    while (records < READ_LIMIT && bases < BASE_LIMIT) {
        rabbit::fq::FastqDataChunk *fqdatachunk;
        fqdatachunk = fqFileReader->readNextChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        vector<neoReference> data;
        rabbit::fq::chunkFormat(fqdatachunk, data, true);
        chunks.push_back(fqdatachunk);
        for (auto item: data) {
            bases += item.lseq;
            mx_len = max(mx_len, (int) item.lseq);
            records++;
            if (records >= READ_LIMIT || bases >= BASE_LIMIT)
                break;
        }
    }
    for (auto item: chunks)
        fastq_data_pool->Release(item);
    delete fastq_data_pool;
    delete fqFileReader;
    return mx_len;
}

string Adapter::AutoDetect(string file_name, int trim_tail) {
    auto *fastq_data_pool = new rabbit::fq::FastqDataPool(128, 1 << 22);
    auto fqFileReader = new rabbit::fq::FastqFileReader(file_name, *fastq_data_pool, "",
            file_name.find(".gz") != string::npos);
    int64_t n_chunks = 0;
    // stat up to 256K reads
    const long READ_LIMIT = 256 * 1024;
    const long BASE_LIMIT = 151 * READ_LIMIT;
    long records = 0;
    long bases = 0;


    vector<neoReference> loadedReads;
    vector<rabbit::fq::FastqDataChunk *> chunks;

    while (records < READ_LIMIT && bases < BASE_LIMIT) {
        rabbit::fq::FastqDataChunk *fqdatachunk;
        fqdatachunk = fqFileReader->readNextChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        vector<neoReference> data;
        rabbit::fq::chunkFormat(fqdatachunk, data, true);
        chunks.push_back(fqdatachunk);
        //        fastq_data_pool->Release(fqdatachunk);
        for (auto item: data) {
            bases += item.lseq;
            loadedReads.push_back(item);
            records++;
            if (records >= READ_LIMIT || bases >= BASE_LIMIT)
                break;
        }
    }


    // we need at least 10000 valid records to evaluate
    if (records < 10000) {
        return "";
    }

    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, trim_tail);

    // why we add trim_tail here? since the last cycle are usually with low quality and should be trimmed
    const int keylen = 10;
    int size = 1 << (keylen * 2);
    unsigned int *counts = new unsigned int[size];
    memset(counts, 0, sizeof(unsigned int) * size);

    for (int i = 0; i < records; i++) {
        neoReference r = loadedReads[i];
        int key = -1;
        for (int pos = 20; pos <= r.lseq - keylen - shiftTail; pos++) {
            key = seq2int(reinterpret_cast<const char *>(r.base + r.pseq), pos, keylen, key);
            if (key >= 0) {
                counts[key]++;
            }
        }
    }
    counts[0] = 0;

    // get the top N
    const int topnum = 10;
    int topkeys[topnum] = {0};
    long total = 0;
    for (int k = 0; k < size; k++) {
        int atcg[4] = {0};
        for (int i = 0; i < keylen; i++) {
            int baseOfBit = (k >> (i * 2)) & 0x03;
            atcg[baseOfBit]++;
        }
        bool lowComplexity = false;
        for (int b = 0; b < 4; b++) {
            if (atcg[b] >= keylen - 4)
                lowComplexity = true;
        }
        if (lowComplexity)
            continue;
        // too many GC
        if (atcg[2] + atcg[3] >= keylen - 2)
            continue;

        // starts with GGGG
        if (k >> 12 == 0xff)
            continue;

        unsigned int val = counts[k];
        total += val;
        for (int t = topnum - 1; t >= 0; t--) {
            // reach the middle
            if (val < counts[topkeys[t]]) {
                if (t < topnum - 1) {
                    for (int m = topnum - 1; m > t + 1; m--) {
                        topkeys[m] = topkeys[m - 1];
                    }
                    topkeys[t + 1] = k;
                }
                break;
            } else if (t == 0) {// reach the top
                for (int m = topnum - 1; m > t; m--) {
                    topkeys[m] = topkeys[m - 1];
                }
                topkeys[t] = k;
            }
        }
    }
    const int FOLD_THRESHOLD = 20;
    for (int t = 0; t < topnum; t++) {
        int key = topkeys[t];
        string seq = int2seq(key, keylen);
        if (key == 0)
            continue;
        long count = counts[key];
        if (count < 10 || count * size < total * FOLD_THRESHOLD)
            break;
        // skip low complexity seq
        int diff = 0;
        for (int s = 0; s < seq.length() - 1; s++) {
            if (seq[s] != seq[s + 1])
                diff++;
        }
        if (diff < 3) {
            continue;
        }
        string adapter = getAdapterWithSeed(key, loadedReads, records, keylen, trim_tail);
        if (!adapter.empty()) {
            delete[] counts;
            return adapter;
            return adapter;
        }
    }

    for (auto item: chunks)
        fastq_data_pool->Release(item);

    delete fastq_data_pool;
    delete fqFileReader;
    delete[] counts;
    return "";
}


string
Adapter::getAdapterWithSeed(int seed, vector<neoReference> loadedReads, long records, int keylen,
        int trim_tail) {
    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, trim_tail);
    NucleotideTree *forwardTree = new NucleotideTree();
    NucleotideTree *backwardTree = new NucleotideTree();
    vector<AdapterSeedInfo> vec;
    //#pragma omp parallel for
    for (int i = 0; i < records; i++) {
        neoReference r = loadedReads[i];
        struct AdapterSeedInfo seedInfo;
        int key = -1;
        for (int pos = 20; pos <= r.lseq - keylen - shiftTail; pos++) {
            key = seq2int(reinterpret_cast<const char *>(r.base + r.pseq), pos, keylen, key);
            if (key == seed) {
                seedInfo.recordsID = i;
                seedInfo.pos = pos;
                //#pragma omp critical
                {
                    vec.push_back(seedInfo);
                }
            }
        }
    }
    vector<AdapterSeedInfo>::iterator it;
    for (it = vec.begin(); it != vec.end(); it++) {

        auto now_it = loadedReads[it->recordsID];
        forwardTree->addSeq(
                string(reinterpret_cast<const char *>(now_it.base + now_it.pseq + it->pos + keylen),
                    now_it.lseq - keylen - shiftTail - it->pos));
        string seq = string(reinterpret_cast<const char *>(now_it.base + now_it.pseq), it->pos);
        string rcseq = Reverse(seq);
        backwardTree->addSeq(rcseq);
    }

    bool reachedLeaf = true;
    string forwardPath = forwardTree->getDominantPath(reachedLeaf);
    string backwardPath = backwardTree->getDominantPath(reachedLeaf);

    string adapter = Reverse(backwardPath) + int2seq(seed, keylen) + forwardPath;
    if (adapter.length() > 60)
        adapter.resize(60);
    delete forwardTree;
    delete backwardTree;
    string matchedAdapter = matchKnownAdapter(adapter);
    if (!matchedAdapter.empty()) {
        map<string, string> knownAdapters = getKnownAdapter();
        //cout << knownAdapters[matchedAdapter] << ": " << matchedAdapter << endl;
        return matchedAdapter;
    } else {
        if (reachedLeaf) {
            //cout << adapter << endl;
            return adapter;
        } else {
            return "";
        }
    }
}


/**
 * @brief Analyze overlap information for pe data
 * @param r1
 * @param r2
 * @param overlap_diff_limit : the max number of different bases in overlap part
 * @param overlap_require : the min length of overlap part
 * @return bool overlaped;int offset;int overlap_len;int diff_num;
 */


#ifdef Vec512

OverlapRes Adapter::AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require) {


    int len1 = r1.lseq;
    int len2 = r2.lseq;

    char *tmpRev = new char[len2];
    const char *str2 = reinterpret_cast<const char *>(r2.base + r2.pseq);

    for (int i = 0; i < len2; i++) {
        tmpRev[i] = rabbit::reMap[str2[len2 - i - 1]];
    }
    str2 = tmpRev;

    const char *str1 = reinterpret_cast<const char *>(r1.base + r1.pseq);


    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff_num = 0;

    // forward
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1 - overlap_require) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);
        diff_num = 0;
        int leng = min(complete_compare_require, overlap_len);
        for (int i = 0; i < leng; i += 64) {
            uint64_t tag = (1ll << min(64, leng - i)) - 1;
            __m512i t1 = _mm512_maskz_loadu_epi8(tag, str1 + offset + i);
            __m512i t2 = _mm512_maskz_loadu_epi8(tag, str2 + i);
            __mmask64 res = _mm512_cmp_epi8_mask(t1, t2, 4);
            diff_num += _mm_popcnt_u64(res);
            if (diff_num > overlap_diff_limit) break;
        }

        if (diff_num <= overlap_diff_limit) {
            delete[] tmpRev;
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
        overlap_len = min(len1, len2 - abs(offset));

        diff_num = 0;

        int leng = min(complete_compare_require, overlap_len);
        for (int i = 0; i < leng; i += 64) {
            uint64_t tag = (1ll << min(64, leng - i)) - 1;
            __m512i t1 = _mm512_maskz_loadu_epi8(tag, str1 + i);
            __m512i t2 = _mm512_maskz_loadu_epi8(tag, str2 - offset + i);
            __mmask64 res = _mm512_cmp_epi8_mask(t1, t2, 4);
            diff_num += _mm_popcnt_u64(res);
            if (diff_num > overlap_diff_limit) break;
        }
        if (diff_num <= overlap_diff_limit) {
            delete[] tmpRev;
            return {true, offset, overlap_len, diff_num};
        }
        offset -= 1;
    }
    delete[] tmpRev;
    return {false, 0, 0, 0};
}

#elif Vec256

OverlapRes Adapter::AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require) {

    int len1 = r1.lseq;
    int len2 = r2.lseq;

    char *tmpRev = new char[len2];
    const char *str2 = reinterpret_cast<const char *>(r2.base + r2.pseq);

    for (int i = 0; i < len2; i++) {
        tmpRev[i] = rabbit::reMap[str2[len2 - i - 1]];
    }
    str2 = tmpRev;

    const char *str1 = reinterpret_cast<const char *>(r1.base + r1.pseq);

    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff_num = 0;

    // forward
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1 - overlap_require) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);
        diff_num = 0;
        int leng = min(complete_compare_require, overlap_len);
        int i = 0;
        for (; i + 32 <= leng; i += 32) {

            __m256i t1 = _mm256_loadu_si256((__m256i const *)(str1 + offset + i));
            __m256i t2 = _mm256_loadu_si256((__m256i const *)(str2 + i));
            __m256i res = _mm256_cmpeq_epi8(t1, t2);
            unsigned mask = _mm256_movemask_epi8(res);
            diff_num += 32 - _mm_popcnt_u32(mask);
            if (diff_num > overlap_diff_limit) break;
        }
        for (; i < leng; i++) {
            diff_num += str1[offset + i] != str2[i];
            if (diff_num > overlap_diff_limit) break;
        }

        if (diff_num <= overlap_diff_limit) {
            delete[] tmpRev;
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
        overlap_len = min(len1, len2 - abs(offset));

        diff_num = 0;
        int leng = min(complete_compare_require, overlap_len);


        int i = 0;
        for (; i + 32 <= leng; i += 32) {
            __m256i t1 = _mm256_loadu_si256((__m256i const *)(str1 + i));
            __m256i t2 = _mm256_loadu_si256((__m256i const *)(str2 - offset + i));
            __m256i res = _mm256_cmpeq_epi8(t1, t2);
            unsigned mask = _mm256_movemask_epi8(res);
            diff_num += 32 - _mm_popcnt_u32(mask);
            if (diff_num > overlap_diff_limit) break;
        }
        for (; i < leng; i++) {
            diff_num += str1[i] != str2[-offset + i];
            if (diff_num > overlap_diff_limit) break;
        }
        if (diff_num <= overlap_diff_limit) {
            delete[] tmpRev;
            return {true, offset, overlap_len, diff_num};
        }
        offset -= 1;
    }
    delete[] tmpRev;
    return {false, 0, 0, 0};
}

#else

OverlapRes
Adapter::AnalyzeOverlap(neoReference &r1, neoReference &r2, int overlap_diff_limit, int overlap_require) {


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
        overlap_len = min(len1 - offset, len2);
        int leng = min(complete_compare_require, overlap_len);
        diff_num = 0;
        for (int i = 0; i < leng; i++) {
            diff_num += str1[offset + i] != rabbit::reMap[str2[len2 - 1 - i]];
            if (diff_num > overlap_diff_limit)
                break;
        }

        if (diff_num <= overlap_diff_limit) {
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
        overlap_len = min(len1, len2 - abs(offset));
        int leng = min(complete_compare_require, overlap_len);
        diff_num = 0;
        for (int i = 0; i < leng; i++) {
            diff_num += str1[i] != rabbit::reMap[str2[len2 - 1 + offset - i]];
            if (diff_num > overlap_diff_limit)
                break;
        }

        if (diff_num <= overlap_diff_limit) {
            return {true, offset, overlap_len, diff_num};
        }

        offset -= 1;
    }

    return {false, 0, 0, 0};
}

#endif

/**
 * @brief Trim adapter by giving adapter sequence
 * @param ref
 * @param adapter_seq
 * @param isR2
 * @return
 */

int Adapter::TrimAdapter(neoReference &ref, string &adapter_seq, bool isR2, int matchReq) {
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


#ifdef Vec512
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        int l = max(0, -pos);
        for (int i = l; i < cmplen; i += 64) {
            uint64_t tag = (1ll << min(64, cmplen - i)) - 1;
            __m512i t1 = _mm512_maskz_loadu_epi8(tag, adata + i);
            __m512i t2 = _mm512_maskz_loadu_epi8(tag, rdata + i + pos);
            __mmask64 res = _mm512_cmp_epi8_mask(t1, t2, 4);
            mismatch += _mm_popcnt_u64(res);
            if (mismatch > allowedMismatch) break;
        }


        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }
#elif Vec256
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        int l = max(0, -pos);
        int i = l;
        for (; i + 32 < cmplen; i += 32) {
            __m256i t1 = _mm256_loadu_si256((__m256i const *)(adata + i));
            __m256i t2 = _mm256_loadu_si256((__m256i const *)(rdata + i + pos));
            __m256i res = _mm256_cmpeq_epi8(t1, t2);
            unsigned mask = _mm256_movemask_epi8(res);
            mismatch += 32 - _mm_popcnt_u32(mask);
            if (mismatch > allowedMismatch) break;
        }
        for (; i < cmplen; i++) {
            mismatch += adata[i] != rdata[i + pos];
            if (mismatch > allowedMismatch) break;
        }
        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }
#elif KK

    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int ok = 1;
        for (int i = max(0, -pos), j = 0; i < cmplen && j < 3; i++, j++) {
            if (adata[i] != rdata[i + pos]) {
                ok = 0;
                break;
            }
        }
        if (ok == 0) {
            for (int i = cmplen - 1, j = 0; i >= max(0, -pos) && j < 3; i--, j++) {
                if (adata[i] != rdata[i + pos]) {
                    ok = 0;
                    break;
                }
            }
        }
        if (ok == 0) {
            for (int i = max(0, -pos) + 3, j = 0; i < cmplen && j < 3; i++, j++) {
                if (adata[i] != rdata[i + pos]) {
                    ok = 0;
                    break;
                }
            }
        }
        if (ok) {
            int mismatch = 0;
            for (int i = max(0, -pos); i < cmplen; i++) {
                mismatch += adata[i] != rdata[i + pos];
                if (mismatch > allowedMismatch) break;
            }
            if (mismatch <= allowedMismatch) {
                found = true;
                break;
            }
        }
    }
#else
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for (int i = max(0, -pos); i < cmplen; i++) {
            mismatch += adata[i] != rdata[i + pos];
            if (mismatch > allowedMismatch) break;
        }
        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }
#endif

    if (found) {
        int res_len = 0;
        if (pos < 0) {
            string adapter = adapter_seq.substr(0, alen + pos);
            ref.lseq = 0;
            ref.lqual = 0;
            res_len = adapter.length();

        } else {
            string new_adapter_seq = string(reinterpret_cast<const char *>(ref.base + ref.pseq + pos),
                    rlen - pos);
            ref.lseq = pos;
            ref.lqual = pos;
            res_len = new_adapter_seq.length();
        }
        return res_len;
    }

    return false;
}


int Adapter::TrimAdapter(neoReference &ref, string &adapter_seq, unordered_map<string, int> &mp, int adapter_len_lim, bool isR2, int matchReq) {
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


#ifdef Vec512
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        int l = max(0, -pos);
        for (int i = l; i < cmplen; i += 64) {
            uint64_t tag = (1ll << min(64, cmplen - i)) - 1;
            __m512i t1 = _mm512_maskz_loadu_epi8(tag, adata + i);
            __m512i t2 = _mm512_maskz_loadu_epi8(tag, rdata + i + pos);
            __mmask64 res = _mm512_cmp_epi8_mask(t1, t2, 4);
            mismatch += _mm_popcnt_u64(res);
            if (mismatch > allowedMismatch) break;
        }


        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }
#elif Vec256
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        int l = max(0, -pos);
        int i = l;
        for (; i + 32 < cmplen; i += 32) {
            __m256i t1 = _mm256_loadu_si256((__m256i const *)(adata + i));
            __m256i t2 = _mm256_loadu_si256((__m256i const *)(rdata + i + pos));
            __m256i res = _mm256_cmpeq_epi8(t1, t2);
            unsigned mask = _mm256_movemask_epi8(res);
            mismatch += 32 - _mm_popcnt_u32(mask);
            if (mismatch > allowedMismatch) break;
        }
        for (; i < cmplen; i++) {
            mismatch += adata[i] != rdata[i + pos];
            if (mismatch > allowedMismatch) break;
        }
        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }
#elif KK

    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int ok = 1;
        for (int i = max(0, -pos), j = 0; i < cmplen && j < 3; i++, j++) {
            if (adata[i] != rdata[i + pos]) {
                ok = 0;
                break;
            }
        }
        if (ok == 0) {
            for (int i = cmplen - 1, j = 0; i >= max(0, -pos) && j < 3; i--, j++) {
                if (adata[i] != rdata[i + pos]) {
                    ok = 0;
                    break;
                }
            }
        }
        if (ok == 0) {
            for (int i = max(0, -pos) + 3, j = 0; i < cmplen && j < 3; i++, j++) {
                if (adata[i] != rdata[i + pos]) {
                    ok = 0;
                    break;
                }
            }
        }
        if (ok) {
            int mismatch = 0;
            for (int i = max(0, -pos); i < cmplen; i++) {
                mismatch += adata[i] != rdata[i + pos];
                if (mismatch > allowedMismatch) break;
            }
            if (mismatch <= allowedMismatch) {
                found = true;
                break;
            }
        }
    }
#else
    for (pos = start; pos < rlen - matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for (int i = max(0, -pos); i < cmplen; i++) {
            mismatch += adata[i] != rdata[i + pos];
            if (mismatch > allowedMismatch) break;
        }
        if (mismatch <= allowedMismatch) {
            found = true;
            break;
        }
    }
#endif


    if (found) {
        int res_len = 0;
        if (pos < 0) {
            string adapter = adapter_seq.substr(0, alen + pos);
            ref.lseq = 0;
            ref.lqual = 0;
            res_len = adapter.length();
            if (adapter.length() >= adapter_len_lim) {
                if (mp.count(adapter) == 0) {
                    mp[adapter] = 1;
                } else {
                    mp[adapter]++;
                }
            }

        } else {
            string new_adapter_seq = string(reinterpret_cast<const char *>(ref.base + ref.pseq + pos),
                    rlen - pos);
            ref.lseq = pos;
            ref.lqual = pos;
            res_len = new_adapter_seq.length();

            if (new_adapter_seq.length() >= adapter_len_lim) {
                if (mp.count(new_adapter_seq) == 0) {
                    mp[new_adapter_seq] = 1;
                } else {
                    mp[new_adapter_seq]++;
                }
            }
        }
        return res_len;
    }
    return false;
}


int Adapter::TrimAdapters(neoReference &ref, vector<string> &adapter_seqs, bool isR2) {
    int matchs = 4;
    if(adapter_seqs.size() > 16)
        matchs = 5;
    if(adapter_seqs.size() > 256)
        matchs = 6;
    int res = 0;
    char *bases0 = reinterpret_cast<char *>(ref.base + ref.pseq); 
    int len0 = ref.lseq;
    for(int i = 0; i < adapter_seqs.size(); i++) {
        res += TrimAdapter(ref, adapter_seqs[i], isR2, matchs);
    }
    if(res) {
        string seq0 = string(bases0, len0);
        string new_adapter_seq = seq0.substr(ref.lseq, len0 - ref.lseq);
        res = new_adapter_seq.length();
    }
    return res;
}

int Adapter::TrimAdapters(neoReference &ref, vector<string> &adapter_seqs, unordered_map<string, int> &mp, int adapter_len_lim, bool isR2) {
    int matchs = 4;
    if(adapter_seqs.size() > 16)
        matchs = 5;
    if(adapter_seqs.size() > 256)
        matchs = 6;
    int res = 0;
    char *bases0 = reinterpret_cast<char *>(ref.base + ref.pseq); 
    int len0 = ref.lseq;
    for(int i = 0; i < adapter_seqs.size(); i++) {
        res += TrimAdapter(ref, adapter_seqs[i], isR2, matchs);
    }
    if(res) {
        string seq0 = string(bases0, len0);
        string new_adapter_seq = seq0.substr(ref.lseq, len0 - ref.lseq);
        res = new_adapter_seq.length();
        if (new_adapter_seq.length() >= adapter_len_lim) {
            if (mp.count(new_adapter_seq) == 0) {
                mp[new_adapter_seq] = 1;
            } else {
                mp[new_adapter_seq]++;
            }
        }
    }
    return res;
}

/**
 * @brief Trim adapter by overlap information
 * @param r1
 * @param r2
 * @param offset
 * @param overlap_len
 * @return true if successfully trim adapter, false else
 */
int Adapter::TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len) {
    //    if(ov.diff<=5 && ov.overlapped && ov.offset < 0 && ol > r1->length()/3)

    if (overlap_len > 0 && offset < 0) {
        string adapter1 = string(reinterpret_cast<const char *>(r1.base + r1.pseq + overlap_len),
                r1.lseq - overlap_len);
        string adapter2 = string(reinterpret_cast<const char *>(r2.base + r2.pseq + overlap_len),
                r2.lseq - overlap_len);
        r1.lseq = overlap_len;
        r1.lqual = overlap_len;
        r2.lseq = overlap_len;
        r2.lqual = overlap_len;

        return adapter1.length() + adapter2.length();
    }
    return false;
}

int Adapter::TrimAdapter(neoReference &r1, neoReference &r2, int offset, int overlap_len,
        unordered_map<string, int> &mp1, unordered_map<string, int> &mp2,
        int adapter_len_lim) {
    //    if(ov.diff<=5 && ov.overlapped && ov.offset < 0 && ol > r1->length()/3)

    if (overlap_len > 0 && offset < 0) {
        string adapter1 = string(reinterpret_cast<const char *>(r1.base + r1.pseq + overlap_len),
                r1.lseq - overlap_len);
        string adapter2 = string(reinterpret_cast<const char *>(r2.base + r2.pseq + overlap_len),
                r2.lseq - overlap_len);
        r1.lseq = overlap_len;
        r1.lqual = overlap_len;
        r2.lseq = overlap_len;
        r2.lqual = overlap_len;
        if (adapter1.length() >= adapter_len_lim) {
            if (mp1.count(adapter1) == 0) {
                mp1[adapter1] = 1;
            } else {
                mp1[adapter1]++;
            }
        }
        if (adapter2.length() >= adapter_len_lim) {
            if (mp2.count(adapter2) == 0) {
                mp2[adapter2] = 1;
            } else {
                mp2[adapter2]++;
            }
        }
        return adapter1.length() + adapter2.length();
    }
    return false;
}


int Adapter::CorrectData(neoReference &r1, neoReference &r2, OverlapRes &overlap_res, bool isPhred64) {
    if (!overlap_res.overlapped)
        return 0;

    int ol = overlap_res.overlap_len;

    int start1 = max(0, overlap_res.offset);
    int start2 = r2.lseq - max(0, -overlap_res.offset) - 1;

    char *seq1 = reinterpret_cast<char *>(r1.base + r1.pseq);
    char *seq2 = reinterpret_cast<char *>(r2.base + r2.pseq);
    char *qual1 = reinterpret_cast<char *>(r1.base + r1.pqual);
    char *qual2 = reinterpret_cast<char *>(r2.base + r2.pqual);

    int phredSub = 33;
    if (isPhred64) phredSub = 64;

    const char GOOD_QUAL = 30 + phredSub;//30
    const char BAD_QUAL = 14 + phredSub; //14



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
            } else if (qual2[p2] >= GOOD_QUAL && qual1[p1] <= BAD_QUAL) {
                // use R2

                seq1[p1] = rabbit::reMap[seq2[p2]];
                qual1[p1] = qual2[p2];
                corrected++;
                r1Corrected = true;
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


    return corrected;
}
