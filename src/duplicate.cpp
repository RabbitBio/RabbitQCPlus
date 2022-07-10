//
// Created by ylf9811 on 2021/7/13.
//

#include "duplicate.h"

#ifdef Vec512

#include <immintrin.h>

#endif

Duplicate::Duplicate(CmdInfo *cmd_info) {
    cmd_info_ = cmd_info;
    key_len_base_ = 12;
    key_len_bit = 1 << (2 * key_len_base_);
    dups_ = new uint64_t[key_len_bit];
    memset(dups_, 0, sizeof(uint64_t) * key_len_bit);
    counts_ = new uint16_t[key_len_bit];
    memset(counts_, 0, sizeof(uint16_t) * key_len_bit);
    gcs_ = new uint8_t[key_len_bit];
    memset(gcs_, 0, sizeof(uint8_t) * key_len_bit);
}

Duplicate::~Duplicate() {
    delete[] dups_;
    delete[] counts_;
    delete[] gcs_;
}


uint64_t Duplicate::seq2int(const char *data, int start, int key_len, bool &valid) {
    uint64_t ret = 0;
    for (int i = 0; i < key_len; i++) {
        ret <<= 2;
        if (valAGCT[data[start + i] & 0x07] == -1) {
            valid = false;
            return 0;
        }
        ret += valAGCT[data[start + i] & 0x07];
        // if it's not the last one, shift it by 2 bits
    }
    return ret;
}


void Duplicate::addRecord(uint32_t key, uint64_t kmer32, uint8_t gc) {
    //    lok.lock();
    //    printf("thread %d is duplicating ...\n", this_thread::get_id());
    //TODO what if kmer1 == kmer2 but gc1 != gc2 (of cause key1 == key2)
    //even if set lock in this function, it is stall thread unsafe.
    //now change code to make it thread safe, but maybe it can be case a logic error.
    if (counts_[key] == 0) {
        counts_[key] = 1;
        dups_[key] = kmer32;
        gcs_[key] = gc;
    } else {
        if (dups_[key] == kmer32) {
            counts_[key]++;
            //add this
            //TODO check it is still logic correct or not
            if (gcs_[key] > gc) gcs_[key] = gc;
        } else if (dups_[key] > kmer32) {
            dups_[key] = kmer32;
            counts_[key] = 1;
            gcs_[key] = gc;
        }
    }
    //    lok.unlock();
}

void Duplicate::statRead(neoReference &ref) {
    if (ref.lseq < 32)
        return;

    int start1 = 0;
    uint64_t zr = 0;
    int start2 = std::max(zr, ref.lseq - 32 - 5);

    const char *data = reinterpret_cast<const char *>(ref.base + ref.pseq);
    bool valid = true;

    uint64_t ret = seq2int(data, start1, key_len_base_, valid);
    uint32_t key = (uint32_t) ret;
    if (!valid)
        return;

    uint64_t kmer32 = seq2int(data, start2, 32, valid);
    if (!valid)
        return;

    int gc = 0;

    // not calculated
    //TODO check correctness
//    if (counts_[key] == 0) {
#ifdef Vec512
    int i = 0;
    __m512i gcV = _mm512_set1_epi32(0);
    __m512i ad1 = _mm512_set1_epi32(1);

    __m128i gcC = _mm_set1_epi8('C');
    __m128i gcT = _mm_set1_epi8('G');
    for (; i + 16 <= ref.lseq; i += 16) {
        __m128i ide = _mm_maskz_loadu_epi8(0xFFFF, data + i);
        __mmask16 mk1 = _mm_cmpeq_epi8_mask(ide, gcC);
        __mmask16 mk2 = _mm_cmpeq_epi8_mask(ide, gcT);
        mk1 = mk1 | mk2;
        gcV = _mm512_mask_add_epi32(gcV, mk1, gcV, ad1);
    }
    for (; i < ref.lseq; i++) {
        if (data[i] == 'C' || data[i] == 'G')
            gc++;
    }
    int cnt[16];
    _mm512_store_epi32(cnt, gcV);
    for (int k = 0; k < 16; k++) gc += cnt[k];
#else
    for (int i = 0; i < ref.lseq; i++) {
        if (data[i] == 'C' || data[i] == 'G')
            gc++;
    }
#endif
    //    }
    gc = int(255.0 * gc / ref.lseq + 0.5);

    addRecord(key, kmer32, (uint8_t) gc);
}

void Duplicate::statPair(neoReference &r1, neoReference &r2) {
    if (r1.lseq < 32 || r2.lseq < 32)
        return;

    const char *data1 = reinterpret_cast<const char *>(r1.base + r1.pseq);
    const char *data2 = reinterpret_cast<const char *>(r2.base + r2.pseq);
    bool valid = true;

    uint64_t ret = seq2int(data1, 0, key_len_base_, valid);
    uint32_t key = (uint32_t) ret;
    if (!valid)
        return;

    uint64_t kmer32 = seq2int(data2, 0, 32, valid);
    if (!valid)
        return;

    int gc = 0;

    // not calculated
    //TODO avx512
    if (counts_[key] == 0) {
        for (int i = 0; i < r1.lseq; i++) {
            if (data1[i] == 'G' || data1[i] == 'C')
                gc++;
        }
        for (int i = 0; i < r2.lseq; i++) {
            if (data2[i] == 'G' || data2[i] == 'C')
                gc++;
        }
    }

    gc = round(255.0 * (double) gc / (double) (r1.lseq + r2.lseq));
    addRecord(key, kmer32, gc);
}

double Duplicate::statAll(int *hist, double *meanGC, int histSize) {
    long totalNum = 0;
    long dupNum = 0;
    int *gcStatNum = new int[histSize];
    memset(gcStatNum, 0, sizeof(int) * histSize);
    for (int key = 0; key < key_len_bit; key++) {
        int count = counts_[key];
        double gc = gcs_[key];

        if (count > 0) {
            totalNum += count;
            dupNum += count - 1;

            if (count >= histSize) {
                hist[histSize - 1]++;
                meanGC[histSize - 1] += gc;
                gcStatNum[histSize - 1]++;
            } else {
                hist[count]++;
                meanGC[count] += gc;
                gcStatNum[count]++;
            }
        }
    }

    for (int i = 0; i < histSize; i++) {
        if (gcStatNum[i] > 0) {
            meanGC[i] = meanGC[i] / 255.0 / gcStatNum[i];
        }
    }

    delete[] gcStatNum;

    if (totalNum == 0)
        return 0.0;
    else
        return (double) dupNum / (double) totalNum;
}