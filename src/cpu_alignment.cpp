#include <cpu_alignment.hpp>
#include <config.hpp>

#include <sequencehelpers.hpp>
#include <hostdevicefunctions.cuh>

#include <vector>
#include <cstring>

#ifdef Vec512
#include <immintrin.h>
#include <emmintrin.h>
#endif

namespace care{
    namespace cpu{
        namespace shd{

#ifdef Vec512
            void Print(__m128i a) {
                unsigned int c[4];
                memcpy(c, &a, sizeof(c));
                printf("%u %u %u %u\n", c[0], c[1], c[2], c[3]);
            }
            inline __m128i mm_bitshift_left1(__m128i x) {
                __m128i carry = _mm_bslli_si128(x, 8);
                carry = _mm_srli_epi64(carry, 63);
                x = _mm_slli_epi64(x, 1);
                return _mm_or_si128(x, carry);
            }
#endif
            AlignmentResult
                cpuShiftedHammingDistancePopcount2BitHiLo(
                        CpuAlignmentHandle& handle,
                        const unsigned int* anchorHiLo,
                        int anchorLength,
                        const unsigned int* candidateHiLo,
                        int candidateLength,
                        int min_overlap,
                        float maxErrorRate,
                        float min_overlap_ratio) noexcept{

                    assert(anchorLength > 0);
                    assert(candidateLength > 0);

                    int cnt1 = 0;
                    int cnt2 = 0;
                    auto popcount = [](auto i){return __builtin_popcount(i);};
                    auto identity = [](auto i){return i;};




                    auto hammingDistanceWithShiftNum = [&](int shift_len, int overlapsize, int max_errors,
                            unsigned int* shiftptr_hi, unsigned int* shiftptr_lo, auto transfunc1,
                            int shiftptr_size,
                            const unsigned int* otherptr_hi, const unsigned int* otherptr_lo,
                            auto transfunc2){

                        if(shift_len < 0) shift_len = -shift_len;

                        {
                            const unsigned int* lhi = shiftptr_hi;
                            const unsigned int* llo = shiftptr_lo;
                            const unsigned int* rhi = otherptr_hi;
                            const unsigned int* rlo = otherptr_lo;
                            int lhi_bitcount = overlapsize;
                            int rhi_bitcount = overlapsize;

                            const int overlap_bitcount = std::min(lhi_bitcount, rhi_bitcount);

                            if(overlap_bitcount == 0)
                                return max_errors+1;

                            //const int partitions = SDIV(overlap_bitcount, (8 * sizeof(unsigned int)));

                            const int partitions = (overlap_bitcount + 32 - 1) >> 5;
                            const int remaining_bitcount = (partitions << 5) - overlap_bitcount;
                            int pos_shift = shift_len >> 5; // 0
                            int bit_shift = shift_len & ((1 << 5) - 1); // 1
                            int result = 0;
                            for(int i = 0; i < partitions - 1 && result < max_errors; i += 1) {
                                unsigned int A1 = lhi[i + pos_shift] << bit_shift;
                                unsigned int A2 = llo[i + pos_shift] << bit_shift;
                                unsigned int B1 = lhi[i + pos_shift + 1] >> (32 - bit_shift);
                                unsigned int B2 = llo[i + pos_shift + 1] >> (32 - bit_shift);
                                if(bit_shift == 0) {
                                    B1 = 0;
                                    B2 = 0;
                                }
                                const unsigned int hixor = (A1 | B1) ^ rhi[i];
                                const unsigned int loxor = (A2 | B2) ^ rlo[i];
                                const unsigned int bits = hixor | loxor;
                                result += popcount(bits);

                            }
                            if(result >= max_errors)
                                return result;

                            const unsigned int mask = remaining_bitcount == 0 ? 0xFFFFFFFF : 0xFFFFFFFF << (remaining_bitcount);
                            unsigned int A1 = lhi[partitions - 1 + pos_shift] << bit_shift;
                            unsigned int A2 = llo[partitions - 1 + pos_shift] << bit_shift;
                            unsigned int B1 = lhi[partitions - 1 + pos_shift + 1] >> (32 - bit_shift);
                            unsigned int B2 = llo[partitions - 1 + pos_shift + 1] >> (32 - bit_shift);
                            if(bit_shift == 0) {
                                B1 = 0;
                                B2 = 0;
                            }
                            const unsigned int hixor = (A1 | B1) ^ rhi[partitions - 1];
                            const unsigned int loxor = (A2 | B2) ^ rlo[partitions - 1];
                            const unsigned int bits = hixor | loxor;
                            result += popcount(bits & mask);
                            return result;

                        }
                    };


                    auto hammingDistanceWithShift = [&](bool doShift, int overlapsize, int max_errors,
                            unsigned int* shiftptr_hi, unsigned int* shiftptr_lo, auto transfunc1,
                            int shiftptr_size,
                            const unsigned int* otherptr_hi, const unsigned int* otherptr_lo,
                            auto transfunc2){
                        if(doShift){
                            //shiftBitArrayLeftBy<1>(shiftptr_hi, shiftptr_size / 2, transfunc1);
                            {
                                unsigned int* array = shiftptr_hi;
                                int size = shiftptr_size / 2;

                                for(int i = 0; i < size - 1; i += 1) {
                                    const unsigned int a = array[i];
                                    const unsigned int b = array[i + 1];
                                    array[i] = (a << 1) | (b >> 31);
                                }

                                array[size - 1] <<= 1;
                            }
                            //shiftBitArrayLeftBy<1>(shiftptr_lo, shiftptr_size / 2, transfunc1);
                            {
                                unsigned int* array = shiftptr_lo;
                                int size = shiftptr_size / 2;

                                for(int i = 0; i < size - 1; i += 1) {
                                    const unsigned int a = array[i];
                                    const unsigned int b = array[i + 1];
                                    array[i] = (a << 1) | (b >> 31);
                                }

                                array[size - 1] <<= 1;
                            }
                        }
                        //const int score = hammingdistanceHiLo(shiftptr_hi,
                        //		shiftptr_lo,
                        //		otherptr_hi,
                        //		otherptr_lo,
                        //		overlapsize,
                        //		overlapsize,
                        //		max_errors,
                        //		transfunc1,
                        //		transfunc2,
                        //		popcount);
                        //return score;

                        {
                            const unsigned int* lhi = shiftptr_hi;
                            const unsigned int* llo = shiftptr_lo;
                            const unsigned int* rhi = otherptr_hi;
                            const unsigned int* rlo = otherptr_lo;
                            int lhi_bitcount = overlapsize;
                            int rhi_bitcount = overlapsize;

                            const int overlap_bitcount = std::min(lhi_bitcount, rhi_bitcount);

                            if(overlap_bitcount == 0)
                                return max_errors+1;

                            //const int partitions = SDIV(overlap_bitcount, (8 * sizeof(unsigned int)));

                            const int partitions = (overlap_bitcount + 32 - 1) >> 5;
                            const int remaining_bitcount = (partitions << 5) - overlap_bitcount;

                            int result = 0;
                            for(int i = 0; i < partitions - 1 && result < max_errors; i += 1) {
                                const unsigned int hixor = lhi[i] ^ rhi[i];
                                const unsigned int loxor = llo[i] ^ rlo[i];
                                const unsigned int bits = hixor | loxor;
                                result += popcount(bits);
                            }
                            if(result >= max_errors)
                                return result;

                            // i == partitions - 1

                            const unsigned int mask = remaining_bitcount == 0 ? 0xFFFFFFFF : 0xFFFFFFFF << (remaining_bitcount);
                            const unsigned int hixor = lhi[partitions - 1] ^ rhi[partitions - 1];
                            const unsigned int loxor = llo[partitions - 1] ^ rlo[partitions - 1];
                            const unsigned int bits = hixor | loxor;
                            result += popcount(bits & mask);

                            return result;
                        }
                    };

                    auto& shiftbuffer = handle.shiftbuffer;

                    const int anchorInts = SequenceHelpers::getEncodedNumInts2BitHiLo(anchorLength);
                    const int candidateInts = SequenceHelpers::getEncodedNumInts2BitHiLo(candidateLength);
                    //const int maxInts = std::max(anchorInts, candidateInts);
                    const int maxInts = (anchorInts + candidateInts) * 2;


                    shiftbuffer.resize(maxInts);
                    for(int i = 0; i < maxInts; i++) {
                        shiftbuffer[i] = 0;
                    }

                    const unsigned int* const anchorBackup_hi = anchorHiLo;
                    const unsigned int* const anchorBackup_lo = anchorHiLo + anchorInts / 2;

                    const unsigned int* const candidateBackup_hi = candidateHiLo;
                    const unsigned int* const candidateBackup_lo = candidateHiLo + candidateInts / 2;

                    const int totalbases = anchorLength + candidateLength;
                    const int minoverlap = std::max(min_overlap, int(float(anchorLength) * min_overlap_ratio));

                    int bestScore = totalbases; // score is number of mismatches
                    int bestShift = -candidateLength; // shift of query relative to anchor. shift < 0 if query begins before anchor

                    auto handle_shift = [&](int shift, int overlapsize,
                            unsigned int* shiftptr_hi, unsigned int* shiftptr_lo, auto transfunc1,
                            int shiftptr_size,
                            const unsigned int* otherptr_hi, const unsigned int* otherptr_lo,
                            auto transfunc2){

                        const int max_errors_excl = std::min(int(float(overlapsize) * maxErrorRate),
                                bestScore - totalbases + 2*overlapsize);

                        if(max_errors_excl > 0){

                            //const int mismatches = hammingDistanceWithShift(shift != 0, overlapsize, max_errors_excl,
                            //const int mismatches = hammingDistanceWithShiftNum(shift, overlapsize, max_errors_excl,
                            //     shiftptr_hi,shiftptr_lo, transfunc1,
                            //     shiftptr_size,
                            //     otherptr_hi, otherptr_lo, transfunc2);

                            int mismatches = 0;
                            int shift_len = shift;
                            if(shift_len < 0) shift_len = -shift_len;

                            const unsigned int* lhi = shiftptr_hi;
                            const unsigned int* llo = shiftptr_lo;
                            const unsigned int* rhi = otherptr_hi;
                            const unsigned int* rlo = otherptr_lo;
                            int max_errors = max_errors_excl;
                            int lhi_bitcount = overlapsize;
                            int rhi_bitcount = overlapsize;

                            const int overlap_bitcount = std::min(lhi_bitcount, rhi_bitcount);

                            if(overlap_bitcount == 0) {
                                mismatches = max_errors + 1;
                            } else {
                                const int partitions = (overlap_bitcount + 32 - 1) >> 5;
                                const int remaining_bitcount = (partitions << 5) - overlap_bitcount;
                                int pos_shift = shift_len >> 5; // 0
                                int bit_shift = shift_len & ((1 << 5) - 1); // 1
                                int result = 0;
                                for(int i = 0; i < partitions - 1 && result < max_errors; i += 1) {
                                    unsigned int A1 = lhi[i + pos_shift] << bit_shift;
                                    unsigned int A2 = llo[i + pos_shift] << bit_shift;
                                    unsigned int B1 = lhi[i + pos_shift + 1] >> (32 - bit_shift);
                                    unsigned int B2 = llo[i + pos_shift + 1] >> (32 - bit_shift);
                                    if(bit_shift == 0) {
                                        B1 = 0;
                                        B2 = 0;
                                    }
                                    const unsigned int hixor = (A1 | B1) ^ rhi[i];
                                    const unsigned int loxor = (A2 | B2) ^ rlo[i];
                                    const unsigned int bits = hixor | loxor;
                                    result += popcount(bits);

                                }
                                if(result < max_errors) {
                                    const unsigned int mask = remaining_bitcount == 0 ? 0xFFFFFFFF : 0xFFFFFFFF << (remaining_bitcount);
                                    unsigned int A1 = lhi[partitions - 1 + pos_shift] << bit_shift;
                                    unsigned int A2 = llo[partitions - 1 + pos_shift] << bit_shift;
                                    unsigned int B1 = lhi[partitions - 1 + pos_shift + 1] >> (32 - bit_shift);
                                    unsigned int B2 = llo[partitions - 1 + pos_shift + 1] >> (32 - bit_shift);
                                    if(bit_shift == 0) {
                                        B1 = 0;
                                        B2 = 0;
                                    }
                                    const unsigned int hixor = (A1 | B1) ^ rhi[partitions - 1];
                                    const unsigned int loxor = (A2 | B2) ^ rlo[partitions - 1];
                                    const unsigned int bits = hixor | loxor;
                                    result += popcount(bits & mask);
                                }
                                mismatches = result;
                            }

                            const int score = (mismatches < max_errors_excl ?
                                    mismatches + totalbases - 2 * overlapsize // non-overlapping regions count as mismatches
                                    : std::numeric_limits<int>::max()); // too many errors, discard

                            if(score < bestScore){
                                bestScore = score;
                                bestShift = shift;
                            }

                            return true;
                        }else{
                            return false;
                        }
                    };


                    std::copy_n(anchorHiLo, anchorInts / 2, shiftbuffer.begin());
                    std::copy_n(anchorHiLo + anchorInts / 2, anchorInts / 2, shiftbuffer.begin() + maxInts / 2);
                    unsigned int* shiftbuffer_hi = shiftbuffer.data();
                    unsigned int* shiftbuffer_lo = shiftbuffer.data() + maxInts / 2;



                    for(int shift = 0; shift < anchorLength - minoverlap + 1; ++shift){
                        const int overlapsize = std::min(anchorLength - shift, candidateLength);


                       	//bool b = handle_shift(shift, overlapsize,
                        //        shiftbuffer_hi, shiftbuffer_lo, identity,
                        //        anchorInts,
                        //        candidateBackup_hi, candidateBackup_lo, identity);
                        //if(!b){
                        //    break;
                        //} 


                        const int max_errors = std::min(int(float(overlapsize) * maxErrorRate), bestScore - totalbases + 2 * overlapsize);
                        if(max_errors <= 0) break;
                        int mismatches = 0;
                        int shift_len = shift;
                        if(shift_len < 0) shift_len = -shift_len;
                        const unsigned int* lhi = shiftbuffer_hi;
                        const unsigned int* llo = shiftbuffer_lo;
                        const unsigned int* rhi = candidateBackup_hi;
                        const unsigned int* rlo = candidateBackup_lo;
                        if(overlapsize == 0) {
                            mismatches = max_errors + 1;
                        } else {
                            const int partitions = (overlapsize + 32 - 1) >> 5;
                            const int remaining_bitcount = (partitions << 5) - overlapsize;
                            int pos_shift = shift_len >> 5; // 0
                            int bit_shift = shift_len & ((1 << 5) - 1); // 1
                            int result = 0;
                            for(int i = 0; i < partitions - 1 && result < max_errors; i += 1) {
                                unsigned int A1 = lhi[i + pos_shift] << bit_shift;
                                unsigned int A2 = llo[i + pos_shift] << bit_shift;

 
                                unsigned int B1 = lhi[i + pos_shift + 1] >> (32 - bit_shift);
                                unsigned int B2 = llo[i + pos_shift + 1] >> (32 - bit_shift);
                                if(bit_shift == 0) {
                                    B1 = 0;
                                    B2 = 0;
                                }
                                const unsigned int hixor = (A1 | B1) ^ rhi[i];
                                const unsigned int loxor = (A2 | B2) ^ rlo[i];
                                const unsigned int bits = hixor | loxor;
                                result += popcount(bits);
                            }
                            if(result < max_errors) {
                                const unsigned int mask = remaining_bitcount == 0 ? 0xFFFFFFFF : 0xFFFFFFFF << (remaining_bitcount);
                                unsigned int A1 = lhi[partitions - 1 + pos_shift] << bit_shift;
                                unsigned int A2 = llo[partitions - 1 + pos_shift] << bit_shift;
                                unsigned int B1 = lhi[partitions - 1 + pos_shift + 1] >> (32 - bit_shift);
                                unsigned int B2 = llo[partitions - 1 + pos_shift + 1] >> (32 - bit_shift);
                                if(bit_shift == 0) {
                                    B1 = 0;
                                    B2 = 0;
                                }
                                const unsigned int hixor = (A1 | B1) ^ rhi[partitions - 1];
                                const unsigned int loxor = (A2 | B2) ^ rlo[partitions - 1];
                                const unsigned int bits = hixor | loxor;
                                result += popcount(bits & mask);
                            }

                            mismatches = result;
                        }
                        const int score = (mismatches < max_errors ?
                                mismatches + totalbases - 2 * overlapsize // non-overlapping regions count as mismatches
                                : std::numeric_limits<int>::max()); // too many errors, discard

                        if(score < bestScore){
                            bestScore = score;
                            bestShift = shift;
                        }

                    }


                    std::copy_n(candidateHiLo, candidateInts / 2, shiftbuffer.begin());
                    std::copy_n(candidateHiLo + candidateInts / 2, candidateInts / 2, shiftbuffer.begin() + maxInts / 2);
                    shiftbuffer_hi = shiftbuffer.data();
                    shiftbuffer_lo = shiftbuffer.data() + maxInts / 2;




                    for(int shift = -1; shift >= -candidateLength + minoverlap; --shift){
                        const int overlapsize = std::min(anchorLength, candidateLength + shift);


						//bool b = handle_shift(shift, overlapsize,
                        //        shiftbuffer_hi, shiftbuffer_lo, identity,
                        //        candidateInts,
                        //        anchorBackup_hi, anchorBackup_lo, identity);

                        //if(!b){
                        //    break;
                        //}


                        const int max_errors = std::min(int(float(overlapsize) * maxErrorRate), bestScore - totalbases + 2 * overlapsize);
                        if(max_errors <= 0) break;
                        int mismatches = 0;
                        int shift_len = shift;
                        if(shift_len < 0) shift_len = -shift_len;
                        const unsigned int* lhi = shiftbuffer_hi;
                        const unsigned int* llo = shiftbuffer_lo;
                        const unsigned int* rhi = anchorBackup_hi;
                        const unsigned int* rlo = anchorBackup_lo;
                        if(overlapsize == 0) {
                            mismatches = max_errors + 1;
                        } else {
                            const int partitions = (overlapsize + 32 - 1) >> 5;
                            const int remaining_bitcount = (partitions << 5) - overlapsize;
                            int pos_shift = shift_len >> 5; // 0
                            int bit_shift = shift_len & ((1 << 5) - 1); // 1
                            int result = 0;
                            for(int i = 0; i < partitions - 1 && result < max_errors; i += 1) {
                                unsigned int A1 = lhi[i + pos_shift] << bit_shift;
                                unsigned int A2 = llo[i + pos_shift] << bit_shift;
                                unsigned int B1 = lhi[i + pos_shift + 1] >> (32 - bit_shift);
                                unsigned int B2 = llo[i + pos_shift + 1] >> (32 - bit_shift);
                                if(bit_shift == 0) {
                                    B1 = 0;
                                    B2 = 0;
                                }
                                const unsigned int hixor = (A1 | B1) ^ rhi[i];
                                const unsigned int loxor = (A2 | B2) ^ rlo[i];
                                const unsigned int bits = hixor | loxor;
                                result += popcount(bits);

                            }
                            if(result < max_errors) {
                                const unsigned int mask = remaining_bitcount == 0 ? 0xFFFFFFFF : 0xFFFFFFFF << (remaining_bitcount);
                                unsigned int A1 = lhi[partitions - 1 + pos_shift] << bit_shift;
                                unsigned int A2 = llo[partitions - 1 + pos_shift] << bit_shift;
                                unsigned int B1 = lhi[partitions - 1 + pos_shift + 1] >> (32 - bit_shift);
                                unsigned int B2 = llo[partitions - 1 + pos_shift + 1] >> (32 - bit_shift);
                                if(bit_shift == 0) {
                                    B1 = 0;
                                    B2 = 0;
                                }
                                const unsigned int hixor = (A1 | B1) ^ rhi[partitions - 1];
                                const unsigned int loxor = (A2 | B2) ^ rlo[partitions - 1];
                                const unsigned int bits = hixor | loxor;
                                result += popcount(bits & mask);
                            }


                            mismatches = result;
                        }
                        const int score = (mismatches < max_errors ?
                                mismatches + totalbases - 2 * overlapsize // non-overlapping regions count as mismatches
                                : std::numeric_limits<int>::max()); // too many errors, discard

                        if(score < bestScore){
                            bestScore = score;
                            bestShift = shift;
                        }
                    }

                    AlignmentResult alignmentresult;
                    alignmentresult.isValid = (bestShift != -candidateLength);

                    const int candidateoverlapbegin_incl = std::max(-bestShift, 0);
                    const int candidateoverlapend_excl = std::min(candidateLength, anchorLength - bestShift);
                    const int overlapsize = candidateoverlapend_excl - candidateoverlapbegin_incl;
                    const int opnr = bestScore - totalbases + 2*overlapsize;

                    alignmentresult.score = bestScore;
                    alignmentresult.overlap = overlapsize;
                    alignmentresult.shift = bestShift;
                    alignmentresult.nOps = opnr;

                    return alignmentresult;
                }






        }
    }
}
