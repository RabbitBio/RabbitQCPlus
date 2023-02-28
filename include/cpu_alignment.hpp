#ifndef CARE_CPU_ALIGNMENT_HPP
#define CARE_CPU_ALIGNMENT_HPP


#include <hostdevicefunctions.cuh>
#include <sequencehelpers.hpp>

#include <config.hpp>

#include <vector>
#include <cassert>
#include <algorithm>

namespace care{
namespace cpu{

    template<class Iter1, class Iter2, class Equal>
    int hammingDistance(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, Equal isEqual){
        int result = 0;

        while(first1 != last1 && first2 != last2){
            result += isEqual(*first1, *first2) ? 0 : 1;

            ++first1;
            ++first2;
        }

        //positions which do not overlap count as mismatch.
        //at least one of the remaining ranges is empty
        result += std::distance(first1, last1);
        result += std::distance(first2, last2);

        return result;
    }

    template<class Iter1, class Iter2>
    int hammingDistance(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2){
        auto isEqual = [](const auto& l, const auto& r){
            return l == r;
        };

        return hammingDistance(first1, last1, first2, last2, isEqual);
    }

    template<class Iter1, class Iter2, class Equal>
    int hammingDistanceOverlap(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, Equal isEqual){
        int result = 0;

        while(first1 != last1 && first2 != last2){
            result += isEqual(*first1, *first2) ? 0 : 1;

            ++first1;
            ++first2;
        }

        //positions which do not overlap are excluded

        return result;
    }

    template<class Iter1, class Iter2>
    int hammingDistanceOverlap(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2){
        auto isEqual = [](const auto& l, const auto& r){
            return l == r;
        };

        return hammingDistanceOverlap(first1, last1, first2, last2, isEqual);
    }



    template<class Iter1, class Iter2, class Equal>
    int longestMatch(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, Equal isEqual){
        int longest = 0;
        int current = 0;

        while(first1 != last1 && first2 != last2){
            if(isEqual(*first1, *first2)){
                current++;
            }else{
                longest = std::max(longest, current);
                current = 0;
            }

            ++first1;
            ++first2;
        }

        longest = std::max(longest, current);

        return longest;
    }

    template<class Iter1, class Iter2>
    int longestMatch(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2){
        auto isEqual = [](const auto& l, const auto& r){
            return l == r;
        };

        return longestMatch(first1, last1, first2, last2, isEqual);
    }





namespace shd{

    enum class ShiftDirection {Left, Right, LeftRight, None};

    struct AlignmentResult{
    	int score;
    	int overlap;
    	int shift;
    	int nOps; //edit distance / number of operations
    	bool isValid;

        int get_score() const { return score;}
        int get_overlap() const { return overlap;}
        int get_shift() const { return shift;}
        int get_nOps() const { return nOps;}
        bool get_isValid() const { return isValid;}

        bool operator==(const AlignmentResult& rhs) const {
            return score == rhs.score && overlap == rhs.overlap && shift == rhs.shift && nOps == rhs.nOps && isValid == rhs.isValid;
        }
        bool operator!=(const AlignmentResult& rhs) const{
            return !(operator==(rhs));
        }
    };

    struct CpuAlignmentHandle{
        std::vector<unsigned int> shiftbuffer;
        std::vector<unsigned int> anchorConversionBuffer;
        std::vector<unsigned int> candidateConversionBuffer;
    };

    AlignmentResult
    cpuShiftedHammingDistancePopcount2BitHiLo(
            CpuAlignmentHandle& handle,
            const unsigned int* anchorHiLo,
            int anchorLength,
            const unsigned int* candidateHiLo,
            int candidateLength,
            int min_overlap,
            float maxErrorRate,
            float min_overlap_ratio) noexcept;


    template<ShiftDirection direction>
    AlignmentResult
    cpuShiftedHammingDistancePopcount2BitHiLoWithDirection(
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

        auto popcount = [](auto i){return __builtin_popcount(i);};
        auto identity = [](auto i){return i;};

        auto hammingDistanceWithShift = [&](bool doShift, int overlapsize, int max_errors,
                                            unsigned int* shiftptr_hi, unsigned int* shiftptr_lo, auto transfunc1,
                                            int shiftptr_size,
                                            const unsigned int* otherptr_hi, const unsigned int* otherptr_lo,
                                            auto transfunc2){
            if(doShift){
                shiftBitArrayLeftBy<1>(shiftptr_hi, shiftptr_size / 2, transfunc1);
                shiftBitArrayLeftBy<1>(shiftptr_lo, shiftptr_size / 2, transfunc1);
            }
            const int score = hammingdistanceHiLo(shiftptr_hi,
                                                    shiftptr_lo,
                                                    otherptr_hi,
                                                    otherptr_lo,
                                                    overlapsize,
                                                    overlapsize,
                                                    max_errors,
                                                    transfunc1,
                                                    transfunc2,
                                                    popcount);

            return score;
        };

        auto& shiftbuffer = handle.shiftbuffer;

        const int anchorInts = SequenceHelpers::getEncodedNumInts2BitHiLo(anchorLength);
        const int candidateInts = SequenceHelpers::getEncodedNumInts2BitHiLo(candidateLength);
        const int maxInts = std::max(anchorInts, candidateInts);


        shiftbuffer.resize(maxInts);

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

                const int mismatches = hammingDistanceWithShift(shift != 0, overlapsize, max_errors_excl,
                                                        shiftptr_hi,shiftptr_lo, transfunc1,
                                                        shiftptr_size,
                                                        otherptr_hi, otherptr_lo, transfunc2);

                const int score = (mismatches < max_errors_excl ?
                                mismatches + totalbases - 2*overlapsize // non-overlapping regions count as mismatches
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

        auto handle_zero_shift = [&](int overlapsize,
                                const unsigned int* ptr_hi, const unsigned int* ptr_lo, auto transfunc1,
                                int /*ptr_size*/,
                                const unsigned int* otherptr_hi, const unsigned int* otherptr_lo,
                                auto transfunc2){

            const int shift = 0;

            const int max_errors_excl = std::min(int(float(overlapsize) * maxErrorRate),
                                            bestScore - totalbases + 2*overlapsize);

            if(max_errors_excl > 0){

                const int mismatches = hammingdistanceHiLo(ptr_hi,
                                                    ptr_lo,
                                                    otherptr_hi,
                                                    otherptr_lo,
                                                    overlapsize,
                                                    overlapsize,
                                                    max_errors_excl,
                                                    transfunc1,
                                                    transfunc2,
                                                    popcount);

                const int score = (mismatches < max_errors_excl ?
                                mismatches + totalbases - 2*overlapsize // non-overlapping regions count as mismatches
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

        //shift == 0
        {
            const int shift = 0;
            const int overlapsize = std::min(anchorLength - shift, candidateLength);

            handle_zero_shift(overlapsize,
                    anchorBackup_hi, anchorBackup_lo, identity,
                    anchorInts,
                    candidateBackup_hi, candidateBackup_lo, identity);

        }

        unsigned int* shiftbuffer_hi;
        unsigned int* shiftbuffer_lo;

        if(direction == ShiftDirection::Right || direction == ShiftDirection::LeftRight){
            //compute alignments with shift to the right
            std::copy_n(anchorHiLo, anchorInts, shiftbuffer.begin());
            shiftbuffer_hi = shiftbuffer.data();
            shiftbuffer_lo = shiftbuffer.data() + anchorInts / 2;

            for(int shift = 1; shift < anchorLength - minoverlap + 1; ++shift){
                const int overlapsize = std::min(anchorLength - shift, candidateLength);
                bool b = handle_shift(shift, overlapsize,
                                    shiftbuffer_hi, shiftbuffer_lo, identity,
                                    anchorInts,
                                    candidateBackup_hi, candidateBackup_lo, identity);
                if(!b){
                    break;
                }
            }
        }

        if(direction == ShiftDirection::Left || direction == ShiftDirection::LeftRight){

            std::copy_n(candidateHiLo, candidateInts, shiftbuffer.begin());
            shiftbuffer_hi = shiftbuffer.data();
            shiftbuffer_lo = shiftbuffer.data() + candidateInts / 2;

            for(int shift = -1; shift >= -candidateLength + minoverlap; --shift){
                const int overlapsize = std::min(anchorLength, candidateLength + shift);

                bool b = handle_shift(shift, overlapsize,
                                    shiftbuffer_hi, shiftbuffer_lo, identity,
                                    candidateInts,
                                    anchorBackup_hi, anchorBackup_lo, identity);

                if(!b){
                    break;
                }
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

    

    template<ShiftDirection direction, class Iter>
    Iter
    cpuShiftedHammingDistancePopcount2BitWithDirection(
            CpuAlignmentHandle& handle,
            Iter destinationBegin,
            const unsigned int* anchor2Bit,
            int anchorLength,
            const unsigned int* candidates2Bit,
            int candidatePitchInInts,
            const int* candidateLengths,
            int numCandidates,
            int min_overlap,
            float maxErrorRate,
            float min_overlap_ratio) noexcept{

        const int newanchorInts = SequenceHelpers::getEncodedNumInts2BitHiLo(anchorLength);

        handle.anchorConversionBuffer.resize(newanchorInts);        

        SequenceHelpers::convert2BitTo2BitHiLo(
            handle.anchorConversionBuffer.data(),
            anchor2Bit,
            anchorLength
        );

        auto curIter = destinationBegin;

        for(int candidateIndex = 0; candidateIndex < numCandidates; candidateIndex++, ++curIter){
            const unsigned int* candidate2Bit = candidates2Bit + candidatePitchInInts * candidateIndex;
            const int candidateLength = candidateLengths[candidateIndex];

            const int candidateIntsHiLo = SequenceHelpers::getEncodedNumInts2BitHiLo(candidateLength);
            handle.candidateConversionBuffer.resize(candidateIntsHiLo);

            SequenceHelpers::convert2BitTo2BitHiLo(
                handle.candidateConversionBuffer.data(),
                candidate2Bit,
                candidateLength
            );

            *curIter = cpuShiftedHammingDistancePopcount2BitHiLoWithDirection<direction>(
                            handle,
                            handle.anchorConversionBuffer.data(),
                            anchorLength,
                            handle.candidateConversionBuffer.data(),
                            candidateLength,
                            min_overlap,
                            maxErrorRate,
                            min_overlap_ratio
                        );
        }

        return curIter;
    }



    template<class Iter>
    Iter
    cpuShiftedHammingDistancePopcount2Bit(
            CpuAlignmentHandle& handle,
            Iter destinationBegin,
            const unsigned int* anchor2Bit,
            int anchorLength,
            const unsigned int* candidates2Bit,
            int candidatePitchInInts,
            const int* candidateLengths,
            int numCandidates,
            int min_overlap,
            float maxErrorRate,
            float min_overlap_ratio) noexcept{

        const int newanchorInts = SequenceHelpers::getEncodedNumInts2BitHiLo(anchorLength);

        handle.anchorConversionBuffer.resize(newanchorInts);

        SequenceHelpers::convert2BitTo2BitHiLo(
            handle.anchorConversionBuffer.data(),
            anchor2Bit,
            anchorLength
        );

        auto curIter = destinationBegin;

        for(int candidateIndex = 0; candidateIndex < numCandidates; candidateIndex++, ++curIter){
            const unsigned int* candidate2Bit = candidates2Bit + candidatePitchInInts * candidateIndex;
            const int candidateLength = candidateLengths[candidateIndex];

            const int candidateIntsHiLo = SequenceHelpers::getEncodedNumInts2BitHiLo(candidateLength);
            handle.candidateConversionBuffer.resize(candidateIntsHiLo);

            SequenceHelpers::convert2BitTo2BitHiLo(
                handle.candidateConversionBuffer.data(),
                candidate2Bit,
                candidateLength
            );

            *curIter = cpuShiftedHammingDistancePopcount2BitHiLo(
                            handle,
                            handle.anchorConversionBuffer.data(),
                            anchorLength,
                            handle.candidateConversionBuffer.data(),
                            candidateLength,
                            min_overlap,
                            maxErrorRate,
                            min_overlap_ratio
                        );
        }

        return curIter;
    }





    template<class Iter>
    Iter
    cpu_multi_shifted_hamming_distance_popcount_updated(Iter destinationbegin,
                                                const unsigned int* anchorHiLo,
                                                int anchorLength,
                                                const std::vector<unsigned int>& querydata,
                                                const std::vector<int>& queryLengths,
                                                int pitchIntsPerSequence,
                                                int min_overlap,
                                                float maxErrorRate,
                                                float min_overlap_ratio) noexcept{

        assert(anchorLength > 0);
        if(queryLengths.size() == 0) return destinationbegin;

        auto popcount = [](auto i){return __builtin_popcount(i);};

        auto identity = [](auto i){return i;};

        auto hammingDistanceWithShift = [&](int shift, int overlapsize, int max_errors,
                                            unsigned int* shiftptr_hi, unsigned int* shiftptr_lo, auto transfunc1,
                                            int shiftptr_size,
                                            const unsigned int* otherptr_hi, const unsigned int* otherptr_lo,
                                            auto transfunc2){

            const int shiftamount = shift == 0 ? 0 : 1;

            shiftBitArrayLeftBy(shiftptr_hi, shiftptr_size / 2, shiftamount, transfunc1);
            shiftBitArrayLeftBy(shiftptr_lo, shiftptr_size / 2, shiftamount, transfunc1);

            const int score = hammingdistanceHiLo(shiftptr_hi,
                                                    shiftptr_lo,
                                                    otherptr_hi,
                                                    otherptr_lo,
                                                    overlapsize,
                                                    overlapsize,
                                                    max_errors,
                                                    transfunc1,
                                                    transfunc2,
                                                    popcount);

            return score;
        };


        const int nQueries = int(queryLengths.size());

        Iter destination = destinationbegin;

        const unsigned int* const anchor = anchorHiLo;
        const int anchorints = SequenceHelpers::getEncodedNumInts2BitHiLo(anchorLength);

        std::vector<unsigned int> shiftbuffer(pitchIntsPerSequence);

        const unsigned int* const anchorBackup_hi = (const unsigned int*)(anchor);
        const unsigned int* const anchorBackup_lo = ((const unsigned int*)anchor) + anchorints / 2;

        for(int index = 0; index < nQueries; index++){
            const unsigned int* const query = querydata.data() + pitchIntsPerSequence * index;
            const int queryLength = queryLengths[index];

            const int queryints = SequenceHelpers::getEncodedNumInts2BitHiLo(queryLength);
            const unsigned int* const queryBackup_hi = query;
            const unsigned int* const queryBackup_lo = query + queryints / 2;

            const int totalbases = anchorLength + queryLength;
            const int minoverlap = std::max(min_overlap, int(float(anchorLength) * min_overlap_ratio));


            int bestScore = totalbases; // score is number of mismatches
            int bestShift = -queryLength; // shift of query relative to anchor. shift < 0 if query begins before anchor

            auto handle_shift = [&](int shift, int overlapsize,
                                    unsigned int* shiftptr_hi, unsigned int* shiftptr_lo, auto transfunc1,
                                    int shiftptr_size,
                                    const unsigned int* otherptr_hi, const unsigned int* otherptr_lo,
                                    auto transfunc2){

                const int max_errors_excl = std::min(int(float(overlapsize) * maxErrorRate),
                                                bestScore - totalbases + 2*overlapsize);

                if(max_errors_excl > 0){

                    const int mismatches = hammingDistanceWithShift(shift, overlapsize, max_errors_excl,
                                                            shiftptr_hi,shiftptr_lo, transfunc1,
                                                            shiftptr_size,
                                                            otherptr_hi, otherptr_lo, transfunc2);

                    const int score = (mismatches < max_errors_excl ?
                                    mismatches + totalbases - 2*overlapsize // non-overlapping regions count as mismatches
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

            std::copy(anchor, anchor + anchorints, shiftbuffer.begin());
            unsigned int* shiftbuffer_hi = shiftbuffer.data();
            unsigned int* shiftbuffer_lo = shiftbuffer.data() + anchorints / 2;


            for(int shift = 0; shift < anchorLength - minoverlap + 1; ++shift){
                const int overlapsize = std::min(anchorLength - shift, queryLength);
                bool b = handle_shift(shift, overlapsize,
                                    shiftbuffer_hi, shiftbuffer_lo, identity,
                                    anchorints,
                                    queryBackup_hi, queryBackup_lo, identity);
                if(!b){
                    break;
                }
            }

            std::copy(query, query + queryints, shiftbuffer.begin());
            shiftbuffer_hi = shiftbuffer.data();
            shiftbuffer_lo = shiftbuffer.data() + queryints / 2;

            for(int shift = -1; shift >= -queryLength + minoverlap; --shift){
                const int overlapsize = std::min(anchorLength, queryLength + shift);

                bool b = handle_shift(shift, overlapsize,
                                    shiftbuffer_hi, shiftbuffer_lo, identity,
                                    queryints,
                                    anchorBackup_hi, anchorBackup_lo, identity);

                if(!b){
                    break;
                }
            }

            AlignmentResult& alignmentresult = *destination;
            alignmentresult.isValid = (bestShift != -queryLength);

            const int queryoverlapbegin_incl = std::max(-bestShift, 0);
            const int queryoverlapend_excl = std::min(queryLength, anchorLength - bestShift);
            const int overlapsize = queryoverlapend_excl - queryoverlapbegin_incl;
            const int opnr = bestScore - totalbases + 2*overlapsize;

            alignmentresult.score = bestScore;
            alignmentresult.overlap = overlapsize;
            alignmentresult.shift = bestShift;
            alignmentresult.nOps = opnr;

            std::advance(destination, 1);
        }

        return destination;
    }




    AlignmentResult
    cpu_shifted_hamming_distance_popcount(const char* anchor,
                                int anchorLength,
                                const char* query,
                                int queryLength,
                                int min_overlap,
                                float maxErrorRate,
                                float min_overlap_ratio) noexcept;


    std::vector<AlignmentResult>
    cpu_multi_shifted_hamming_distance_popcount(const char* anchor_charptr,
                                int anchorLength,
                                const std::vector<char>& querydata,
                                const std::vector<int>& queryLengths,
                                int max_sequence_bytes,
                                int min_overlap,
                                float maxErrorRate,
                                float min_overlap_ratio) noexcept;


} //namespace shd

using SHDResult = shd::AlignmentResult;


} //namespace care::cpu

} //namespace care




#endif
