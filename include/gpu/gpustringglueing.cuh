#ifndef CARE_GPU_STRING_GLUEING_HPP
#define CARE_GPU_STRING_GLUEING_HPP

#include <hpc_helpers.cuh>
#include <mystringview.hpp>
#include <hostdevicefunctions.cuh>

#include <cub/cub.cuh>

#include <cassert>

namespace care{
namespace gpu{

template<class StringViewS1, class StringViewS2, class StringViewQ1, class StringViewQ2>
struct GlueDecisionWithQuality{
    /*
        strings s1 and s2 will be combined to produce a string of length resultlength
        in the resulting string, s1[0] will be at index s1FirstResultIndex, 
        and s2[0] will be at index s2FirstResultIndex
    */
    bool valid{};
    int s1FirstResultIndex{};
    int s2FirstResultIndex{};
    int resultlength{};
    StringViewS1 s1{};
    StringViewS2 s2{};
    StringViewQ1 q1{};
    StringViewQ2 q2{};
};

template<int blocksize>
struct MismatchRatioGlueDecider{
public:
    using BlockReduce = cub::BlockReduce<int, blocksize>;

    struct TempStorage{
        typename BlockReduce::TempStorage reduce;
        int broadcast_int;
    };

    __device__
    MismatchRatioGlueDecider(TempStorage& tmp, int overlapLowerBound_, float mismatchRatioUpperBound_)
        : overlapLowerBound(overlapLowerBound_), 
        mismatchRatioUpperBound(mismatchRatioUpperBound_),
        tempStorage(tmp)
    {

    }

    // Given strings s1 and s2, check if it is possible to overlap them and glue them together to produce string of length resultlength
    //Overlap must have a length of at least overlapLowerBound, and must contain at most mismatchRatioUpperBound * length mismatches
    // If this is not possible, the result is empty
    template<class StringViewS1, class StringViewS2, class StringViewQ1, class StringViewQ2>
    __device__
    GlueDecisionWithQuality<StringViewS1, StringViewS2, StringViewQ1, StringViewQ2> 
    operator()(StringViewS1 s1, StringViewS2 s2, int resultlength, StringViewQ1 q1, StringViewQ2 q2) const{
        assert(s1.size() == q1.size());
        assert(s2.size() == q2.size());

        return compute(s1, s2, resultlength, q1, q2);
    }

private:

    template<class StringViewS1, class StringViewS2, class StringViewQ1, class StringViewQ2>
    //q1 and q2 should be optional<StringViewQ1 / Q2>
    __device__
    GlueDecisionWithQuality<StringViewS1, StringViewS2, StringViewQ1, StringViewQ2> 
    compute(StringViewS1 s1, StringViewS2 s2, int resultlength, StringViewQ1 q1, StringViewQ2 q2) const{
        const auto s1begin = s1.begin();
        const auto s1end = s1.end();
        const auto s2begin = s2.begin();
        const auto s2end = s2.end();

        const int s1length = std::distance(s1begin, s1end);
        const int s2length = std::distance(s2begin, s2end);

        //the last position of rlString should be at position x in the combined string
        const int x = resultlength - 1;

        const int s2BeginInCombined = x - s2length + 1; //theoretical position of first character of s2 in the combined string
        const int s1EndInCombined = std::min(x, s1length - 1); //inclusive
        const int overlapSize = std::min(
            std::max(0, s1EndInCombined - s2BeginInCombined + 1), // overlap can be at most the positions between begin of s2 and end of s1
            std::min(
                std::min(s1length, s2length), // overlap cannot be longer than any of both strings
                resultlength //overlap cannot be longer than specified segment length
            )
        );

        if(overlapSize >= overlapLowerBound){
            const int s2Start = std::max(0, -s2BeginInCombined);
            const int numPositionsToCheck = min(
                int(s1.size()) - ((s1EndInCombined+1) - overlapSize),
                int(s2.size()) - s2Start
            );
            int ham = 0;
            for(int i = threadIdx.x; i < numPositionsToCheck; i += blocksize){
                const char c1 = s1[(s1EndInCombined+1) - overlapSize + i];
                const char c2 = s2[s2Start + i];
                ham += (c1 != c2) ? 1 : 0;
            }

            ham = BlockReduce(tempStorage.reduce).Sum(ham);
            __syncthreads();
            if(threadIdx.x == 0){
                tempStorage.broadcast_int = ham;
            }
            __syncthreads();
            ham = tempStorage.broadcast_int;
            __syncthreads();
            const float mismatchRatio = float(ham) / float(overlapSize);

            if(fleq(mismatchRatio, mismatchRatioUpperBound)){
                GlueDecisionWithQuality<StringViewS1, StringViewS2, StringViewQ1, StringViewQ2> decision;
                decision.valid = true;

                decision.s1FirstResultIndex = 0;
                decision.s2FirstResultIndex = s2BeginInCombined >= 0 ? s2BeginInCombined : 0;
                decision.resultlength = resultlength;
                decision.s1 = s1;
                decision.s2 = s2;
                decision.s1.remove_suffix(s1length - (s1EndInCombined + 1));
                assert(s2Start <= int(decision.s2.size()));
                decision.s2.remove_prefix(s2Start);

                //if(q1.has_value()){
                    assert(q1.size() == s1.size());
                    //decision.q1 = *q1;
                    decision.q1 = q1;
                    decision.q1.remove_suffix(s1length - (s1EndInCombined + 1));
                //}

                //if(q2.has_value()){
                    assert(q2.size() == s2.size());
                    //decision.q2 = *q2;
                    decision.q2 = q2;
                    decision.q2.remove_prefix(s2Start);
                //}

                return decision;
            }else{
                GlueDecisionWithQuality<StringViewS1, StringViewS2, StringViewQ1, StringViewQ2> decision;
                decision.valid = false; //no result
                return decision;
            }
        }else{
            GlueDecisionWithQuality<StringViewS1, StringViewS2, StringViewQ1, StringViewQ2> decision;
            decision.valid = false; //no result
            return decision;
        }
    }

    int overlapLowerBound;
    float mismatchRatioUpperBound;
    TempStorage& tempStorage;
};



struct QualityWeightedGapGluer{
public:
    __device__
    QualityWeightedGapGluer(int origLength1, int origLength2) : originalReadLength1(origLength1), originalReadLength2(origLength2){}

    template<class Decision, class Group>
    __device__
    void operator()(Group group, const Decision& g, char* resultSequence, char* resultQualities) const{
        assert(g.s1FirstResultIndex == 0);
        assert(int(g.s1.size()) <= g.resultlength);
        assert(int(g.s2.size()) <= g.resultlength);

        const int numToCopys1 = min(
            g.resultlength,
            min(int(g.s1.size()), originalReadLength1)
        );
        const int numToCopys2 = min(
            g.resultlength,
            min(int(g.s2.size()), originalReadLength2)
        );

        const int gapbegin = numToCopys1; 
        const int gapend = g.resultlength - numToCopys2;

        for(int i = group.thread_rank(); i < numToCopys1; i += group.size()){
            resultSequence[i] = g.s1[i];
        }

        for(int i = group.thread_rank(); i < numToCopys2; i += group.size()){
            resultSequence[gapend + i] = g.s2[g.s2.size() - numToCopys2 + i];
        }

        assert(int(g.q1.size()) >= numToCopys1);

        for(int i = group.thread_rank(); i < numToCopys1; i += group.size()){
            resultQualities[i] = g.q1[i];
        }
        for(int i = group.thread_rank(); i < numToCopys2; i += group.size()){
            resultQualities[gapend + i] = g.q2[g.q2.size() - numToCopys2 + i];
        }


        //fill the gap. for each position, use base of sequence with highest quality
        if(gapbegin < gapend){

            auto getweight = [&](const auto& sequence, int origlength, const auto& quality, int pos){
                assert(sequence.size() == quality.size());
                if(pos < 0 || pos >= int(sequence.size())){
                    return 0.0f;
                }else{
                    // original positions have weight 1
                    if(pos < origlength){
                        return 1.0f;
                    }else{
                        return getQualityWeight(quality[pos]);
                    }
                }
            };

            const int gapsize = gapend - gapbegin;
            const int firstGapPos = gapbegin;

            for(int i = group.thread_rank(); i < gapsize; i += group.size()){
                const int positionInS1 = firstGapPos + i;
                const int positionInS2 = firstGapPos + i - g.s2FirstResultIndex;

                const float w1 = getweight(g.s1, originalReadLength1, g.q1, positionInS1);
                const float w2 = getweight(g.s2, originalReadLength2, g.q2, positionInS2);
                assert((w1 != 0.0f) || (w2 != 0.0f));

                if(fgeq(w1, w2)){
                    assert(positionInS1 < int(g.s1.size()));
                    resultSequence[gapbegin + i] = g.s1[positionInS1];
                    resultQualities[gapbegin + i] = g.q1[positionInS1];
                }else{
                    assert(positionInS2 < int(g.s2.size()));
                    resultSequence[gapbegin + i] = g.s2[positionInS2];
                    resultQualities[gapbegin + i] = g.q2[positionInS2];
                }
            }
        }
    }
private:
    int originalReadLength1{};
    int originalReadLength2{};
};















} //namespace gpu
} //namespace care

#endif