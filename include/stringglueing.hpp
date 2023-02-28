#ifndef CARE_STRING_GLUEING_HPP
#define CARE_STRING_GLUEING_HPP

#include <cpu_alignment.hpp>
#include <hostdevicefunctions.cuh>

#include <algorithm>
#include <optional>
#include <string>
#include <string_view>
#include <cassert>

namespace care{


struct GlueDecision{
    /*
        strings s1 and s2 will be combined to produce a string of length resultlength
        in the resulting string, s1[0] will be at index s1FirstResultIndex, 
        and s2[0] will be at index s2FirstResultIndex
    */
    int s1FirstResultIndex{};
    int s2FirstResultIndex{};
    int resultlength{};
    std::string_view s1{};
    std::string_view s2{};
};

struct GlueDecisionWithQuality{
    /*
        strings s1 and s2 will be combined to produce a string of length resultlength
        in the resulting string, s1[0] will be at index s1FirstResultIndex, 
        and s2[0] will be at index s2FirstResultIndex
    */
    int s1FirstResultIndex{};
    int s2FirstResultIndex{};
    int resultlength{};
    std::string_view s1{};
    std::string_view s2{};
    std::string_view q1{};
    std::string_view q2{};
};


struct MismatchRatioGlueDecider{
public:
    MismatchRatioGlueDecider(int overlapLowerBound_, float mismatchRatioUpperBound_)
        : overlapLowerBound(overlapLowerBound_), mismatchRatioUpperBound(mismatchRatioUpperBound_)
    {

    }

    // Given strings s1 and s2, check if it is possible to overlap them and glue them together to produce string of length resultlength
    //Overlap must have a length of at least overlapLowerBound, and must contain at most mismatchRatioUpperBound * length mismatches
    // If this is not possible, the result is empty
    std::optional<GlueDecision> operator()(std::string_view s1, std::string_view s2, int resultlength) const{
        std::optional<GlueDecisionWithQuality> qresult = compute(s1, s2, resultlength, std::nullopt, std::nullopt);

        if(qresult.has_value()){
            GlueDecision result;
            result.s1FirstResultIndex = qresult->s1FirstResultIndex;
            result.s2FirstResultIndex = qresult->s2FirstResultIndex;
            result.resultlength = qresult->resultlength;
            result.s1 = qresult->s1;
            result.s2 = qresult->s2;

            return result;
        }else{
            return std::nullopt;
        }
    }

    std::optional<GlueDecisionWithQuality> operator()(std::string_view s1, std::string_view s2, int resultlength, std::string_view q1, std::string_view q2) const{
        assert(s1.size() == q1.size());
        assert(s2.size() == q2.size());

        std::optional<GlueDecisionWithQuality> qresult = compute(s1, s2, resultlength, q1, q2);

        return qresult;
    }

private:
    std::optional<GlueDecisionWithQuality> compute(std::string_view s1, std::string_view s2, int resultlength, std::optional<std::string_view> q1, std::optional<std::string_view> q2) const{
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
            const int ham = cpu::hammingDistanceOverlap(
                s1begin + (s1EndInCombined+1) - overlapSize, s1end, 
                s2begin + s2Start, s2end
            );
            const float mismatchRatio = float(ham) / float(overlapSize);

            if(fleq(mismatchRatio, mismatchRatioUpperBound)){
                GlueDecisionWithQuality decision;

                decision.s1FirstResultIndex = 0;
                decision.s2FirstResultIndex = s2BeginInCombined >= 0 ? s2BeginInCombined : 0;
                decision.resultlength = resultlength;
                decision.s1 = s1;
                decision.s2 = s2;
                decision.s1.remove_suffix(s1length - (s1EndInCombined + 1));
                assert(s2Start <= int(decision.s2.size()));
                decision.s2.remove_prefix(s2Start);

                if(q1.has_value()){
                    assert(q1->size() == s1.size());
                    decision.q1 = *q1;
                    decision.q1.remove_suffix(s1length - (s1EndInCombined + 1));
                }

                if(q2.has_value()){
                    assert(q2->size() == s2.size());
                    decision.q2 = *q2;
                    decision.q2.remove_prefix(s2Start);
                }

                return decision;
            }else{
                return {}; //no result
            }
        }else{
            return {}; //no result
        }
    }

    int overlapLowerBound;
    float mismatchRatioUpperBound;
};



struct MatchLengthGlueDecider{
public:
    MatchLengthGlueDecider(int overlapLowerBound_, float matchLengthLowerBound_)
        : overlapLowerBound(overlapLowerBound_), matchLengthLowerBound(matchLengthLowerBound_)
    {

    }

    // Given strings s1 and s2, check if it is possible to overlap them and glue them together to produce string of length resultlength
    // overlap must have a length of at least overlapLowerBound, and must contain a match of at least matchLengthLowerBound consecutive bases
    // If this is not possible, the result is empty
    std::optional<GlueDecision> operator()(std::string_view s1, std::string_view s2, int resultlength) const{
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
            const int l = cpu::longestMatch(
                s1begin + (s1EndInCombined+1) - overlapSize, s1end, 
                s2begin + s2Start, s2end
            );

            if(l >= matchLengthLowerBound){
                GlueDecision decision;

                decision.s1FirstResultIndex = 0;
                decision.s2FirstResultIndex = s2BeginInCombined >= 0 ? s2BeginInCombined : 0;
                decision.resultlength = resultlength;
                decision.s1 = s1;
                decision.s2 = s2;
                decision.s1.remove_suffix(s1length - (s1EndInCombined + 1));
                assert(s2Start <= int(decision.s2.size()));
                decision.s2.remove_prefix(s2Start);

                return decision;
            }else{
                return {}; //no result
            }
        }else{
            return {}; //no result
        }
    }

private:
    int overlapLowerBound{};
    int matchLengthLowerBound{};
};


/*
    Very naive gluer. Take as much from s2 as possible, fill remaining positions at the left end with s1
*/
struct NaiveGluer{
public:

    std::string operator()(const GlueDecision& g) const{
        assert(g.s1FirstResultIndex == 0);
        assert(int(g.s1.size()) <= g.resultlength);
        assert(int(g.s2.size()) <= g.resultlength);

        std::string result(g.resultlength, '#');

        const int s11remaining = g.resultlength - g.s2.size();
        const auto it = std::copy_n(g.s1.begin(), s11remaining, result.begin());

        auto end = std::copy(g.s2.begin(), g.s2.end(), it);
        if(end != result.end())
            throw std::runtime_error("Error glueing");

        return result;
    }

};


/*
    result string will begin with prefix s1[0 : min(s1length, originalReadLength)]
    result string will end with suffix s2[max(0, s2length - originalReadLength) : min(s2length, originalReadLength)]
    Gap is filled with either s1 or s2, depending on the weight of the position
*/
struct WeightedGapGluer{
public:
    WeightedGapGluer(int origLength) : originalReadLength(origLength){}

    std::string operator()(const GlueDecision& g) const{
        assert(g.s1FirstResultIndex == 0);
        assert(int(g.s1.size()) <= g.resultlength);
        assert(int(g.s2.size()) <= g.resultlength);

        std::string result(g.resultlength, '#');

        const int numToCopys1 = std::min(
            g.resultlength,
            std::min(int(g.s1.size()), originalReadLength)
        );
        const int numToCopys2 = std::min(
            g.resultlength,
            std::min(int(g.s2.size()), originalReadLength)
        );

        const auto gapbegin = std::copy_n(g.s1.begin(), numToCopys1, result.begin());      
        const auto gapend = result.begin() + result.size() - numToCopys2;
        std::copy_n(g.s2.begin() + g.s2.size() - numToCopys2, numToCopys2, gapend);

        //fill the gap
        if(std::distance(result.begin(), gapbegin) < std::distance(result.begin(), gapend)){

            auto getweight = [&](const auto& s, int pos){
                if(pos < 0 || pos >= int(s.size())){
                    return 0.0f;
                }else{
                    // original positions have weight 1, weight of other positions decreases with distance to original
                    if(pos < originalReadLength){
                        return 1.0f;
                    }else{
                        //linear interpolation
                        const int unoriginalpositions = s.size() - originalReadLength;
                        const float f = 1.0f / (unoriginalpositions + 1);

                        return 1.0f - f * (pos + 1 - originalReadLength);
                    }
                }
            };

            const int gapsize = std::distance(gapbegin, gapend);
            const int firstGapPos = std::distance(result.begin(), gapbegin);

            auto iter = gapbegin;
            for(int i = 0; i < gapsize; i++){
                const int positionInS1 = firstGapPos + i;
                const int positionInS2 = firstGapPos + i - g.s2FirstResultIndex;

                const float w1 = getweight(g.s1, positionInS1);
                const float w2 = getweight(g.s2, positionInS2);
                assert(!feq(w1, 0.0f) || !feq(w2, 0.0f));

                if(fgeq(w1, w2)){
                    assert(positionInS1 < int(g.s1.size()));
                    *iter = g.s1[positionInS1];
                }else{
                    assert(positionInS2 < int(g.s2.size()));
                    *iter = g.s2[positionInS2];
                }
                //*iter = fgeq(getweight(g.s1, positionInS1), getweight(g.s2, positionInS2)) ? g.s1[positionInS1] : g.s2[positionInS2];

                ++iter;
            }
        }


        return result;
    }
private:
    int originalReadLength{};
};










/*
    result string will begin with prefix s1[0 : min(s1length, originalReadLength)]
    result string will end with suffix s2[max(0, s2length - originalReadLength) : min(s2length, originalReadLength)]
    Gap is filled with either s1 or s2, depending on the weight of the position
*/
struct QualityWeightedGapGluer{
public:
    QualityWeightedGapGluer(int origLength1, int origLength2) : originalReadLength1(origLength1), originalReadLength2(origLength2){}

    std::pair<std::string, std::string> operator()(const GlueDecisionWithQuality& g) const{
        assert(g.s1FirstResultIndex == 0);
        assert(int(g.s1.size()) <= g.resultlength);
        assert(int(g.s2.size()) <= g.resultlength);

        std::string result(g.resultlength, '#');
        std::string resultquality(g.resultlength, '!');

        const int numToCopys1 = std::min(
            g.resultlength,
            std::min(int(g.s1.size()), originalReadLength1)
        );
        const int numToCopys2 = std::min(
            g.resultlength,
            std::min(int(g.s2.size()), originalReadLength2)
        );

        const auto gapbegin = std::copy_n(g.s1.begin(), numToCopys1, result.begin());      
        const auto gapend = result.begin() + result.size() - numToCopys2;
        std::copy_n(g.s2.begin() + g.s2.size() - numToCopys2, numToCopys2, gapend);

        assert(int(g.q1.size()) >= numToCopys1);
        const auto gapbeginq = std::copy_n(g.q1.begin(), numToCopys1, resultquality.begin());      
        const auto gapendq = resultquality.begin() + resultquality.size() - numToCopys2;
        std::copy_n(g.q2.begin() + g.q2.size() - numToCopys2, numToCopys2, gapendq);


        //fill the gap
        if(std::distance(result.begin(), gapbegin) < std::distance(result.begin(), gapend)){

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

            const int gapsize = std::distance(gapbegin, gapend);
            const int firstGapPos = std::distance(result.begin(), gapbegin);

            auto iter = gapbegin;
            auto iterq = gapbeginq;

            for(int i = 0; i < gapsize; i++){
                const int positionInS1 = firstGapPos + i;
                const int positionInS2 = firstGapPos + i - g.s2FirstResultIndex;

                const float w1 = getweight(g.s1, originalReadLength1, g.q1, positionInS1);
                const float w2 = getweight(g.s2, originalReadLength2, g.q2, positionInS2);
                assert(!feq(w1, 0.0f) || !feq(w2, 0.0f));

                if(fgeq(w1, w2)){
                    assert(positionInS1 < int(g.s1.size()));
                    *iter = g.s1[positionInS1];
                    *iterq = g.q1[positionInS1];
                }else{
                    assert(positionInS2 < int(g.s2.size()));
                    *iter = g.s2[positionInS2];
                    *iterq = g.q2[positionInS2];
                }

                ++iter;
                ++iterq;
            }
        }


        return std::make_pair(std::move(result), std::move(resultquality));
    }
private:
    int originalReadLength1{};
    int originalReadLength2{};
};








} //namespace care

#endif