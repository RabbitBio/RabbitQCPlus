 
#ifndef CARE_GPU_MSA_CUH
#define CARE_GPU_MSA_CUH

#ifdef __NVCC__

#include <config.hpp>
#include <hostdevicefunctions.cuh>

#include <alignmentorientation.hpp>

#include <sequencehelpers.hpp>
#include <hpc_helpers.cuh>

#include <cassert>
#include <cstdint>
#include <type_traits>

#include <cub/cub.cuh>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include <cooperative_groups.h>

namespace cg = cooperative_groups;

namespace care{
namespace gpu{

    struct MSAColumnProperties{
        int anchorColumnsBegin_incl;
        int anchorColumnsEnd_excl;
        int firstColumn_incl;
        int lastColumn_excl;
    };

    struct GpuMSAProperties{
        float avg_support;
        float min_support;
        int min_coverage;
        int max_coverage;
    };

    struct GpuSingleMSA{
    public:

        __device__ __forceinline__
        void printCounts(int printbegin = -1, int printend = -1){
            if(printbegin == -1){
                printbegin = 0;
            }
            if(printend == -1){
                printend = columnPitchInElements;
            }
            printf("%d - %d, first %d, last %d, anchor %d\n", 
                printbegin, printend, columnProperties->firstColumn_incl, columnProperties->lastColumn_excl, columnProperties->anchorColumnsBegin_incl);
            printf("Counts A:\n");
            for(int i = printbegin; i < printend; i++){
                printf("%d ", counts[0 * columnPitchInElements + i]);
            }
            printf("\n");
            printf("Counts C:\n");
            for(int i = printbegin; i < printend; i++){
                printf("%d ", counts[1 * columnPitchInElements + i]);
            }
            printf("\n");
            printf("Counts G:\n");
            for(int i = printbegin; i < printend; i++){
                printf("%d ", counts[2 * columnPitchInElements + i]);
            }
            printf("\n");
            printf("Counts T:\n");
            for(int i = printbegin; i < printend; i++){
                printf("%d ", counts[3 * columnPitchInElements + i]);
            }
            printf("\n");
        }

        __device__ __forceinline__
        void printWeights(int printbegin = -1, int printend = -1){
            if(printbegin == -1){
                printbegin = 0;
            }
            if(printend == -1){
                printend = columnPitchInElements;
            }
            printf("%d - %d, first %d, last %d, anchor %d\n", 
                printbegin, printend, columnProperties->firstColumn_incl, columnProperties->lastColumn_excl, columnProperties->anchorColumnsBegin_incl);
            printf("Weights A:\n");
            for(int i = printbegin; i < printend; i++){
                printf("%.7f ", weights[0 * columnPitchInElements + i]);
            }
            printf("\n");
            printf("Weights C:\n");
            for(int i = printbegin; i < printend; i++){
                printf("%.7f ", weights[1 * columnPitchInElements + i]);
            }
            printf("\n");
            printf("Weights G:\n");
            for(int i = printbegin; i < printend; i++){
                printf("%.7f ", weights[2 * columnPitchInElements + i]);
            }
            printf("\n");
            printf("Weights T:\n");
            for(int i = printbegin; i < printend; i++){
                printf("%.7f ", weights[3 * columnPitchInElements + i]);
            }
            printf("\n");
        }

        template<class ThreadGroup>
        __device__ __forceinline__
        void checkAfterBuild(ThreadGroup& group, int anchorIndex = -1, int line = -1, read_number anchorReadId = 0){
    
            const int firstColumn_incl = columnProperties->firstColumn_incl;
            const int lastColumn_excl = columnProperties->lastColumn_excl;
    
            for(int column = firstColumn_incl + group.thread_rank(); 
                    column < lastColumn_excl; 
                    column += group.size()){

                const int* const mycounts = counts + column;
                const float* const myweights = weights + column;
                float sumOfWeights = 0.0f;
    
                #pragma unroll
                for(int k = 0; k < 4; k++){
                    const int count = mycounts[k * columnPitchInElements];
                    const float weight = myweights[k * columnPitchInElements];

                    if(count > 0 && fleq(weight, 0.0f)){
                        // printf("msa check failed! firstColumn_incl %d, lastColumn_excl %d, anchorIndex %d, anchorReadId %u, column %d, base %d, count %d, weight %.20f, line %d\n",
                        //     firstColumn_incl, lastColumn_excl, anchorIndex, anchorReadId, column, k, count, weight, line);
                        assert(false);
                    }
    
                    if(count < 0){
                        // printf("msa check failed! firstColumn_incl %d, lastColumn_excl %d, anchorIndex %d, anchorReadId %u, column %d, base %d, count %d, weight %.20f, line %d\n",
                        //     firstColumn_incl, lastColumn_excl, anchorIndex, anchorReadId, column, k, count, weight, line);
                        assert(false);
                    }

                    if(count == 0 && (weight - 1e-5) > 0.0f){
                        // printf("msa check failed! firstColumn_incl %d, lastColumn_excl %d, anchorIndex %d, anchorReadId %u, column %d, base %d, count %d, weight %.20f, line %d\n",
                        //     firstColumn_incl, lastColumn_excl, anchorIndex, anchorReadId, column, k, count, weight, line);
                        assert(false);
                    }
    
                    sumOfWeights += weight;
                }

                if(sumOfWeights == 0){
                    //printf("anchorIndex %d, column %d, codeline %d\n", anchorIndex, column, line);
                    assert(false);
                }

                const int cov = coverages[column];
                if(cov <= 0){
                    // printf("msa check failed! firstColumn_incl %d, lastColumn_excl %d, anchorIndex %d, anchorReadId %u, column %d, cov %d, line %d\n", 
                    //     firstColumn_incl, lastColumn_excl, anchorIndex, anchorReadId, column, cov, line);
                    assert(false);
                }
    
            }

            group.sync();
        }


        //only group.thread_rank() == 0 returns the correct MSAProperties
        template<
            class ThreadGroup,
            class GroupReduceFloatSum,
            class GroupReduceFloatMin,
            class GroupReduceIntMin,
            class GroupReduceIntMax
        >
        __device__ __forceinline__
        GpuMSAProperties getMSAProperties(
            ThreadGroup& group,
            GroupReduceFloatSum& groupReduceFloatSum,
            GroupReduceFloatMin& groupReduceFloatMin,
            GroupReduceIntMin& groupReduceIntMin,
            GroupReduceIntMax& groupReduceIntMax,
            int firstCol,
            int lastCol //exclusive
        ) const {
            #ifndef NDEBUG
            const int firstColumn_incl = columnProperties->firstColumn_incl;
            const int lastColumn_excl = columnProperties->lastColumn_excl;
            assert(firstColumn_incl <= firstCol && firstCol <= lastColumn_excl);
            assert(firstColumn_incl <= lastColumn_excl && lastColumn_excl <= lastColumn_excl);
            #endif

            float avg_support = 0;
            float min_support = 1.0f;
            int min_coverage = std::numeric_limits<int>::max();
            int max_coverage = std::numeric_limits<int>::min();

            for(int i = firstCol + group.thread_rank(); 
                    i < lastCol; 
                    i += group.size()){
                
                assert(firstColumn_incl <= i && i < lastColumn_excl);

                avg_support += support[i];
                min_support = min(support[i], min_support);
                min_coverage = min(coverages[i], min_coverage);
                max_coverage = max(coverages[i], max_coverage);
            }

            avg_support = groupReduceFloatSum(avg_support);
            avg_support /= (lastCol - firstCol);

            min_support = groupReduceFloatMin(min_support);

            min_coverage = groupReduceIntMin(min_coverage);

            max_coverage = groupReduceIntMax(max_coverage);


            GpuMSAProperties msaProperties;

            msaProperties.min_support = min_support;
            msaProperties.avg_support = avg_support;
            msaProperties.min_coverage = min_coverage;
            msaProperties.max_coverage = max_coverage;

            return msaProperties;
        }

        //computes column properties. only return value of thread_rank 0 is valid
        template<
            class ThreadGroup, 
            class GroupReduceIntMin, 
            class GroupReduceIntMax
        >
        __device__ __forceinline__
        static MSAColumnProperties computeColumnProperties(
                ThreadGroup& group,
                GroupReduceIntMin& groupReduceIntMin,
                GroupReduceIntMax& groupReduceIntMax,
                const int* __restrict__ goodCandidateIndices,
                int numGoodCandidates,
                const int* __restrict__ shifts,
                const int anchorLength,
                const int* __restrict__ candidateLengths
        ){

            int startindex = 0;
            int endindex = anchorLength;

            for(int k = group.thread_rank(); k < numGoodCandidates; k += group.size()) {
                const int localCandidateIndex = goodCandidateIndices[k];

                const int shift = shifts[localCandidateIndex];
                const int queryLength = candidateLengths[localCandidateIndex];

                const int queryEndsAt = queryLength + shift;
                startindex = min(startindex, shift);
                endindex = max(endindex, queryEndsAt);
            }

            startindex = groupReduceIntMin(startindex);
            endindex = groupReduceIntMax(endindex);

            group.sync();

            MSAColumnProperties my_columnproperties;

            my_columnproperties.anchorColumnsBegin_incl = max(-startindex, 0);
            my_columnproperties.anchorColumnsEnd_excl = my_columnproperties.anchorColumnsBegin_incl + anchorLength;
            my_columnproperties.firstColumn_incl = 0;
            my_columnproperties.lastColumn_excl = endindex - startindex;

            return my_columnproperties;
        }


        template<
            class ThreadGroup, 
            class GroupReduceIntMin, 
            class GroupReduceIntMax
        >
        __device__ __forceinline__
        void initColumnProperties(
                ThreadGroup& group,
                GroupReduceIntMin& groupReduceIntMin,
                GroupReduceIntMax& groupReduceIntMax,
                const int* __restrict__ goodCandidateIndices,
                int numGoodCandidates,
                const int* __restrict__ shifts,
                const AlignmentOrientation* __restrict__ alignmentFlags,
                const int anchorLength,
                const int* __restrict__ candidateLengths
        ){

            int startindex = 0;
            int endindex = anchorLength;

            for(int k = group.thread_rank(); k < numGoodCandidates; k += group.size()) {
                const int localCandidateIndex = goodCandidateIndices[k];

                #ifndef NDEBUG
                const AlignmentOrientation flag = alignmentFlags[localCandidateIndex];
                assert(flag != AlignmentOrientation::None);
                #endif

                const int shift = shifts[localCandidateIndex];
                const int queryLength = candidateLengths[localCandidateIndex];


                const int queryEndsAt = queryLength + shift;
                startindex = min(startindex, shift);
                endindex = max(endindex, queryEndsAt);
            }

            startindex = groupReduceIntMin(startindex);
            endindex = groupReduceIntMax(endindex);

            group.sync();


            if(group.thread_rank() == 0) {
                MSAColumnProperties my_columnproperties;

                my_columnproperties.anchorColumnsBegin_incl = max(-startindex, 0);
                my_columnproperties.anchorColumnsEnd_excl = my_columnproperties.anchorColumnsBegin_incl + anchorLength;
                my_columnproperties.firstColumn_incl = 0;
                my_columnproperties.lastColumn_excl = endindex - startindex;

                *columnProperties = my_columnproperties;
            }
        }


        template<
            class ThreadGroup, 
            class GroupReduceIntMin, 
            class GroupReduceIntMax
        >
        __device__ __forceinline__
        void initColumnProperties(
                ThreadGroup& group,
                GroupReduceIntMin& groupReduceIntMin,
                GroupReduceIntMax& groupReduceIntMax,
                int numGoodCandidates,
                const int* __restrict__ shifts,
                const AlignmentOrientation* __restrict__ alignmentFlags,
                const int anchorLength,
                const int* __restrict__ candidateLengths
        ){

            int startindex = 0;
            int endindex = anchorLength;

            for(int k = group.thread_rank(); k < numGoodCandidates; k += group.size()) {
                const int localCandidateIndex = k;

                const int shift = shifts[localCandidateIndex];
                const AlignmentOrientation flag = alignmentFlags[localCandidateIndex];
                const int queryLength = candidateLengths[localCandidateIndex];

                assert(flag != AlignmentOrientation::None);

                const int queryEndsAt = queryLength + shift;
                startindex = min(startindex, shift);
                endindex = max(endindex, queryEndsAt);
            }

            startindex = groupReduceIntMin(startindex);
            endindex = groupReduceIntMax(endindex);

            group.sync();


            if(group.thread_rank() == 0) {
                MSAColumnProperties my_columnproperties;

                my_columnproperties.anchorColumnsBegin_incl = max(-startindex, 0);
                my_columnproperties.anchorColumnsEnd_excl = my_columnproperties.anchorColumnsBegin_incl + anchorLength;
                my_columnproperties.firstColumn_incl = 0;
                my_columnproperties.lastColumn_excl = endindex - startindex;

                *columnProperties = my_columnproperties;
            }
        }

        template<class ThreadGroup>
        __device__ __forceinline__
        bool updateColumnProperties(ThreadGroup& group){

            const int firstColumn_incl = columnProperties->firstColumn_incl;
            const int lastColumn_excl = columnProperties->lastColumn_excl;
            const int numColumnsToCheck = lastColumn_excl - firstColumn_incl;

            int newFirstColumn_incl = -1;
            int newLastColumn_excl = -1;

            bool error = false;

            for(int i = group.thread_rank(); i < numColumnsToCheck-1; i += group.size()){
                const int column = firstColumn_incl + i;

                const int thisCoverage = coverages[column];
                const int nextCoverage = coverages[column+1];

                if(thisCoverage < 0 || nextCoverage < 0){
                    error = true;
                }

                if(thisCoverage == 0 && nextCoverage > 0){
                    newFirstColumn_incl = column+1;
                }

                if(thisCoverage > 0 && nextCoverage == 0){
                    newLastColumn_excl = column+1;
                }
            }

            //there can be at most one thread for which this is true
            if(newFirstColumn_incl != -1){
                columnProperties->firstColumn_incl = newFirstColumn_incl;
            }
            //there can be at most one thread for which this is true
            if(newLastColumn_excl != -1){
                columnProperties->lastColumn_excl = newLastColumn_excl;
            }
            
            group.sync();

            return error;
        }

        template<bool doAdd, class ThreadGroup>
        __device__ 
        void msaAddOrDeleteASequence2BitAtomic(
            ThreadGroup& group,
            const unsigned int* __restrict__ sequence, //not transposed
            int sequenceLength, 
            bool isForward,
            int columnStart,
            float overlapweight,
            const char* __restrict__ quality, //not transposed
            bool canUseQualityScores,
            bool isPairedCandidate
        ){
            constexpr char defaultQualityScore = 'I';

            auto defaultweightfunc = [&](char q){
                return canUseQualityScores ? getQualityWeight(q) * overlapweight : overlapweight;
            };

            // auto weightfuncIncreasedOverlapweightForPairedCandidate = [&](char q){
            //     constexpr float increasefactor = 2.0f;

            //     const float newoverlapweight = isPairedCandidate ? min(1.0f, overlapweight * increasefactor) : overlapweight;
            //     return canUseQualityScores ? getQualityWeight(q) * newoverlapweight : newoverlapweight;
            // };

            // auto weightfuncOverlap1ForPairedCandidate = [&](char q){

            //     const float newoverlapweight = isPairedCandidate ? 1.0f : overlapweight;
            //     return canUseQualityScores ? getQualityWeight(q) * newoverlapweight : newoverlapweight;
            // };

            // auto weightfuncIncreasedQualityweightForPairedCandidate = [&](char q){
            //     constexpr float increasefactor = 2.0f;

            //     if(canUseQualityScores){
            //         if(isPairedCandidate){
            //             return min(1.0f, getQualityWeight(q) * increasefactor) * overlapweight;
            //         }else{
            //             return getQualityWeight(q) * overlapweight;
            //         }
            //     }else{
            //         return overlapweight;
            //     }
            // };

            // auto weightfuncQuality1ForPairedCandidate = [&](char q){
   
            //     if(canUseQualityScores){
            //         if(isPairedCandidate){
            //             return 1.0f * overlapweight;
            //         }else{
            //             return getQualityWeight(q) * overlapweight;
            //         }
            //     }else{
            //         return overlapweight;
            //     }
            // };

            // auto perfectWeightForPairedCandidate = [&](char q){
            //     if(isPairedCandidate){
            //         return 1.0f;
            //     }else{
            //         return defaultweightfunc(q);
            //     }
            // };

            auto weightfunc = defaultweightfunc;

            constexpr int nucleotidesPerInt2Bit = SequenceHelpers::basesPerInt2Bit();
            const int fullInts = sequenceLength / nucleotidesPerInt2Bit;

            for(int intIndex = group.thread_rank(); intIndex < fullInts; intIndex += group.size()){
                const unsigned int currentDataInt = sequence[intIndex];

                #pragma unroll
                for(int k = 0; k < 4; k++){
                    alignas(4) char currentFourQualities[4];

                    if(canUseQualityScores){
                        *((int*)&currentFourQualities[0]) = ((const int*)quality)[intIndex * 4 + k];
                    }else{
                        #pragma unroll
                        for(int l = 0; l < 4; l++){
                            currentFourQualities[l] = defaultQualityScore;
                        }
                    }

                    #pragma unroll
                    for(int l = 0; l < 4; l++){
                        const int posInInt = k * 4 + l;

                        std::int8_t encodedBaseAsInt = SequenceHelpers::getEncodedNucFromInt2Bit(currentDataInt, posInInt);
                        if(!isForward){
                            encodedBaseAsInt = SequenceHelpers::complementBase2Bit(encodedBaseAsInt);
                        }

                        const float weight = weightfunc(currentFourQualities[l]);

                        assert(weight != 0);
                        const int rowOffset = encodedBaseAsInt * columnPitchInElements;
                        const int columnIndex = columnStart 
                                + (isForward ? (intIndex * nucleotidesPerInt2Bit + posInInt) : sequenceLength - 1 - (intIndex * nucleotidesPerInt2Bit + posInInt));
                        
                        atomicAdd(counts + rowOffset + columnIndex, doAdd ? 1 : -1);
                        float n = atomicAdd(weights + rowOffset + columnIndex, doAdd ? weight : -weight);
                        atomicAdd(coverages + columnIndex, doAdd ? 1 : -1);
                    }
                }
            }

            //add remaining positions
            if(sequenceLength % nucleotidesPerInt2Bit != 0){
                const unsigned int currentDataInt = sequence[fullInts];
                const int maxPos = sequenceLength - fullInts * nucleotidesPerInt2Bit;

                for(int posInInt = group.thread_rank(); posInInt < maxPos; posInInt += group.size()){
                    std::int8_t encodedBaseAsInt = SequenceHelpers::getEncodedNucFromInt2Bit(currentDataInt, posInInt);
                    if(!isForward){
                        //reverse complement
                        encodedBaseAsInt = SequenceHelpers::complementBase2Bit(encodedBaseAsInt);
                    }
                    char q = defaultQualityScore;
                    if(canUseQualityScores){
                        q = quality[fullInts * nucleotidesPerInt2Bit + posInInt];
                    }
                    const float weight = weightfunc(q);

                    assert(weight != 0);
                    const int rowOffset = encodedBaseAsInt * columnPitchInElements;
                    const int columnIndex = columnStart 
                        + (isForward ? (fullInts * nucleotidesPerInt2Bit + posInInt) : sequenceLength - 1 - (fullInts * nucleotidesPerInt2Bit + posInInt));

                    atomicAdd(counts + rowOffset + columnIndex, doAdd ? 1 : -1);
                    atomicAdd(weights + rowOffset + columnIndex, doAdd ? weight : -weight);
                    atomicAdd(coverages + columnIndex, doAdd ? 1 : -1);
                } 
            }
            group.sync();
        }


        template<bool doAdd, class ThreadGroup>
        __device__ 
        void msaAddOrDeleteASequence2BitNormal(
            ThreadGroup& group,
            const unsigned int* __restrict__ sequence, //not transposed
            int sequenceLength, 
            bool isForward,
            int columnStart,
            float overlapweight,
            const char* __restrict__ quality, //not transposed
            bool canUseQualityScores,
            bool isPairedCandidate
        ){

            auto defaultweightfunc = [&](char q){
                return canUseQualityScores ? getQualityWeight(q) * overlapweight : overlapweight;
            };
            auto weightfunc = defaultweightfunc;

            for(int p = group.thread_rank(); p < sequenceLength; p += group.size()){
                std::int8_t encodedBaseAsInt = SequenceHelpers::getEncodedNuc2Bit(sequence,sequenceLength,p);
                if(!isForward){
                    encodedBaseAsInt = SequenceHelpers::complementBase2Bit(encodedBaseAsInt);
                }
                const char qual = quality[p];
                const float weight = weightfunc(qual);

                assert(weight != 0);
                const int rowOffset = encodedBaseAsInt * columnPitchInElements;
                const int columnIndex = columnStart + (isForward ? p : sequenceLength - 1 - p);
                if(doAdd){
                    counts[rowOffset + columnIndex] += 1;
                    weights[rowOffset + columnIndex] += weight;
                    coverages[columnIndex] += 1;
                }else{
                    counts[rowOffset + columnIndex] -= 1;
                    weights[rowOffset + columnIndex] -= weight;
                    coverages[columnIndex] -= 1;
                }
            }
            group.sync();
        }

        template<class ThreadGroup>
        __device__ __forceinline__
        void setWeightsToZeroIfCountIsZero(ThreadGroup& group){
            const int firstColumn_incl = columnProperties->firstColumn_incl;
            const int lastColumn_excl = columnProperties->lastColumn_excl;
    
            for(int column = firstColumn_incl + group.thread_rank(); 
                    column < lastColumn_excl; 
                    column += group.size()){
    
                #pragma unroll
                for(int k = 0; k < 4; k++){
                    const int count = counts[k * columnPitchInElements + column];
                    if(count == 0){
                        weights[k * columnPitchInElements + column] = 0;
                    }
                }
            }
            group.sync();
        }

        template<class ThreadGroup>
        __device__ __forceinline__
        void constructFromSequences(
            ThreadGroup& group,
            const int* __restrict__ myShifts,
            const int* __restrict__ myOverlaps,
            const int* __restrict__ myNops,
            const AlignmentOrientation* __restrict__ myAlignmentFlags,
            const unsigned int* __restrict__ myAnchorSequenceData,
            const char* __restrict__ myAnchorQualityData,
            const unsigned int* __restrict__ myCandidateSequencesData,
            const char* __restrict__ myCandidateQualities,
            const int* __restrict__ myCandidateLengths,
            const bool* __restrict__ myPairedCandidateFlags,
            const int* __restrict__ myIndices,
            int numIndices,
            bool canUseQualityScores, 
            size_t encodedSequencePitchInInts,
            size_t qualityPitchInBytes,
            float desiredAlignmentMaxErrorRate,
            int anchorIndex
        ){   
            for(int column = group.thread_rank(); column < columnPitchInElements; column += group.size()){
                for(int i = 0; i < 4; i++){
                    counts[i * columnPitchInElements + column] = 0;
                    weights[i * columnPitchInElements + column] = 0.0f;
                }

                coverages[column] = 0;
            }
            
            group.sync();
            
            const int anchorColumnsBegin_incl = columnProperties->anchorColumnsBegin_incl;
            const int anchorColumnsEnd_excl = columnProperties->anchorColumnsEnd_excl;

            const int anchorLength = anchorColumnsEnd_excl - anchorColumnsBegin_incl;
            const unsigned int* const anchor = myAnchorSequenceData;
            const char* const anchorQualityScore = myAnchorQualityData;
                            
            //add anchor
            msaAddOrDeleteASequence2BitNormal<true>(
                group,
                anchor, 
                anchorLength, 
                true,
                anchorColumnsBegin_incl,
                1.0f,
                anchorQualityScore,
                canUseQualityScores,
                false
            );

            //add candidates
            for(int indexInList = 0; indexInList < numIndices; indexInList += 1){

                const int localCandidateIndex = myIndices[indexInList];
                assert(localCandidateIndex != -1);
                const int shift = myShifts[localCandidateIndex];
                const AlignmentOrientation flag = myAlignmentFlags[localCandidateIndex];

                const int queryLength = myCandidateLengths[localCandidateIndex];
                const unsigned int* const query = myCandidateSequencesData 
                    + std::size_t(localCandidateIndex) * encodedSequencePitchInInts;

                const char* const queryQualityScore = myCandidateQualities 
                    + std::size_t(localCandidateIndex) * qualityPitchInBytes;

                const int query_alignment_overlap = myOverlaps[localCandidateIndex];
                const int query_alignment_nops = myNops[localCandidateIndex];

                const float overlapweight = calculateOverlapWeight(
                    anchorLength, 
                    query_alignment_nops, 
                    query_alignment_overlap,
                    desiredAlignmentMaxErrorRate
                );

                assert(overlapweight <= 1.0f);
                assert(overlapweight >= 0.0f);
                assert(flag != AlignmentOrientation::None); // indices should only be pointing to valid alignments

                const int defaultcolumnoffset = anchorColumnsBegin_incl + shift;

                const bool isForward = flag == AlignmentOrientation::Forward;

                msaAddOrDeleteASequence2BitNormal<true>(
                    group,
                    query, 
                    queryLength, 
                    isForward,
                    defaultcolumnoffset,
                    overlapweight,
                    queryQualityScore,
                    canUseQualityScores,
                    myPairedCandidateFlags[localCandidateIndex]
                );
            }
        }

        template<class ThreadGroup, class Selector>
        __device__ __forceinline__
        void removeCandidates(
            ThreadGroup& group,
            Selector shouldBeRemoved, //remove candidate myIndices[i] if shouldBeRemoved(i) == true
            const int* __restrict__ myShifts,
            const int* __restrict__ myOverlaps,
            const int* __restrict__ myNops,
            const AlignmentOrientation* __restrict__ myAlignmentFlags,
            const unsigned int* __restrict__ myCandidateSequencesData, //not transposed
            const char* __restrict__ myCandidateQualities, //not transposed
            const int* __restrict__ myCandidateLengths,
            const int* __restrict__ myIndices,
            const bool* __restrict__ myPairedCandidateFlags,
            int numIndices,
            bool canUseQualityScores, 
            size_t encodedSequencePitchInInts,
            size_t qualityPitchInBytes,
            float desiredAlignmentMaxErrorRate
        ){  
    

            const int anchorColumnsBegin_incl = columnProperties->anchorColumnsBegin_incl;
            const int anchorColumnsEnd_excl = columnProperties->anchorColumnsEnd_excl;

            for(int indexInList = 0; indexInList < numIndices; indexInList += 1){

                if(shouldBeRemoved(indexInList)){

                    const int localCandidateIndex = myIndices[indexInList];
                    assert(localCandidateIndex != -1);
                    const int shift = myShifts[localCandidateIndex];
                    const AlignmentOrientation flag = myAlignmentFlags[localCandidateIndex];

                    const int anchorLength = anchorColumnsEnd_excl - anchorColumnsBegin_incl;
                    const int queryLength = myCandidateLengths[localCandidateIndex];
                    const unsigned int* const query = myCandidateSequencesData + localCandidateIndex * encodedSequencePitchInInts;

                    const char* const queryQualityScore = myCandidateQualities + std::size_t(localCandidateIndex) * qualityPitchInBytes;

                    const int query_alignment_overlap = myOverlaps[localCandidateIndex];
                    const int query_alignment_nops = myNops[localCandidateIndex];

                    const float overlapweight = calculateOverlapWeight(
                        anchorLength, 
                        query_alignment_nops, 
                        query_alignment_overlap,
                        desiredAlignmentMaxErrorRate
                    );

                    assert(overlapweight <= 1.0f);
                    assert(overlapweight >= 0.0f);
                    assert(flag != AlignmentOrientation::None); // indices should only be pointing to valid alignments

                    const int defaultcolumnoffset = anchorColumnsBegin_incl + shift;

                    const bool isForward = flag == AlignmentOrientation::Forward;

                    msaAddOrDeleteASequence2BitNormal<false>(
                        group,
                        query, 
                        queryLength, 
                        isForward,
                        defaultcolumnoffset,
                        overlapweight,
                        queryQualityScore,
                        canUseQualityScores,
                        myPairedCandidateFlags[localCandidateIndex]
                    );
                }
            }
        }

        template<class ThreadGroup>
        __device__ __forceinline__
        void findConsensus_impl_char4(
            ThreadGroup& group,
            const unsigned int* __restrict__ myAnchorSequenceData, 
            int encodedSequencePitchInInts,
            int anchorIndex = -1
        ){

            const int anchorColumnsBegin_incl = columnProperties->anchorColumnsBegin_incl;
            const int anchorColumnsEnd_excl = columnProperties->anchorColumnsEnd_excl;
            const int firstColumn_incl = columnProperties->firstColumn_incl;
            const int lastColumn_excl = columnProperties->lastColumn_excl;

            assert(lastColumn_excl <= columnPitchInElements);

            const int anchorLength = anchorColumnsEnd_excl - anchorColumnsBegin_incl;
            const unsigned int* const anchor = myAnchorSequenceData;

            //set columns to zero which are not in range firstColumn_incl <= column && column < lastColumn_excl

            for(int column = group.thread_rank(); 
                    column < firstColumn_incl; 
                    column += group.size()){

                support[column] = 0;
                origWeights[column] = 0;
                origCoverages[column] = 0;
            }

            const int leftoverRight = columnPitchInElements - lastColumn_excl;

            for(int i = group.thread_rank(); i < leftoverRight; i += group.size()){
                const int column = lastColumn_excl + i;

                support[column] = 0;
                origWeights[column] = 0;
                origCoverages[column] = 0;
            }

            for(int column = group.thread_rank(); 
                    column < firstColumn_incl; 
                    column += group.size()){
                    
                consensus[column] = std::uint8_t{5};
            }

            for(int i = group.thread_rank(); i < leftoverRight; i += group.size()){
                const int column = lastColumn_excl + i;

                consensus[column] = std::uint8_t{5};
            }

            const int* const myCountsA = counts + 0 * columnPitchInElements;
            const int* const myCountsC = counts + 1 * columnPitchInElements;
            const int* const myCountsG = counts + 2 * columnPitchInElements;
            const int* const myCountsT = counts + 3 * columnPitchInElements;

            const float* const myWeightsA = weights + 0 * columnPitchInElements;
            const float* const myWeightsC = weights + 1 * columnPitchInElements;
            const float* const myWeightsG = weights + 2 * columnPitchInElements;
            const float* const myWeightsT = weights + 3 * columnPitchInElements;

            const int numOuterIters = SDIV(lastColumn_excl, 4);

            for(int outerIter = group.thread_rank(); outerIter < numOuterIters; outerIter += group.size()){

                alignas(4) std::uint8_t consensusArray[4];

                #pragma unroll 
                for(int i = 0; i < 4; i++){
                    const int column = outerIter * 4 + i;

                    if(firstColumn_incl <= column && column < lastColumn_excl){

                        const int ca = myCountsA[column];
                        const int cc = myCountsC[column];
                        const int cg = myCountsG[column];
                        const int ct = myCountsT[column];
                        const float wa = myWeightsA[column];
                        const float wc = myWeightsC[column];
                        const float wg = myWeightsG[column];
                        const float wt = myWeightsT[column];

                        std::uint8_t cons = 5;
                        float consWeight = 0.0f;
                        if(wa > consWeight){
                            cons = std::uint8_t{0};
                            consWeight = wa;
                        }
                        if(wc > consWeight){
                            cons = std::uint8_t{1};
                            consWeight = wc;
                        }
                        if(wg > consWeight){
                            cons = std::uint8_t{2};
                            consWeight = wg;
                        }
                        if(wt > consWeight){
                            cons = std::uint8_t{3};
                            consWeight = wt;
                        }

                        consensusArray[i] = cons;

                        const float columnWeight = wa + wc + wg + wt;
                        assert(columnWeight != 0);

                        support[column] = consWeight / columnWeight;

                        if(anchorColumnsBegin_incl <= column && column < anchorColumnsEnd_excl){
                            constexpr std::uint8_t A_enc = SequenceHelpers::encodedbaseA();
                            constexpr std::uint8_t C_enc = SequenceHelpers::encodedbaseC();
                            constexpr std::uint8_t G_enc = SequenceHelpers::encodedbaseG();
                            constexpr std::uint8_t T_enc = SequenceHelpers::encodedbaseT();

                            const int localIndex = column - anchorColumnsBegin_incl;
                            const std::uint8_t encNuc = SequenceHelpers::getEncodedNuc2Bit(anchor, anchorLength, localIndex);

                            if(encNuc == A_enc){
                                origWeights[column] = wa;
                                origCoverages[column] = ca;
                            }else if(encNuc == C_enc){
                                origWeights[column] = wc;
                                origCoverages[column] = cc;
                            }else if(encNuc == G_enc){
                                origWeights[column] = wg;
                                origCoverages[column] = cg;
                            }else if(encNuc == T_enc){
                                origWeights[column] = wt;
                                origCoverages[column] = ct;
                            }
                        }
                    }
                }

                *((char4*)(consensus + 4*outerIter)) = *((const char4*)(&consensusArray[0]));
            }
        }

        template<class ThreadGroup>
        __device__ __forceinline__
        void findConsensus_impl_char_gmem(
            ThreadGroup& group,
            const unsigned int* __restrict__ myAnchorSequenceData, 
            int encodedSequencePitchInInts,
            int anchorIndex = -1
        ){

            const int anchorColumnsBegin_incl = columnProperties->anchorColumnsBegin_incl;
            const int anchorColumnsEnd_excl = columnProperties->anchorColumnsEnd_excl;
            const int firstColumn_incl = columnProperties->firstColumn_incl;
            const int lastColumn_excl = columnProperties->lastColumn_excl;

            assert(lastColumn_excl <= columnPitchInElements);

            const int anchorLength = anchorColumnsEnd_excl - anchorColumnsBegin_incl;
            const unsigned int* const anchor = myAnchorSequenceData;

            //set columns to zero which are not in range firstColumn_incl <= column && column < lastColumn_excl

            for(int column = group.thread_rank(); 
                    column < firstColumn_incl; 
                    column += group.size()){

                support[column] = 0;
                origWeights[column] = 0;
                origCoverages[column] = 0;
            }

            const int leftoverRight = columnPitchInElements - lastColumn_excl;

            for(int i = group.thread_rank(); i < leftoverRight; i += group.size()){
                const int column = lastColumn_excl + i;

                support[column] = 0;
                origWeights[column] = 0;
                origCoverages[column] = 0;
            }

            for(int column = group.thread_rank(); 
                    column < firstColumn_incl; 
                    column += group.size()){
                    
                consensus[column] = std::uint8_t{5};
            }

            for(int i = group.thread_rank(); i < leftoverRight; i += group.size()){
                const int column = lastColumn_excl + i;

                consensus[column] = std::uint8_t{5};
            }

            const int* const myCountsA = counts + 0 * columnPitchInElements;
            const int* const myCountsC = counts + 1 * columnPitchInElements;
            const int* const myCountsG = counts + 2 * columnPitchInElements;
            const int* const myCountsT = counts + 3 * columnPitchInElements;

            const float* const myWeightsA = weights + 0 * columnPitchInElements;
            const float* const myWeightsC = weights + 1 * columnPitchInElements;
            const float* const myWeightsG = weights + 2 * columnPitchInElements;
            const float* const myWeightsT = weights + 3 * columnPitchInElements;

            const int numOuterIters = SDIV(lastColumn_excl, group.size());

            for(int outerIter = 0; outerIter < numOuterIters; outerIter += 1){
                const int begin = outerIter * group.size();
                const int num = std::min(int(group.size()), lastColumn_excl - begin);

                if(group.thread_rank() < num){
                    const int column = begin + group.thread_rank();
                    if(firstColumn_incl <= column && column < lastColumn_excl){

                        const float wa = myWeightsA[column];
                        const float wc = myWeightsC[column];
                        const float wg = myWeightsG[column];
                        const float wt = myWeightsT[column];

                        std::uint8_t cons = 5;
                        float consWeight = 0.0f;
                        if(wa > consWeight){
                            cons = std::uint8_t{0};
                            consWeight = wa;
                        }
                        if(wc > consWeight){
                            cons = std::uint8_t{1};
                            consWeight = wc;
                        }
                        if(wg > consWeight){
                            cons = std::uint8_t{2};
                            consWeight = wg;
                        }
                        if(wt > consWeight){
                            cons = std::uint8_t{3};
                            consWeight = wt;
                        }

                        consensus[column] = cons;

                        const float columnWeight = wa + wc + wg + wt;
                        assert(columnWeight != 0);

                        support[column] = consWeight / columnWeight;

                        if(anchorColumnsBegin_incl <= column && column < anchorColumnsEnd_excl){
                            constexpr std::uint8_t A_enc = SequenceHelpers::encodedbaseA();
                            constexpr std::uint8_t C_enc = SequenceHelpers::encodedbaseC();
                            constexpr std::uint8_t G_enc = SequenceHelpers::encodedbaseG();
                            constexpr std::uint8_t T_enc = SequenceHelpers::encodedbaseT();

                            const int localIndex = column - anchorColumnsBegin_incl;
                            const std::uint8_t encNuc = SequenceHelpers::getEncodedNuc2Bit(anchor, anchorLength, localIndex);

                            if(encNuc == A_enc){
                                origWeights[column] = wa;
                                origCoverages[column] = myCountsA[column];
                            }else if(encNuc == C_enc){
                                origWeights[column] = wc;
                                origCoverages[column] = myCountsC[column];
                            }else if(encNuc == G_enc){
                                origWeights[column] = wg;
                                origCoverages[column] = myCountsG[column];
                            }else if(encNuc == T_enc){
                                origWeights[column] = wt;
                                origCoverages[column] = myCountsT[column];
                            }
                        }
                    }
                }
            }
        }

        template<class ThreadGroup>
        __device__ __forceinline__
        void findConsensus(
            ThreadGroup& group,
            const unsigned int* __restrict__ myAnchorSequenceData, 
            int encodedSequencePitchInInts,
            int anchorIndex = -1
        ){

            // findConsensus_impl_char4(
            //     group,
            //     myAnchorSequenceData,
            //     encodedSequencePitchInInts,
            //     anchorIndex
            // );

            findConsensus_impl_char_gmem(
                group,
                myAnchorSequenceData,
                encodedSequencePitchInInts,
                anchorIndex
            );
        }

        template<
            class ThreadGroup, 
            class GroupReduceBool, 
            class GroupReduceInt2
        >
        __device__ __forceinline__
        int flagCandidatesOfDifferentRegion(
            ThreadGroup& group,
            GroupReduceBool& groupReduceBool,
            GroupReduceInt2& groupReduceInt2,
            int* __restrict__ myNewIndicesPtr,
            const unsigned int* __restrict__ myAnchorSequenceData,
            const int anchorLength,
            const unsigned int* __restrict__ myCandidateSequencesData,
            const int* __restrict__ myCandidateLengths,
            const AlignmentOrientation* myAlignmentFlags,
            const int* __restrict__ myShifts,
            const int* __restrict__ myNops,
            const int* __restrict__ myOverlaps,
            bool* __restrict__ myShouldBeKept,
            float desiredAlignmentMaxErrorRate,
            int anchorIndex,
            int encodedSequencePitchInInts,
            const int* __restrict__ myIndices,
            const int myNumIndices,
            int dataset_coverage
        ){
            
            static_assert(std::is_same<ThreadGroup, cg::thread_block>::value, 
                "Can only use thread_block for msa::flagCandidatesOfDifferentRegion"); //because of shared per group               

            auto is_significant_count = [](int count, int coverage){
                if(int(coverage * 0.3f) <= count)
                    return true;
                return false;
            };        

            __shared__ bool broadcastbufferbool;
            __shared__ int broadcastbufferint4[4];
            __shared__ int smemcounts[1];

            const unsigned int* const anchorptr = myAnchorSequenceData;

            const int anchorColumnsBegin_incl = columnProperties->anchorColumnsBegin_incl;
            const int anchorColumnsEnd_excl = columnProperties->anchorColumnsEnd_excl;

            //check if anchor and consensus differ at at least one position

            bool hasMismatchToConsensus = false;

            for(int pos = group.thread_rank(); pos < anchorLength && !hasMismatchToConsensus; pos += group.size()){
                const int column = anchorColumnsBegin_incl + pos;
                const std::uint8_t consbase = consensus[column];
                const std::uint8_t anchorbase = SequenceHelpers::getEncodedNuc2Bit(anchorptr, anchorLength, pos);

                hasMismatchToConsensus |= (consbase != anchorbase);
            }

            hasMismatchToConsensus = groupReduceBool(hasMismatchToConsensus, [](auto l, auto r){
                return l || r;
            });

            if(group.thread_rank() == 0){
                broadcastbufferbool = hasMismatchToConsensus;
            }
            group.sync();

            hasMismatchToConsensus = broadcastbufferbool;

            //if anchor and consensus differ at at least one position, check columns in msa

            if(hasMismatchToConsensus){
                int col = std::numeric_limits<int>::max();
                bool foundColumn = false;
                std::uint8_t foundBase = 5;

                const int* const myCountsA = counts + 0 * columnPitchInElements;
                const int* const myCountsC = counts + 1 * columnPitchInElements;
                const int* const myCountsG = counts + 2 * columnPitchInElements;
                const int* const myCountsT = counts + 3 * columnPitchInElements;

                for(int columnindex = anchorColumnsBegin_incl + group.thread_rank(); 
                        columnindex < anchorColumnsEnd_excl && !foundColumn; 
                        columnindex += group.size()){

                    int regcounts[4];
                    regcounts[0] = myCountsA[columnindex];
                    regcounts[1] = myCountsC[columnindex];
                    regcounts[2] = myCountsG[columnindex];
                    regcounts[3] = myCountsT[columnindex];

                    const std::uint8_t consbase = consensus[columnindex];

                    assert(consbase < 4);

                    //find out if there is a non-consensus base with significant coverage
                    int significantBaseIndex = -1;

                    #pragma unroll
                    for(int i = 0; i < 4; i++){
                        if(i != consbase){
                            const bool significant = is_significant_count(regcounts[i], dataset_coverage);

                            significantBaseIndex = significant ? i : significantBaseIndex;
                        }
                    }

                    if(significantBaseIndex != -1){
                        foundColumn = true;
                        col = columnindex;
                        foundBase = significantBaseIndex;
                    }
                }

                int2 packed{col, foundBase};
                //find packed value with smallest col
                packed = groupReduceInt2(packed, [](auto l, auto r){
                    if(l.x < r.x){
                        return l;
                    }else{
                        return r;
                    }
                });

                if(group.thread_rank() == 0){
                    if(packed.x != std::numeric_limits<int>::max()){
                        broadcastbufferint4[0] = 1;
                        broadcastbufferint4[1] = packed.x;
                        broadcastbufferint4[2] = packed.y;
                    }else{
                        broadcastbufferint4[0] = 0;
                    }
                }

                group.sync();

                foundColumn = (1 == broadcastbufferint4[0]);
                col = broadcastbufferint4[1];
                foundBase = broadcastbufferint4[2];

                if(foundColumn){
                    
                    auto discard_rows = [&](bool keepMatching){
                        bool veryGoodAlignment = false;

                        for(int k = group.thread_rank(); k < myNumIndices; k += group.size()){
                            const int localCandidateIndex = myIndices[k];
                            const unsigned int* const candidateptr = myCandidateSequencesData + std::size_t(localCandidateIndex) * encodedSequencePitchInInts;
                            const int candidateLength = myCandidateLengths[localCandidateIndex];
                            const int shift = myShifts[localCandidateIndex];
                            const AlignmentOrientation alignmentFlag = myAlignmentFlags[localCandidateIndex];

                            //check if row is affected by column col
                            const int row_begin_incl = anchorColumnsBegin_incl + shift;
                            const int row_end_excl = row_begin_incl + candidateLength;
                            const bool notAffected = (col < row_begin_incl || row_end_excl <= col);
                            std::uint8_t base = 5;
                            if(!notAffected){
                                if(alignmentFlag == AlignmentOrientation::Forward){
                                    base = SequenceHelpers::getEncodedNuc2Bit(candidateptr, candidateLength, (col - row_begin_incl));
                                }else{
                                    //candidates cannot have AlignmentOrientation::None
                                    assert(alignmentFlag == AlignmentOrientation::ReverseComplement); 

                                    const std::uint8_t forwardbaseEncoded = SequenceHelpers::getEncodedNuc2Bit(candidateptr, candidateLength, row_end_excl-1 - col);
                                    base = SequenceHelpers::complementBase2Bit(forwardbaseEncoded);
                                }
                            }

                            bool shouldBeKept = false;
                            if(notAffected){
                                shouldBeKept = true;
                            }else if(keepMatching && (base == foundBase)){
                                //keep candidates which match the found base
                                shouldBeKept = true;
                            }else if(!keepMatching && (base != foundBase)){
                                //keep candidates which do not match the found base
                                shouldBeKept = true;
                            }

                            myShouldBeKept[k] = shouldBeKept;

                            if(!shouldBeKept && !veryGoodAlignment){
                                const int nOps = myNops[localCandidateIndex];
                                const int overlapsize = myOverlaps[localCandidateIndex];
                                const float overlapweight = calculateOverlapWeight(
                                    anchorLength, 
                                    nOps, 
                                    overlapsize,
                                    desiredAlignmentMaxErrorRate
                                );

                                assert(overlapweight <= 1.0f);
                                assert(overlapweight >= 0.0f);

                                if(fgeq(overlapweight, 0.90f)){
                                    veryGoodAlignment = true;
                                }
                            }
                        }
                        #if 1
                        //check that no candidate which should be removed has very good alignment.
                        //if there is such a candidate, none of the candidates will be removed.

                        veryGoodAlignment = groupReduceBool(veryGoodAlignment, [](auto l, auto r){
                            return l || r;
                        });

                        if(group.thread_rank() == 0){
                            broadcastbufferbool = veryGoodAlignment;
                        }
                        group.sync();

                        veryGoodAlignment = broadcastbufferbool;

                        if(veryGoodAlignment){
                            for(int k = group.thread_rank(); k < myNumIndices; k += group.size()){
                                myShouldBeKept[k] = true;
                            }
                        }
                        #endif

                        //select indices of candidates to keep and write them to new indices
                        if(group.thread_rank() == 0){
                            smemcounts[0] = 0;
                        }
                        group.sync();

                        assert(group.size() <= 256);

                        const int limit = SDIV(myNumIndices, group.size()) * group.size();
                        for(int k = group.thread_rank(); k < limit; k += group.size()){
                            bool keep = false;
                            if(k < myNumIndices){
                                keep = myShouldBeKept[k];
                            }                               
                
                            //warp time-sliced compaction to make it a stable compaction
                            for(int x = 0; x < group.size() / 32; x++){
                                if(group.thread_rank() / 32 == x){
                                    if(keep){
                                        cg::coalesced_group g = cg::coalesced_threads();
                                        int outputPos;
                                        if (g.thread_rank() == 0) {
                                            outputPos = atomicAdd(&smemcounts[0], g.size());
                                        }
                                        outputPos = g.thread_rank() + g.shfl(outputPos, 0);
                                        myNewIndicesPtr[outputPos] = myIndices[k];
                                    }
                                }
                                group.sync();
                            }
                        }

                        group.sync();

                        return smemcounts[0];
                    };

                    //compare found base to original base
                    const std::uint8_t originalbase = SequenceHelpers::getEncodedNuc2Bit(anchorptr, anchorLength, col - anchorColumnsBegin_incl);

                    if(originalbase == foundBase){
                        //discard all candidates whose base in column col differs from foundBase
                        return discard_rows(true);
                    }else{
                        //discard all candidates whose base in column col matches foundBase
                        return discard_rows(false);
                    }

                }else{
                    //did not find a significant columns

                    for(int k = group.thread_rank(); k < myNumIndices; k += group.size()){
                        myShouldBeKept[k] = true;
                    }

                    //remove no candidate
                    for(int k = group.thread_rank(); k < myNumIndices; k += group.size()){
                        myNewIndicesPtr[k] = myIndices[k];
                    }
                    group.sync();

                    return myNumIndices;
                }

            }else{
                //no mismatch between consensus and anchor

                //remove no candidate
                for(int k = group.thread_rank(); k < myNumIndices; k += group.size()){
                    myShouldBeKept[k] = true;
                }

                for(int k = group.thread_rank(); k < myNumIndices; k += group.size()){
                    myNewIndicesPtr[k] = myIndices[k];
                }
                group.sync();

                return myNumIndices;
            }
        }

        __device__ __forceinline__
        auto getConsensusQualityIterator() const{
            return thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                [this](int pos){
                    return getConsensusQualityOfPosition(pos);
                }
            );
        }

        __device__ __forceinline__
        auto getDecodedConsensusIterator() const{
            return thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                [this](int pos){
                    return getDecodedConsensusOfPosition(pos);
                }
            );
        }


        template<class ThreadGroup>
        __device__ __forceinline__
        void computeConsensusQuality(
            ThreadGroup& group,
            char* quality,
            int maxlength
        ) const {
            const int begin = columnProperties->firstColumn_incl;
            const int end = columnProperties->lastColumn_excl;

            auto consensusQualityIterator = getConsensusQualityIterator();

            for(int i = begin + group.thread_rank(); i < end; i += group.size()){
                if(i - begin < maxlength){
                    const int outpos = i - begin;
                    quality[outpos] = consensusQualityIterator[i];
                }
            }
        }

        template<class ThreadGroup>
        __device__ __forceinline__
        void computeDecodedConsensus(
            ThreadGroup& group,
            char* decodedConsensus,
            int maxlength
        ) const {
            const int begin = columnProperties->firstColumn_incl;
            const int end = columnProperties->lastColumn_excl;

            auto decodedConsensusIterator = getDecodedConsensusIterator();

            for(int i = begin + group.thread_rank(); i < end; i += group.size()){
                if(i - begin < maxlength){
                    const int outpos = i - begin;
                    decodedConsensus[outpos] = decodedConsensusIterator[i];
                }
            }
        }

        __device__ __forceinline__
        int computeSize() const {
            const int begin = columnProperties->firstColumn_incl;
            const int end = columnProperties->lastColumn_excl;

            return end - begin;
        }

    private:
        __device__ __forceinline__
        char getConsensusQualityOfPosition(int pos) const{
            const float sup = support[pos];
            return getQualityChar(sup);
        }

        __device__ __forceinline__
        char getDecodedConsensusOfPosition(int pos) const{
            const std::uint8_t encoded = consensus[pos];
            char decoded = 'F';
            if(encoded == std::uint8_t{0}){
                decoded = 'A';
            }else if(encoded == std::uint8_t{1}){
                decoded = 'C';
            }else if(encoded == std::uint8_t{2}){
                decoded = 'G';
            }else if(encoded == std::uint8_t{3}){
                decoded = 'T';
            }
            return decoded;
        };

    public:
        int columnPitchInElements;
        std::uint8_t* consensus;
        int* counts;
        int* coverages;
        int* origCoverages;
        float* weights;
        float* support;
        float* origWeights;
        MSAColumnProperties* columnProperties;
    };

    struct GPUMultiMSA{
    public:
        HOSTDEVICEQUALIFIER
        GpuSingleMSA getSingleMSA(int msaIndex) const{
            GpuSingleMSA msa;

            msa.columnPitchInElements = columnPitchInElements;
            msa.counts = getCountsOfMSA(msaIndex);
            msa.weights = getWeightsOfMSA(msaIndex);
            msa.coverages = getCoveragesOfMSA(msaIndex);
            msa.consensus = getConsensusOfMSA(msaIndex);
            msa.support = getSupportOfMSA(msaIndex);
            msa.origWeights = getOrigWeightsOfMSA(msaIndex);
            msa.origCoverages = getOrigCoveragesOfMSA(msaIndex);
            msa.columnProperties = getColumnPropertiesOfMSA(msaIndex);

            return msa;
        }

        HOSTDEVICEQUALIFIER
        int* getCountsOfMSA(int msaIndex) const{
            return counts + std::size_t(columnPitchInElements) * 4 * msaIndex;
        }

        HOSTDEVICEQUALIFIER
        float* getWeightsOfMSA(int msaIndex) const{
            return weights + std::size_t(columnPitchInElements) * 4 * msaIndex;
        }

        HOSTDEVICEQUALIFIER
        int* getCoveragesOfMSA(int msaIndex) const{
            return coverages + std::size_t(columnPitchInElements) * msaIndex;
        }

        HOSTDEVICEQUALIFIER
        std::uint8_t* getConsensusOfMSA(int msaIndex) const{
            return consensus + std::size_t(columnPitchInElements) * msaIndex;
        }

        HOSTDEVICEQUALIFIER
        float* getSupportOfMSA(int msaIndex) const{
            return support + std::size_t(columnPitchInElements) * msaIndex;
        }

        HOSTDEVICEQUALIFIER
        float* getOrigWeightsOfMSA(int msaIndex) const{
            return origWeights + std::size_t(columnPitchInElements) * msaIndex;
        }

        HOSTDEVICEQUALIFIER
        int* getOrigCoveragesOfMSA(int msaIndex) const{
            return origCoverages + std::size_t(columnPitchInElements) * msaIndex;
        }

        HOSTDEVICEQUALIFIER
        MSAColumnProperties* getColumnPropertiesOfMSA(int msaIndex) const{
            return columnProperties + msaIndex;
        }

    public:
        int numMSAs;
        int columnPitchInElements;
        std::uint8_t* consensus;
        int* counts;
        int* coverages;
        int* origCoverages;
        float* weights;
        float* support;
        float* origWeights;
        MSAColumnProperties* columnProperties;
    };


    
}
}


#endif

#endif