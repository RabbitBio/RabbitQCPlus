#ifndef CARE_MSA_SPLITS_HPP
#define CARE_MSA_SPLITS_HPP

#include <hpc_helpers.cuh>
#include <hostdevicefunctions.cuh>
#include <msa.hpp>
#include <sequencehelpers.hpp>

#include <algorithm>
#include <limits>

#ifdef __CUDACC__
#include <cub/cub.cuh>
#endif

namespace care{


struct CheckAmbiguousColumns{

    const int* countsA;
    const int* countsC;
    const int* countsG;
    const int* countsT;
    const int* coverages;

    CheckAmbiguousColumns(const MultipleSequenceAlignment& msa)
        : countsA(msa.countsA.data()), countsC(msa.countsC.data()), countsG(msa.countsG.data()), countsT(msa.countsT.data()), coverages(msa.coverage.data()){}

    CheckAmbiguousColumns(const int* cA, const int* cC, const int* cG, const int* cT, const int* cov) 
        : countsA(cA), countsC(cC), countsG(cG), countsT(cT), coverages(cov){}

    struct SplitInfo{
        char nuc;
        int column;
        //float ratio;

        SplitInfo() = default;
        SplitInfo(char c, int n) : nuc(c), column(n){}
        SplitInfo(char c, int n, float /*f*/) : nuc(c), column(n)/*, ratio(f)*/{}
    };

    struct SplitInfos{       
        static constexpr int maxSplitInfos = 64;

        int numSplitInfos = 0;
        SplitInfo splitInfos[maxSplitInfos];
    };

    SplitInfos getSplitInfos(int begin, int end, float relativeCountLowerBound, float relativeCountUpperBound) const{
        SplitInfos result{};

        //find columns where counts/coverage is in range [relativeCountLowerBound, relativeCountUpperBound]
        //if there are exactly 2 different bases for which this is true in a column, SplitInfos of both cases are saved in result.splitInfos
        for(int col = begin; col < end; col++){

            SplitInfo myResults[2];
            int myNumResults = 0;

            auto checkNuc = [&](const auto& counts, const char nuc){
                const float ratio = float(counts[col]) / float(coverages[col]);
                if(counts[col] >= 2 && fgeq(ratio, relativeCountLowerBound) && fleq(ratio, relativeCountUpperBound)){
                    if(myNumResults < 2){
                        myResults[myNumResults] = {nuc, col, ratio};
                        myNumResults++;
                    }                       
                }
            };

            checkNuc(countsA, 'A');
            checkNuc(countsC, 'C');
            checkNuc(countsG, 'G');
            checkNuc(countsT, 'T');

            if(myNumResults == 2){
                if(result.numSplitInfos < SplitInfos::maxSplitInfos - 2){
                    result.splitInfos[result.numSplitInfos + 0] = myResults[0];
                    result.splitInfos[result.numSplitInfos + 1] = myResults[1];
                    result.numSplitInfos += 2;
                }else{
                    break;
                }
            }
        }

        return result;
    }

    int getNumberOfSplits(
        const SplitInfos& splitInfos, 
        const MultipleSequenceAlignment& msa
    ){

        // 64-bit flags. 2 bit per column. 
        // 00 -> nuc does not match any of both splitInfos
        // 10 -> nuc  matches first splitInfos
        // 11 -> nuc  matches second splitInfos
        constexpr int numPossibleColumnsPerFlag = 32; 
        //

        if(splitInfos.numSplitInfos == 0) return 1;
        if(splitInfos.numSplitInfos == 2) return 2;

        static constexpr int maxNumEncodedRows = 128;
        int numEncodedRows = 0;
        std::uint64_t encodedRows[maxNumEncodedRows];

        const int anchorColumnsBegin_incl = msa.anchorColumnsBegin_incl;
        const int numColumnsToCheck = std::min(numPossibleColumnsPerFlag, splitInfos.numSplitInfos / 2);
        const int maxCandidatesToCheck = msa.nCandidates;

        for(int c = 0; c < maxCandidatesToCheck; c += 1){
            if(numEncodedRows < maxNumEncodedRows){
                std::uint64_t flags = 0;

                const char* const myCandidate = msa.inputData.candidates + c * msa.inputData.candidatesPitch;
                const int candidateShift = msa.inputData.candidateShifts[c];
                const int candidateLength = msa.inputData.candidateLengths[c];

                for(int k = 0; k < numColumnsToCheck; k++){
                    flags <<= 2;

                    const SplitInfo psc0 = splitInfos.splitInfos[2*k+0];
                    const SplitInfo psc1 = splitInfos.splitInfos[2*k+1];
                    assert(psc0.column == psc1.column);

                    const int candidateColumnsBegin_incl = candidateShift + anchorColumnsBegin_incl;
                    const int candidateColumnsEnd_excl = candidateLength + candidateColumnsBegin_incl;
                    
                    //column range check for row
                    if(candidateColumnsBegin_incl <= psc0.column && psc0.column < candidateColumnsEnd_excl){                        
                        const int positionInCandidate = psc0.column - candidateColumnsBegin_incl;
                        const char nuc = myCandidate[positionInCandidate];

                        if(nuc == psc0.nuc){
                            flags = flags | 0b10;
                        }else if(nuc == psc1.nuc){
                            flags = flags | 0b11;
                        }else{
                            flags = flags | 0b00;
                        } 

                    }else{
                        flags = flags | 0b00;
                    } 

                }

                encodedRows[numEncodedRows] = flags;
                numEncodedRows++;
            }
        }


        //sort the computed encoded rows, and make them unique
        std::sort(&encodedRows[0], &encodedRows[numEncodedRows]);
        auto uniqueEnd = std::unique(&encodedRows[0], &encodedRows[numEncodedRows]);

        numEncodedRows = std::distance(&encodedRows[0], uniqueEnd);

        int numNoEqualEncodedRow = 1;
        for(int i = 0; i < numEncodedRows - 1; i++){
            const std::uint64_t encodedRow = encodedRows[i];

            std::uint64_t mask = 0;
            for(int s = 0; s < numColumnsToCheck; s++){
                if(encodedRow >> (2*s+1) & 1){
                    mask = mask | (0x03 << (2*s));
                }
            }

            bool any = false;

            for(int k = i+1; k < numEncodedRows; k++){
                //check if encodedRow is equal to another flag masked with mask. if yes, it can be removed
                if((encodedRows[k] & mask) == encodedRow){
                    any = true;
                    break;
                }
            }

            if(!any){
                numNoEqualEncodedRow++;
            }
        }

        return numNoEqualEncodedRow;
    }

};


namespace gpu{

#ifdef __CUDACC__

    template<int blocksize>
    struct CheckAmbiguousColumnsGpu{

        using BlockReduce = cub::BlockReduce<int, blocksize>;
        //using BlockScan = cub::BlockScan<int, blocksize>;
        using BlockSort = cub::BlockRadixSort<std::uint64_t, blocksize, 1>;
        using BlockDiscontinuity = cub::BlockDiscontinuity<std::uint64_t, blocksize>;
        using MyBlockSelect = BlockSelect<std::uint64_t, blocksize>;

        const int* countsA;
        const int* countsC;
        const int* countsG;
        const int* countsT;
        const int* coverages;

        void* tempstorage;

        __host__ __device__
        CheckAmbiguousColumnsGpu(const int* cA, const int* cC, const int* cG, const int* cT, const int* cov) 
            : countsA(cA), countsC(cC), countsG(cG), countsT(cT), coverages(cov){}

        //thread 0 returns number of ambiguous columns in given range. Block-wide algorithm
        __device__
        int getAmbiguousColumnCount(int begin, int end, typename BlockReduce::TempStorage& temp) const{ 

            int myCount = 0;

            for(int col = threadIdx.x; col < end; col += blocksize){
                if(col >= begin){

                    int numNucs = 0;

                    auto checkNuc = [&](const auto& counts, const char nuc){
                        const float ratio = float(counts[col]) / float(coverages[col]);
                        if(counts[col] >= 2 && fgeq(ratio, 0.4f) && fleq(ratio, 0.6f)){
                            numNucs++;                                
                        }
                    };

                    checkNuc(countsA, 'A');
                    checkNuc(countsC, 'C');
                    checkNuc(countsG, 'G');
                    checkNuc(countsT, 'T');

                    if(numNucs > 0){
                        myCount++;
                    }
                }
            }

            myCount = BlockReduce(temp).Sum(myCount);

            return myCount;
        }

        struct SplitInfo{
            char nuc;
            int column;
            //float ratio;

            SplitInfo() = default;

            __host__ __device__
            SplitInfo(char c, int n) : nuc(c), column(n){}

            __host__ __device__
            SplitInfo(char c, int n, float f) : nuc(c), column(n)/*, ratio(f)*/{}
        };

        struct SplitInfos{       
            static constexpr int maxSplitInfos = 64;

            int numSplitInfos;
            SplitInfo splitInfos[maxSplitInfos];
        };

        struct TempStorage{
            int broadcastint;
            union{
                typename BlockReduce::TempStorage blockreducetemp;
                //typename BlockScan::TempStorage blockscantemp;
                typename BlockSort::TempStorage blocksorttemp;
                typename BlockDiscontinuity::TempStorage blockdiscontinuity;
                typename MyBlockSelect::TempStorage blockselect;
            } cub;

            static constexpr int maxNumEncodedRows = 128;
            int numEncodedRows;
            std::uint64_t encodedRows[maxNumEncodedRows];
            int flags[maxNumEncodedRows];
        };



        __device__
        void getSplitInfos(int begin, int end, float relativeCountLowerBound, float relativeCountUpperBound, SplitInfos& result) const{
            using PSC = MultipleSequenceAlignment::PossibleSplitColumn;

            if(threadIdx.x == 0){
                result.numSplitInfos = 0;
            }
            __syncthreads();

            //find columns where counts/coverage is in range [relativeCountLowerBound, relativeCountUpperBound]
            //if there are exactly 2 different bases for which this is true in a column, SplitInfos of both cases are saved in result.splitInfos
            for(int col = threadIdx.x; col < end; col += blocksize){
                if(col >= begin){

                    SplitInfo myResults[2];
                    int myNumResults = 0;

                    auto checkNuc = [&](const auto& counts, const char nuc){
                        const float ratio = float(counts[col]) / float(coverages[col]);
                        if(counts[col] >= 2 && fgeq(ratio, relativeCountLowerBound) && fleq(ratio, relativeCountUpperBound)){

                            #pragma unroll
                            for(int k = 0; k < 2; k++){
                                if(myNumResults == k){
                                    myResults[k] = {nuc, col, ratio};
                                    myNumResults++;
                                    break;
                                }
                            }
                              
                        }
                    };

                    checkNuc(countsA, 'A');
                    checkNuc(countsC, 'C');
                    checkNuc(countsG, 'G');
                    checkNuc(countsT, 'T');

                    if(myNumResults == 2){
                        if(result.numSplitInfos < SplitInfos::maxSplitInfos){
                            int firstPos = atomicAdd(&result.numSplitInfos, 2);
                            if(firstPos <= SplitInfos::maxSplitInfos - 2){
                                result.splitInfos[firstPos + 0] = myResults[0];
                                result.splitInfos[firstPos + 1] = myResults[1];
                            }
                        }
                    }
                }
            }

            __syncthreads();
            if(threadIdx.x == 0){
                if(result.numSplitInfos > SplitInfos::maxSplitInfos){
                    result.numSplitInfos = SplitInfos::maxSplitInfos;
                }
            }
            __syncthreads();
        }

        __device__
        int getNumberOfSplits(
            const SplitInfos& splitInfos, 
            const gpu::MSAColumnProperties& msaColumnProperties,
            int numCandidates, 
            const int* candidateShifts,
            const AlignmentOrientation* bestAlignmentFlags,
            const int* candidateLengths,
            const unsigned int* encodedCandidates, 
            int encodedSequencePitchInInts, 
            TempStorage& temp
        ){
            assert(TempStorage::maxNumEncodedRows == blocksize);

            // 64-bit flags. 2 bit per column. 
            // 00 -> nuc does not match any of both splitInfos
            // 10 -> nuc  matches first splitInfos
            // 11 -> nuc  matches second splitInfos
            constexpr int numPossibleColumnsPerFlag = 32; 
            //

            if(splitInfos.numSplitInfos == 0) return 1;
            if(splitInfos.numSplitInfos == 2) return 2;

            if(threadIdx.x == 0){
                temp.numEncodedRows = 0;
            }
            __syncthreads();

            const int anchorColumnsBegin_incl = msaColumnProperties.anchorColumnsBegin_incl;
            const int numColumnsToCheck = std::min(numPossibleColumnsPerFlag, splitInfos.numSplitInfos / 2);
            const int maxCandidatesToCheck = std::min(blocksize, numCandidates);

            #if 0
            if(threadIdx.x == 0){                
                if(splitInfos.numSplitInfos > 0){
                    printf("numSplitInfos %d\n", splitInfos.numSplitInfos);
                    for(int i = 0; i < splitInfos.numSplitInfos; i++){
                        printf("(%c,%d, %f) ", 
                            splitInfos.splitInfos[i].nuc, 
                            splitInfos.splitInfos[i].column,
                            splitInfos.splitInfos[i].ratio);
                    }
                    printf("\n");

                    for(int c = 0; c < maxCandidatesToCheck; c++){
                        const unsigned int* const myCandidate = encodedCandidates + c * encodedSequencePitchInInts;
                        const int candidateShift = candidateShifts[c];
                        const int candidateLength = candidateLengths[c];
                        const AlignmentOrientation alignmentFlag = bestAlignmentFlags[c];

                        for(int i = 0; i < candidateShift; i++){
                            printf("0");
                        }
                        for(int i = 0; i < candidateLength; i++){
                            char nuc = 'F';
                            if(alignmentFlag == AlignmentOrientation::Forward){
                                const int positionInCandidate = i;
                                std::uint8_t encodedCandidateNuc = SequenceHelpers::getEncodedNuc2Bit(myCandidate, candidateLength, positionInCandidate);
                                nuc = SequenceHelpers::decodeBase(encodedCandidateNuc);
                            }else{
                                const int positionInCandidate = candidateLength - 1 - i;
                                std::uint8_t encodedCandidateNuc = SequenceHelpers::getEncodedNuc2Bit(myCandidate, candidateLength, positionInCandidate);
                                nuc = SequenceHelpers::complementBaseDecoded(SequenceHelpers::decodeBase(encodedCandidateNuc));
                            }
                            printf("%c", nuc);
                        }
                        printf("\n");
                    }
                }
            }

            __syncthreads(); //DEBUG
            #endif

            for(int c = threadIdx.x; c < maxCandidatesToCheck; c += blocksize){
                if(temp.numEncodedRows < TempStorage::maxNumEncodedRows){
                    std::uint64_t flags = 0;

                    const unsigned int* const myCandidate = encodedCandidates + c * encodedSequencePitchInInts;
                    const int candidateShift = candidateShifts[c];
                    const int candidateLength = candidateLengths[c];
                    const AlignmentOrientation alignmentFlag = bestAlignmentFlags[c];

                    for(int k = 0; k < numColumnsToCheck; k++){
                        flags <<= 2;

                        const SplitInfo psc0 = splitInfos.splitInfos[2*k+0];
                        const SplitInfo psc1 = splitInfos.splitInfos[2*k+1];
                        assert(psc0.column == psc1.column);

                        const int candidateColumnsBegin_incl = candidateShift + anchorColumnsBegin_incl;
                        const int candidateColumnsEnd_excl = candidateLength + candidateColumnsBegin_incl;
                        
                        //column range check for row
                        if(candidateColumnsBegin_incl <= psc0.column && psc0.column < candidateColumnsEnd_excl){                        

                            char nuc = 'F';
                            if(alignmentFlag == AlignmentOrientation::Forward){
                                const int positionInCandidate = psc0.column - candidateColumnsBegin_incl;
                                std::uint8_t encodedCandidateNuc = SequenceHelpers::getEncodedNuc2Bit(myCandidate, candidateLength, positionInCandidate);
                                nuc = SequenceHelpers::decodeBase(encodedCandidateNuc);
                            }else{
                                const int positionInCandidate = candidateLength - 1 - (psc0.column - candidateColumnsBegin_incl);
                                std::uint8_t encodedCandidateNuc = SequenceHelpers::getEncodedNuc2Bit(myCandidate, candidateLength, positionInCandidate);
                                nuc = SequenceHelpers::complementBaseDecoded(SequenceHelpers::decodeBase(encodedCandidateNuc));
                            }

                            //printf("cand %d col %d %c\n", c, psc0.column, nuc);

                            if(nuc == psc0.nuc){
                                flags = flags | 0b10;
                            }else if(nuc == psc1.nuc){
                                flags = flags | 0b11;
                            }else{
                                flags = flags | 0b00;
                            } 

                        }else{
                            flags = flags | 0b00;
                        } 

                    }

                    const int tempPos = atomicAdd(&temp.numEncodedRows, 1);
                    if(tempPos < TempStorage::maxNumEncodedRows){
                        temp.encodedRows[tempPos] = flags;
                    }
                }
            }
        
            __syncthreads();
            if(threadIdx.x == 0){
                if(temp.numEncodedRows > TempStorage::maxNumEncodedRows){
                    temp.numEncodedRows = TempStorage::maxNumEncodedRows;
                }
            }
            __syncthreads();

            {
                //sort the computed encoded rows, and make them unique

                std::uint64_t encodedRow[1];

                if(threadIdx.x < temp.numEncodedRows){
                    encodedRow[0] = temp.encodedRows[threadIdx.x];
                }else{
                    encodedRow[0] = std::numeric_limits<std::uint64_t>::max();
                }

                BlockSort(temp.cub.blocksorttemp).Sort(encodedRow);
                __syncthreads();

                if(threadIdx.x < temp.numEncodedRows){
                    temp.encodedRows[threadIdx.x] = encodedRow[0];
                    temp.flags[threadIdx.x] = 0;
                }

                __syncthreads();

                int headflag[1];

                if(threadIdx.x < temp.numEncodedRows){
                    encodedRow[0] = temp.encodedRows[threadIdx.x];
                }else{
                    encodedRow[0] = temp.encodedRows[temp.numEncodedRows-1];
                }

                BlockDiscontinuity(temp.cub.blockdiscontinuity).FlagHeads(headflag, encodedRow, cub::Inequality());

                __syncthreads();

                int numselected = MyBlockSelect(temp.cub.blockselect).Flagged(encodedRow, headflag, &temp.encodedRows[0], temp.numEncodedRows);

                __syncthreads();
                if(threadIdx.x == 0){
                    temp.numEncodedRows = numselected;
                }
                __syncthreads();
            }

            // if(threadIdx.x == 0){                
            //     if(splitInfos.numSplitInfos > 0){
            //         printf("numEncodedRows %d\n", temp.numEncodedRows);
            //         for(int i = 0; i < temp.numEncodedRows; i++){
            //             printf("%lu ", temp.encodedRows[i]);
            //         }
            //         printf("\n");
            //     }
            // }

            // __syncthreads(); //DEBUG

            for(int i = 0; i < temp.numEncodedRows - 1; i++){
                const std::uint64_t encodedRow = temp.encodedRows[i];

                std::uint64_t mask = 0;
                for(int s = 0; s < numColumnsToCheck; s++){
                    if(encodedRow >> (2*s+1) & 1){
                        mask = mask | (0x03 << (2*s));
                    }
                }
                
                if(i < threadIdx.x && threadIdx.x < temp.numEncodedRows){
                    //check if encodedRow is equal to another flag masked with mask. if yes, it can be removed
                    if((temp.encodedRows[threadIdx.x] & mask) == encodedRow){
                        // if(encodedRow == 172){
                        //     printf("i = %d, thread=%d temp.encodedRows[threadIdx.x] = %lu, mask = %lu", i, threadIdx.x, temp.encodedRows[threadIdx.x], mask)
                        // }
                        atomicAdd(&temp.flags[i], 1);
                    }
                }
            }

            __syncthreads();

            //count number of remaining row flags
            int count = 0;
            if(threadIdx.x < temp.numEncodedRows){
                count = temp.flags[threadIdx.x] == 0 ? 1 : 0;
            }
            count = BlockReduce(temp.cub.blockreducetemp).Sum(count);
            if(threadIdx.x == 0){
                temp.broadcastint = count;
                // if(splitInfos.numSplitInfos > 0){
                //     printf("count = %d\n", count);
                // }
            }
            __syncthreads();

            count = temp.broadcastint;

            return count;
        }

    };



#endif

}


}

#endif