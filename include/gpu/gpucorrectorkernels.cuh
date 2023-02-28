#ifndef CARE_GPUCORRECTORKERNELS_CUH
#define CARE_GPUCORRECTORKERNELS_CUH


#include <config.hpp>
#include <correctedsequence.hpp>

#include <gpu/groupmemcpy.cuh>
#include <gpu/cuda_block_select.cuh>
#include <thrust/binary_search.h>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

namespace care{
namespace gpu{
    
namespace gpucorrectorkernels{

    __global__
    void copyCorrectionInputDeviceData(
        int* __restrict__ output_numAnchors,
        int* __restrict__ output_numCandidates,
        read_number* __restrict__ output_anchor_read_ids,
        unsigned int* __restrict__ output_anchor_sequences_data,
        int* __restrict__ output_anchor_sequences_lengths,
        read_number* __restrict__ output_candidate_read_ids,
        int* __restrict__ output_candidates_per_anchor,
        int* __restrict__ output_candidates_per_anchor_prefixsum,
        const int encodedSequencePitchInInts,
        const int input_numAnchors,
        const int input_numCandidates,
        const read_number* __restrict__ input_anchor_read_ids,
        const unsigned int* __restrict__ input_anchor_sequences_data,
        const int* __restrict__ input_anchor_sequences_lengths,
        const read_number* __restrict__ input_candidate_read_ids,
        const int* __restrict__ input_candidates_per_anchor,
        const int* __restrict__ input_candidates_per_anchor_prefixsum
    );

    __global__ 
    void copyMinhashResultsKernel(
        int* __restrict__ d_numCandidates,
        int* __restrict__ h_numCandidates,
        read_number* __restrict__ h_candidate_read_ids,
        const int* __restrict__ d_candidates_per_anchor_prefixsum,
        const read_number* __restrict__ d_candidate_read_ids,
        const int numAnchors
    );

    __global__
    void setAnchorIndicesOfCandidateskernel(
        int* __restrict__ d_anchorIndicesOfCandidates,
        const int* __restrict__ numAnchorsPtr,
        const int* __restrict__ d_candidates_per_anchor,
        const int* __restrict__ d_candidates_per_anchor_prefixsum
    );


    template<int blocksize, class Flags>
    __global__
    void selectIndicesOfFlagsOneBlock(
        int* __restrict__ selectedIndices,
        int* __restrict__ numSelectedIndices,
        const Flags flags,
        const int* __restrict__ numFlagsPtr
    ){
        constexpr int ITEMS_PER_THREAD = 4;
        constexpr int itemsPerIteration = blocksize * ITEMS_PER_THREAD;

        using MyBlockSelect = BlockSelect<int, blocksize>;

        __shared__ typename MyBlockSelect::TempStorage temp_storage;

        int aggregate = 0;
        const int numFlags = *numFlagsPtr;
        const int iters = SDIV(numFlags, blocksize * ITEMS_PER_THREAD);
        const int threadoffset = ITEMS_PER_THREAD * threadIdx.x;

        int remainingItems = numFlags;

        for(int iter = 0; iter < iters; iter++){
            const int validItems = min(remainingItems, itemsPerIteration);

            int data[ITEMS_PER_THREAD];

            const int iteroffset = itemsPerIteration * iter;

            #pragma unroll
            for(int k = 0; k < ITEMS_PER_THREAD; k++){
                if(iteroffset + threadoffset + k < numFlags){
                    data[k] = int(flags[iteroffset + threadoffset + k]);
                }else{
                    data[k] = 0;
                }
            }

            #pragma unroll
            for(int k = 0; k < ITEMS_PER_THREAD; k++){
                if(iteroffset + threadoffset + k < numFlags){
                    data[k] = data[k] != 0 ? 1 : 0;
                }
            }

            const int numSelected = MyBlockSelect(temp_storage).ForEachFlaggedPosition(data, validItems,
                [&](const auto& flaggedPosition, const int& outputpos){
                    selectedIndices[aggregate + outputpos] = iteroffset + flaggedPosition;
                }
            );

            aggregate += numSelected;
            remainingItems -= validItems;

            __syncthreads();
        }

        if(threadIdx.x == 0){
            *numSelectedIndices = aggregate;

            // for(int i = 0; i < aggregate; i++){
            //     printf("%d ", selectedIndices[i]);
            // }
            // printf("\n");
        }

    }

    __global__ 
    void initArraysBeforeCandidateCorrectionKernel(
        int maxNumCandidates,
        const int* __restrict__ d_numAnchors,
        int* __restrict__ d_num_corrected_candidates_per_anchor,
        bool* __restrict__ d_candidateCanBeCorrected
    );

    __global__
    void compactEditsKernel(
        const care::EncodedCorrectionEdit* __restrict__ d_inputEdits,
        care::EncodedCorrectionEdit* __restrict__ d_outputEdits,
        const int* __restrict__ d_editsOutputOffsets,
        const int* __restrict__ d_numSequences,
        const int* __restrict__ d_numEditsPerSequence,
        int doNotUseEditsValue,
        std::size_t editsPitchInBytes
    );

    template<int groupsize, class IndexIterator>
    __global__
    void compactCorrectedSequencesKernel(
        const char* __restrict__ d_inputCorrectedSequences,
        char* __restrict__ d_outputCorrectedSequences,
        std::size_t decodedSequencePitchInBytes,
        const int* __restrict__ d_numSequences,
        const int* __restrict__ d_numEditsPerSequence,
        int doNotUseEditsValue,
        const int* __restrict__ d_outputOffsets,
        IndexIterator d_indicesOfCorrectedSequences
    ){
        const int N = *d_numSequences;

        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        const int numGroups = (blockDim.x * gridDim.x) / groupsize;

        for(int c = groupId; c < N; c += numGroups){
            const int indexOfCorrectedSequence = d_indicesOfCorrectedSequences[c];
            const int numEdits = d_numEditsPerSequence[c];

            if(numEdits == doNotUseEditsValue){
                const int outputOffset = c == 0 ? 0 : d_outputOffsets[c];

                char* outputPtr = d_outputCorrectedSequences + outputOffset;
                const char* inputPtr = d_inputCorrectedSequences + indexOfCorrectedSequence * decodedSequencePitchInBytes;

                care::gpu::memcpy<int4>(group, outputPtr, inputPtr, decodedSequencePitchInBytes);
            }
        }
    }

    template<int blocksize, int smemSizeBytes>
    __global__
    void flagPairedCandidatesKernel(
        int numPairedAnchors,
        const int* __restrict__ d_numCandidatesPerAnchor,
        const int* __restrict__ d_numCandidatesPerAnchorPrefixSum,
        const read_number* __restrict__ d_candidateReadIds,
        bool* __restrict__ d_isPairedCandidate
    ){

        constexpr int numSharedElements = SDIV(smemSizeBytes, sizeof(int));

        __shared__ read_number sharedElements[numSharedElements];

        //search elements of array1 in array2. if found, set output element to true
        //array1 and array2 must be sorted
        auto process = [&](
            const read_number* array1,
            int numElements1,
            const read_number* array2,
            int numElements2,
            bool* output
        ){
            const int numIterations = SDIV(numElements2, numSharedElements);

            for(int iteration = 0; iteration < numIterations; iteration++){

                const int begin = iteration * numSharedElements;
                const int end = min((iteration+1) * numSharedElements, numElements2);
                const int num = end - begin;

                for(int i = threadIdx.x; i < num; i += blocksize){
                    sharedElements[i] = array2[begin + i];
                }

                __syncthreads();

                //TODO in iteration > 0, we may skip elements at the beginning of first range

                for(int i = threadIdx.x; i < numElements1; i += blocksize){
                    if(!output[i]){
                        const read_number readId = array1[i];
                        const read_number readIdToFind = readId % 2 == 0 ? readId + 1 : readId - 1;

                        const bool found = thrust::binary_search(thrust::seq, sharedElements, sharedElements + num, readIdToFind);
                        if(found){
                            output[i] = true;
                        }
                    }
                }

                __syncthreads();
            }
        };

        for(int a = blockIdx.x; a < numPairedAnchors; a += gridDim.x){
            const int firstTask = 2*a;
            //const int secondTask = firstTask + 1;

            //check for pairs in current candidates
            const int rangeBegin = d_numCandidatesPerAnchorPrefixSum[firstTask];                        
            const int rangeMid = d_numCandidatesPerAnchorPrefixSum[firstTask + 1];
            const int rangeEnd = rangeMid + d_numCandidatesPerAnchor[firstTask + 1];

            process(
                d_candidateReadIds + rangeBegin,
                rangeMid - rangeBegin,
                d_candidateReadIds + rangeMid,
                rangeEnd - rangeMid,
                d_isPairedCandidate + rangeBegin
            );

            process(
                d_candidateReadIds + rangeMid,
                rangeEnd - rangeMid,
                d_candidateReadIds + rangeBegin,
                rangeMid - rangeBegin,
                d_isPairedCandidate + rangeMid
            );
        }
    }


} //namespace gpucorrectorkernels   

} //namespace gpu
} //namespace care


#endif