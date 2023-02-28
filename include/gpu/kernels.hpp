#ifndef CARE_GPU_KERNELS_HPP
#define CARE_GPU_KERNELS_HPP

#include <hpc_helpers.cuh>
#include <gpu/gpumsa.cuh>
#include <gpu/cuda_block_select.cuh>
#include <cub/cub.cuh>

#include <alignmentorientation.hpp>
#include <correctedsequence.hpp>
#include <gpu/forest_gpu.cuh>

#include <config.hpp>
#include <gpu/sequenceconversionkernels.cuh>

#include <map>

#include <rmm/mr/device/per_device_resource.hpp>

namespace care {
namespace gpu {


#ifdef __NVCC__



struct AnchorHighQualityFlag{
    char data;

    __host__ __device__
    bool hq() const{
        return data == 1;
    }

    __host__ __device__
    void hq(bool isHq){
        data = isHq ? 1 : 0;
    }
};

struct IsHqAnchor{
    DEVICEQUALIFIER
    bool operator() (const AnchorHighQualityFlag& flag) const{
        return flag.hq();
    }
};



void callShiftedHammingDistanceKernel(
    int* d_alignment_overlaps,
    int* d_alignment_shifts,
    int* d_alignment_nOps,
    AlignmentOrientation* d_alignment_best_alignment_flags,
    const unsigned int* d_anchorSequencesData,
    const unsigned int* d_candidateSequencesData,
    const int* d_anchorSequencesLength,
    const int* d_candidateSequencesLength,
    const int* d_anchorIndicesOfCandidates,
    int numAnchors,
    int numCandidates,
    const bool* d_anchorContainsN,
    bool removeAmbiguousAnchors,
    const bool* d_candidateContainsN,
    bool removeAmbiguousCandidates,
    int maximumSequenceLengthAnchor,
    int maximumSequenceLengthCandidate,
    std::size_t encodedSequencePitchInInts2BitAnchor,
    std::size_t encodedSequencePitchInInts2BitCandidate,
    int min_overlap,
    float maxErrorRate,
    float min_overlap_ratio,
    float estimatedNucleotideErrorRate,
    cudaStream_t stream
);

void callRightShiftedHammingDistanceKernel(
    int* d_alignment_overlaps,
    int* d_alignment_shifts,
    int* d_alignment_nOps,
    AlignmentOrientation* d_alignment_best_alignment_flags,
    const unsigned int* d_anchorSequencesData,
    const unsigned int* d_candidateSequencesData,
    const int* d_anchorSequencesLength,
    const int* d_candidateSequencesLength,
    const int* d_anchorIndicesOfCandidates,
    int numAnchors,
    int numCandidates,
    const bool* d_anchorContainsN,
    bool removeAmbiguousAnchors,
    const bool* d_candidateContainsN,
    bool removeAmbiguousCandidates,
    int maximumSequenceLengthAnchor,
    int maximumSequenceLengthCandidate,
    std::size_t encodedSequencePitchInInts2BitAnchor,
    std::size_t encodedSequencePitchInInts2BitCandidate,
    int min_overlap,
    float maxErrorRate,
    float min_overlap_ratio,
    float estimatedNucleotideErrorRate,
    cudaStream_t stream
);

void callSelectIndicesOfGoodCandidatesKernelAsync(
            int* d_indicesOfGoodCandidates,
            int* d_numIndicesPerAnchor,
            int* d_totalNumIndices,
            const AlignmentOrientation* d_alignmentFlags,
            const int* d_candidates_per_anchor,
            const int* d_candidates_per_anchor_prefixsum,
            const int* d_anchorIndicesOfCandidates,
            int numAnchors,
            int numCandidates,
            cudaStream_t stream);


void call_cuda_filter_alignments_by_mismatchratio_kernel_async(
            AlignmentOrientation* d_bestAlignmentFlags,
            const int* d_nOps,
            const int* d_overlaps,
            const int* d_candidates_per_anchor_prefixsum,
            int numAnchors,
            int numCandidates,
            float mismatchratioBaseFactor,
            float goodAlignmentsCountThreshold,
            cudaStream_t stream);

// msa construction kernels

void callComputeMaximumMsaWidthKernel(
    int* d_result,
    const int* d_shifts,
    const int* d_anchorLengths,
    const int* d_candidateLengths,
    const int* d_indices,
    const int* d_indices_per_anchor,
    const int* d_candidatesPerAnchorPrefixSum,
    const int numAnchors,
    cudaStream_t stream
);

void callComputeMsaConsensusQualityKernel(
    char* d_consensusQuality,
    int consensusQualityPitchInBytes,
    GPUMultiMSA multiMSA,
    cudaStream_t stream
);

void callComputeDecodedMsaConsensusKernel(
    char* d_consensus,
    int consensusPitchInBytes,
    GPUMultiMSA multiMSA,
    cudaStream_t stream
);

void callComputeMsaSizesKernel(
    int* d_sizes,
    GPUMultiMSA multiMSA,
    cudaStream_t stream
);

void callConstructAndRefineMultipleSequenceAlignmentsKernel(
    int* __restrict__ d_newIndices,
    int* __restrict__ d_newNumIndicesPerAnchor,
    int* __restrict__ d_newNumIndices,
    GPUMultiMSA multiMSA,
    const int* __restrict__ overlaps,
    const int* __restrict__ shifts,
    const int* __restrict__ nOps,
    const AlignmentOrientation* __restrict__ bestAlignmentFlags,
    const int* __restrict__ anchorLengths,
    const int* __restrict__ candidateLengths,
    int* __restrict__ indices,
    int* __restrict__ indices_per_anchor,
    const int* __restrict__ candidatesPerAnchorPrefixSum,            
    const unsigned int* __restrict__ anchorSequencesData,
    const unsigned int* __restrict__ candidateSequencesData,
    const bool* __restrict__ d_isPairedCandidate,
    const char* __restrict__ anchorQualities,
    const char* __restrict__ candidateQualities,
    bool* __restrict__ d_shouldBeKept,
    int numAnchors,
    float desiredAlignmentMaxErrorRate,
    bool canUseQualityScores,
    int encodedSequencePitchInInts,
    size_t qualityPitchInBytes,
    int dataset_coverage,
    int numRefinementIterations,
    cudaStream_t stream
);

void callConstructMultipleSequenceAlignmentsKernel_async(
    GPUMultiMSA multiMSA,
    const int* d_overlaps,
    const int* d_shifts,
    const int* d_nOps,
    const AlignmentOrientation* d_bestAlignmentFlags,
    const int* d_anchorLengths,
    const int* d_candidateLengths,
    const int* d_indices,
    const int* d_indices_per_anchor,
    const int* d_candidatesPerAnchorPrefixSum,            
    const unsigned int* d_anchorSequencesData,
    const unsigned int* d_candidateSequencesData,
    const bool* d_isPairedCandidate,
    const char* d_anchorQualities,
    const char* d_candidateQualities,
    float desiredAlignmentMaxErrorRate,
    int numAnchors,
    bool canUseQualityScores,
    int encodedSequencePitchInInts,
    size_t qualityPitchInBytes,
    cudaStream_t stream
);

void callMsaCandidateRefinementKernel_multiiter_async(
    int* d_newIndices,
    int* d_newNumIndicesPerAnchor,
    int* d_newNumIndices,
    GPUMultiMSA multiMSA,
    const AlignmentOrientation* d_bestAlignmentFlags,
    const int* d_shifts,
    const int* d_nOps,
    const int* d_overlaps,
    const unsigned int* d_anchorSequencesData,
    const unsigned int* d_candidateSequencesData,
    const bool* d_isPairedCandidate,
    const int* d_anchorSequencesLength,
    const int* d_candidateSequencesLength,
    const char* d_anchorQualities,
    const char* d_candidateQualities,
    bool* d_shouldBeKept,
    const int* d_candidates_per_anchor_prefixsum,
    float desiredAlignmentMaxErrorRate,
    int numAnchors,
    bool canUseQualityScores,
    size_t encodedSequencePitchInInts,
    size_t qualityPitchInBytes,
    int* d_indices,
    int* d_indices_per_anchor,
    int dataset_coverage,
    int numIterations,
    cudaStream_t stream
);



// correction kernels


void call_msaCorrectAnchorsKernel_async(
    char* d_correctedAnchors,
    bool* d_anchorIsCorrected,
    AnchorHighQualityFlag* d_isHighQualityAnchor,
    GPUMultiMSA multiMSA,
    const unsigned int* d_anchorSequencesData,
    const int* d_indices_per_anchor,
    const int* d_numAnchors,
    int maxNumAnchors,
    int encodedSequencePitchInInts,
    size_t sequence_pitch,
    float estimatedErrorrate,
    float avg_support_threshold,
    float min_support_threshold,
    float min_coverage_threshold,
    int maximum_sequence_length,
    cudaStream_t stream
);


void callFlagCandidatesToBeCorrectedKernel_async(
    bool* d_candidateCanBeCorrected,
    int* d_numCorrectedCandidatesPerAnchor,
    GPUMultiMSA multiMSA,
    const int* d_alignmentShifts,
    const int* d_candidateSequencesLengths,
    const int* d_anchorIndicesOfCandidates,
    const AnchorHighQualityFlag* d_hqflags,
    const int* d_candidatesPerAnchorPrefixsum,
    const int* d_localGoodCandidateIndices,
    const int* d_numLocalGoodCandidateIndicesPerAnchor,
    const int* d_numAnchors,
    const int* d_numCandidates,
    float min_support_threshold,
    float min_coverage_threshold,
    int new_columns_to_correct,
    cudaStream_t stream
);

void callFlagCandidatesToBeCorrectedWithExcludeFlagsKernel(
    bool* d_candidateCanBeCorrected,
    int* d_numCorrectedCandidatesPerAnchor,
    GPUMultiMSA multiMSA,
    const bool* d_excludeFlags, //candidates with flag == true will not be considered
    const int* d_alignmentShifts,
    const int* d_candidateSequencesLengths,
    const int* d_anchorIndicesOfCandidates,
    const AnchorHighQualityFlag* d_hqflags,
    const int* d_candidatesPerAnchorPrefixsum,
    const int* d_localGoodCandidateIndices,
    const int* d_numLocalGoodCandidateIndicesPerAnchor,
    const int* d_numAnchors,
    const int* d_numCandidates,
    float min_support_threshold,
    float min_coverage_threshold,
    int new_columns_to_correct,
    cudaStream_t stream
);

void callCorrectCandidatesKernel(
    char* __restrict__ correctedCandidates,
    GPUMultiMSA multiMSA,
    const int* __restrict__ shifts,
    const AlignmentOrientation* __restrict__ bestAlignmentFlags,
    const unsigned int* __restrict__ candidateSequencesData,
    const int* __restrict__ candidateSequencesLengths,
    const bool* __restrict__ d_candidateContainsN,
    const int* __restrict__ candidateIndicesOfCandidatesToBeCorrected,
    const int* __restrict__ numCandidatesToBeCorrected,
    const int* __restrict__ anchorIndicesOfCandidates,
    int numCandidates,
    int encodedSequencePitchInInts,
    size_t decodedSequencePitchInBytes,
    int maximum_sequence_length,
    cudaStream_t stream
);

void callCorrectCandidatesAndComputeEditsKernel(
    char* __restrict__ correctedCandidates,
    EncodedCorrectionEdit* __restrict__ d_editsPerCorrectedCandidate,
    int* __restrict__ d_numEditsPerCorrectedCandidate,
    GPUMultiMSA multiMSA,
    const int* __restrict__ shifts,
    const AlignmentOrientation* __restrict__ bestAlignmentFlags,
    const unsigned int* __restrict__ candidateSequencesData,
    const int* __restrict__ candidateSequencesLengths,
    const bool* __restrict__ d_candidateContainsN,
    const int* __restrict__ candidateIndicesOfCandidatesToBeCorrected,
    const int* __restrict__ numCandidatesToBeCorrected,
    const int* __restrict__ anchorIndicesOfCandidates,
    const int* d_numAnchors,
    const int* d_numCandidates,
    int doNotUseEditsValue,
    int numEditsThreshold,
    int encodedSequencePitchInInts,
    size_t decodedSequencePitchInBytes,
    size_t editsPitchInBytes,
    int maximum_sequence_length,
    cudaStream_t stream
);

void callMsaCorrectAnchorsWithForestKernel(
    char* d_correctedAnchors,
    bool* d_anchorIsCorrected,
    AnchorHighQualityFlag* d_isHighQualityAnchor,
    GPUMultiMSA multiMSA,
    GpuForest::Clf gpuForest,
    float forestThreshold,
    const unsigned int* d_anchorSequencesData,
    const int* d_indices_per_anchor,
    const int numAnchors,
    int encodedSequencePitchInInts,
    size_t decodedSequencePitchInBytes,
    int maximumSequenceLength,
    float estimatedErrorrate,
    float estimatedCoverage,
    float avg_support_threshold,
    float min_support_threshold,
    float min_coverage_threshold,
    cudaStream_t stream
);

void callMsaCorrectCandidatesWithForestKernel(
    char* d_correctedCandidates,
    GPUMultiMSA multiMSA,
    GpuForest::Clf gpuForest,
    float forestThreshold,
    float estimatedCoverage,
    const int* d_shifts,
    const AlignmentOrientation* d_bestAlignmentFlags,
    const unsigned int* d_candidateSequencesData,
    const int* d_candidateSequencesLengths,
    const int* d_candidateIndicesOfCandidatesToBeCorrected,
    const int* d_anchorIndicesOfCandidates,
    const int numCandidates,        
    int encodedSequencePitchInInts,
    size_t decodedSequencePitchInBytes,
    int maximum_sequence_length,
    cudaStream_t stream
);

void callConstructSequenceCorrectionResultsKernel(
    EncodedCorrectionEdit* d_edits,
    int* d_numEditsPerCorrection,
    int doNotUseEditsValue,
    const int* d_indicesOfUncorrectedSequences,
    const int* d_numIndices,
    const bool* d_readContainsN,
    const unsigned int* d_uncorrectedEncodedSequences,
    const int* d_sequenceLengths,
    const char* d_correctedSequences,
    const int numCorrectedSequencesUpperBound, // >= *d_numIndices. d_edits must be large enought to store the edits of this many sequences
    bool isCompactCorrection,
    int numEditsThreshold,
    size_t encodedSequencePitchInInts,
    size_t decodedSequencePitchInBytes,
    size_t editsPitchInBytes,        
    cudaStream_t stream
);



/*
    If toFind[s] exists in segment s, remove it from this segment s 
    by shifting remaining elements to the left.
    segmentSizes[s] will be updated.
    */
    template<class T, int BLOCKSIZE, int ITEMS_PER_THREAD>
    __global__ 
    void findAndRemoveFromSegmentKernel(
        const T* __restrict__ toFind,
        T* items,
        int numSegments,
        int* __restrict__ segmentSizes,
        const int* __restrict__ segmentBeginOffsets
    ){
        constexpr int itemsPerIteration = ITEMS_PER_THREAD * BLOCKSIZE;

        assert(BLOCKSIZE == blockDim.x);

        using BlockLoad = cub::BlockLoad<T, BLOCKSIZE, ITEMS_PER_THREAD, cub::BLOCK_LOAD_WARP_TRANSPOSE>;
        using MyBlockSelect = BlockSelect<T, BLOCKSIZE>;

        __shared__ union{
            typename BlockLoad::TempStorage load;
            typename MyBlockSelect::TempStorage select;
        } temp_storage;

        for(int s = blockIdx.x; s < numSegments; s += gridDim.x){
            const int segmentsize = segmentSizes[s];
            const int beginOffset = segmentBeginOffsets[s];
            const T idToRemove = toFind[s];

            const int numIterations = SDIV(segmentsize, itemsPerIteration);
            T myitems[ITEMS_PER_THREAD];
            int flags[ITEMS_PER_THREAD];

            int numSelectedTotal = 0;
            int remainingItems = segmentsize;
            const T* inputdata = items + beginOffset;
            T* outputdata = items + beginOffset;

            for(int iter = 0; iter < numIterations; iter++){
                const int validItems = min(remainingItems, itemsPerIteration);
                BlockLoad(temp_storage.load).Load(inputdata, myitems, validItems);

                #pragma unroll
                for(int i = 0; i < ITEMS_PER_THREAD; i++){
                    if(threadIdx.x * ITEMS_PER_THREAD + i < validItems && myitems[i] != idToRemove){
                        flags[i] = 1;
                    }else{
                        flags[i] = 0;
                    }
                }

                __syncthreads();

                const int numSelected = MyBlockSelect(temp_storage.select).ForEachFlagged(myitems, flags, validItems,
                    [&](const auto& item, const int& pos){
                        outputdata[pos] = item;
                    }
                );
                assert(numSelected <= validItems);

                numSelectedTotal += numSelected;
                outputdata += numSelected;
                inputdata += validItems;
                remainingItems -= validItems;

                __syncthreads();
            }

            assert(segmentsize >= numSelectedTotal);

            //update segment size
            if(numSelectedTotal != segmentsize){
                if(threadIdx.x == 0){
                    segmentSizes[s] = numSelectedTotal;
                }
            }
        }
    }

    template<class T, int BLOCKSIZE, int ITEMS_PER_THREAD>
    void callFindAndRemoveFromSegmentKernel(
        const T* d_toFind,
        T* d_items,
        int numSegments,
        int* d_segmentSizes,
        const int* d_segmentBeginOffsets,
        cudaStream_t stream
    ){
        if(numSegments <= 0){
            return;
        }

        dim3 block = BLOCKSIZE;
        dim3 grid = numSegments;

        findAndRemoveFromSegmentKernel<T, BLOCKSIZE, ITEMS_PER_THREAD><<<grid, block, 0, stream>>>
            (d_toFind, d_items, numSegments, d_segmentSizes, d_segmentBeginOffsets);
    }



#endif //ifdef __NVCC__

}
}


#endif
