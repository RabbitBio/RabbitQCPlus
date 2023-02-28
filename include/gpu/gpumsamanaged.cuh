#ifndef CARE_GPU_MSA_MANAGED_CUH
#define CARE_GPU_MSA_MANAGED_CUH

#include <gpu/gpumsa.cuh>
#include <gpu/cubvector.cuh>
#include <gpu/kernels.hpp>
#include <hpc_helpers.cuh>
#include <gpu/cudaerrorcheck.cuh>
#include <memorymanagement.hpp>
#include <cub/cub.cuh>
#include <rmm/mr/device/device_memory_resource.hpp>
#include <rmm/device_uvector.hpp>
#include <rmm/device_scalar.hpp>
#include <gpu/rmm_utilities.cuh>


#include <cstdint>
#include <limits>

namespace care{
namespace gpu{


    struct MSAColumnCount{
    public:
        constexpr MSAColumnCount(int c) noexcept : value(c){}

        constexpr int count() const noexcept{
            return value;
        }

        static constexpr int unknown() noexcept{
            return std::numeric_limits<int>::max();
        }
    private:
        int value;
    };

    bool operator==(MSAColumnCount l, MSAColumnCount r) noexcept{
        return l.count() == r.count();
    }

    bool operator!=(MSAColumnCount l, MSAColumnCount r) noexcept{
        return !(operator==(l,r));
    }


    struct ManagedGPUMultiMSA{
    public:
        ManagedGPUMultiMSA(cudaStream_t stream, rmm::mr::device_memory_resource* mr_, int* h_tempstorage = nullptr) 
            : mr(mr_),
            d_consensusEncoded(0, stream, mr),
            d_counts(0, stream, mr),
            d_coverages(0, stream, mr),
            d_origCoverages(0, stream, mr),
            d_weights(0, stream, mr),
            d_support(0, stream, mr),
            d_origWeights(0, stream, mr),
            d_columnProperties(0, stream, mr){

            if(h_tempstorage != nullptr){
                tempvalue = h_tempstorage;
            }else{
                pinnedValue.resize(1);
                tempvalue = pinnedValue.data();
            }
        }

        void constructAndRefine(
            int* d_newIndices,
            int* d_newNumIndicesPerAnchor,
            int* d_newNumIndices,
            const int* d_overlaps,
            const int* d_shifts,
            const int* d_nOps,
            const AlignmentOrientation* d_bestAlignmentFlags,
            const int* d_anchorLengths,
            const int* d_candidateLengths,
            int* d_candidatePositionsInSegments,
            int* d_numCandidatePositionsInSegments,
            const int* d_segmentBeginOffsets,            
            const unsigned int* d_anchorSequencesData,
            const unsigned int* d_candidateSequencesData,
            const bool* d_isPairedCandidate,
            const char* d_anchorQualities,
            const char* d_candidateQualities,
            int numAnchors,
            int maxNumCandidates,
            float desiredAlignmentMaxErrorRate,
            bool canUseQualityScores,
            int encodedSequencePitchInInts,
            size_t qualityPitchInBytes,
            int dataset_coverage,
            int numRefinementIterations,
            MSAColumnCount maximumMsaWidth, // upper bound for number of columns in a single msa. must be large enough to actually fit the data.
            cudaStream_t stream
        ){
            initializeBuffers(
                maximumMsaWidth, 
                numAnchors, 
                d_shifts,
                d_anchorLengths,
                d_candidateLengths,
                d_candidatePositionsInSegments,
                d_numCandidatePositionsInSegments,
                d_segmentBeginOffsets,
                stream
            );

            rmm::device_uvector<bool> d_temp(maxNumCandidates, stream, mr);

            callConstructAndRefineMultipleSequenceAlignmentsKernel(
                d_newIndices,
                d_newNumIndicesPerAnchor,
                d_newNumIndices,
                multiMSA,
                d_overlaps,
                d_shifts,
                d_nOps,
                d_bestAlignmentFlags,
                d_anchorLengths,
                d_candidateLengths,
                d_candidatePositionsInSegments,
                d_numCandidatePositionsInSegments,
                d_segmentBeginOffsets,            
                d_anchorSequencesData,
                d_candidateSequencesData,
                d_isPairedCandidate,
                d_anchorQualities,
                d_candidateQualities,
                d_temp.data(),
                numAnchors,
                desiredAlignmentMaxErrorRate,
                canUseQualityScores,
                encodedSequencePitchInInts,
                qualityPitchInBytes,
                dataset_coverage,
                numRefinementIterations,
                stream
            );
        }

        void construct(
            const int* d_alignment_overlaps,
            const int* d_alignment_shifts,
            const int* d_alignment_nOps,
            const AlignmentOrientation* d_alignment_best_alignment_flags,
            const int* d_candidatePositionsInSegments,
            const int* d_numCandidatePositionsInSegments,
            const int* d_segmentBeginOffsets,
            const int* d_anchorSequencesLength,
            const unsigned int* d_anchorSequences,
            const char* d_anchorQualities,
            int numAnchors,
            const int* d_candidateSequencesLength,
            const unsigned int* d_candidateSequences,
            const char* d_candidateQualities,
            const bool* d_isPairedCandidate,
            std::size_t encodedSequencePitchInInts,
            std::size_t qualityPitchInBytes,
            bool useQualityScores,
            float desiredAlignmentMaxErrorRate,
            MSAColumnCount maximumMsaWidth, // upper bound for number of columns in a single msa. must be large enough to actually fit the data.
            cudaStream_t stream
        ){
            initializeBuffers(
                maximumMsaWidth, 
                numAnchors, 
                d_alignment_shifts,
                d_anchorSequencesLength,
                d_candidateSequencesLength,
                d_candidatePositionsInSegments,
                d_numCandidatePositionsInSegments,
                d_segmentBeginOffsets,
                stream
            );

            callConstructMultipleSequenceAlignmentsKernel_async(
                multiMSA,
                d_alignment_overlaps,
                d_alignment_shifts,
                d_alignment_nOps,
                d_alignment_best_alignment_flags,
                d_anchorSequencesLength,
                d_candidateSequencesLength,
                d_candidatePositionsInSegments,
                d_numCandidatePositionsInSegments,
                d_segmentBeginOffsets,
                d_anchorSequences,
                d_candidateSequences,
                d_isPairedCandidate,
                d_anchorQualities,
                d_candidateQualities,
                desiredAlignmentMaxErrorRate,
                numAnchors,
                useQualityScores,
                encodedSequencePitchInInts,
                qualityPitchInBytes,
                stream
            );
        }

        void refine(
            int* d_newCandidatePositionsInSegments,
            int* d_newNumCandidatePositionsInSegments,
            int* d_newNumCandidates,
            const int* d_alignment_overlaps,
            const int* d_alignment_shifts,
            const int* d_alignment_nOps,
            const AlignmentOrientation* d_alignment_best_alignment_flags,
            int* d_candidatePositionsInSegments,
            int* d_numCandidatePositionsInSegments,
            const int* d_segmentBeginOffsets,
            const int* d_anchorSequencesLength,
            const unsigned int* d_anchorSequences,
            const char* d_anchorQualities,
            int numAnchors,
            const int* d_candidateSequencesLength,
            const unsigned int* d_candidateSequences,
            const char* d_candidateQualities,
            const bool* d_isPairedCandidate,
            int maxNumCandidates,
            std::size_t encodedSequencePitchInInts,
            std::size_t qualityPitchInBytes,
            bool useQualityScores,
            float desiredAlignmentMaxErrorRate,
            int dataset_coverage,
            int numIterations,
            cudaStream_t stream
        ){
            //std::cerr << "thread " << std::this_thread::get_id() << " msa refine, stream " << stream << "\n";

            rmm::device_uvector<bool> d_temp(maxNumCandidates, stream, mr);

            callMsaCandidateRefinementKernel_multiiter_async(
                d_newCandidatePositionsInSegments,
                d_newNumCandidatePositionsInSegments,
                d_newNumCandidates,
                multiMSA,
                d_alignment_best_alignment_flags,
                d_alignment_shifts,
                d_alignment_nOps,
                d_alignment_overlaps,
                d_anchorSequences,
                d_candidateSequences,
                d_isPairedCandidate,
                d_anchorSequencesLength,
                d_candidateSequencesLength,
                d_anchorQualities,
                d_candidateQualities,
                d_temp.data(),
                d_segmentBeginOffsets,
                desiredAlignmentMaxErrorRate,
                numAnchors,
                useQualityScores,
                encodedSequencePitchInInts,
                qualityPitchInBytes,
                d_candidatePositionsInSegments,
                d_numCandidatePositionsInSegments,
                dataset_coverage,
                numIterations,
                stream
            );
        }

        void computeConsensusQuality(
            char* d_consensusQuality,
            int consensusQualityPitchInBytes,
            cudaStream_t stream
        ){
            callComputeMsaConsensusQualityKernel(
                d_consensusQuality,
                consensusQualityPitchInBytes,
                multiMSA,
                stream
            );
        }

        void computeConsensus(
            char* d_consensus,
            int consensusPitchInBytes,
            cudaStream_t stream
        ){
            callComputeDecodedMsaConsensusKernel(
                d_consensus,
                consensusPitchInBytes,
                multiMSA,
                stream
            );
        }

        void computeMsaSizes(
            int* d_sizes,
            cudaStream_t stream
        ){
            callComputeMsaSizesKernel(
                d_sizes,
                multiMSA,
                stream
            );
        }

        MemoryUsage getMemoryInfo() const{
            int deviceId;
            CUDACHECK(cudaGetDevice(&deviceId));

            MemoryUsage info{};
            // auto handleHost = [&](const auto& h){
            //     info.host += h.sizeInBytes();
            // };
            auto handleDevice = [&](const auto& d){
                using ElementType = typename std::remove_reference<decltype(d)>::type::value_type;
                info.device[deviceId] += d.size() * sizeof(ElementType);
            };

            handleDevice(d_consensusEncoded);
            handleDevice(d_counts);
            handleDevice(d_coverages);
            handleDevice(d_origCoverages);
            handleDevice(d_weights);
            handleDevice(d_support);
            handleDevice(d_origWeights);
            handleDevice(d_columnProperties);

            return info;
        }

        void destroy(cudaStream_t stream){
            ::destroy(d_consensusEncoded, stream);
            ::destroy(d_counts, stream);
            ::destroy(d_coverages, stream);
            ::destroy(d_origCoverages, stream);
            ::destroy(d_weights, stream);
            ::destroy(d_support, stream);
            ::destroy(d_origWeights, stream);
            ::destroy(d_columnProperties, stream);

            multiMSA.numMSAs = 0;
            multiMSA.columnPitchInElements = 0;
            multiMSA.counts = nullptr;
            multiMSA.weights = nullptr;
            multiMSA.coverages = nullptr;
            multiMSA.consensus = nullptr;
            multiMSA.support = nullptr;
            multiMSA.origWeights = nullptr;
            multiMSA.origCoverages = nullptr;
            multiMSA.columnProperties = nullptr;
        }

        GPUMultiMSA multiMSAView() const{
            return multiMSA;
        }

        int getMaximumMsaWidth() const{
            return columnPitchInElements;
        }

    private:
        void initializeBuffers(
            MSAColumnCount maximumMsaWidth, 
            int numAnchors, 
            const int* d_alignment_shifts,
            const int* d_anchorSequencesLength,
            const int* d_candidateSequencesLength,
            const int* d_candidatePositionsInSegments,
            const int* d_numCandidatePositionsInSegments,
            const int* d_segmentBeginOffsets,
            cudaStream_t stream
        ){
            if(maximumMsaWidth == MSAColumnCount::unknown()){
                rmm::device_scalar<int> d_maxMsaWidth(stream, mr);

                gpu::callComputeMaximumMsaWidthKernel(
                    d_maxMsaWidth.data(),
                    d_alignment_shifts,
                    d_anchorSequencesLength,
                    d_candidateSequencesLength,
                    d_candidatePositionsInSegments,
                    d_numCandidatePositionsInSegments,
                    d_segmentBeginOffsets,
                    numAnchors,
                    stream
                );

                CUDACHECK(cudaMemcpyAsync(
                    pinnedValue.data(),
                    d_maxMsaWidth.data(),
                    sizeof(int),
                    D2H,
                    stream
                ));

                CUDACHECK(cudaStreamSynchronize(stream));

                columnPitchInElements = *pinnedValue;
            }else{
                columnPitchInElements = maximumMsaWidth.count();
            }

            numMSAs = numAnchors;

            //pad to 128
            columnPitchInElements = SDIV(columnPitchInElements, 128) * 128;

            resizeUninitialized(d_consensusEncoded, columnPitchInElements * numAnchors, stream);
            resizeUninitialized(d_counts, 4 * columnPitchInElements * numAnchors, stream);
            resizeUninitialized(d_coverages, columnPitchInElements * numAnchors, stream);
            resizeUninitialized(d_origCoverages, columnPitchInElements * numAnchors, stream);
            resizeUninitialized(d_weights, 4 * columnPitchInElements * numAnchors, stream);
            resizeUninitialized(d_support, columnPitchInElements * numAnchors, stream);
            resizeUninitialized(d_origWeights, columnPitchInElements * numAnchors, stream);
            resizeUninitialized(d_columnProperties, numAnchors, stream);
            // CUDACHECK(cudaMemsetAsync(d_consensusEncoded.data(), 0, sizeof(std::uint8_t) * columnPitchInElements * numAnchors, stream));
            // CUDACHECK(cudaMemsetAsync(d_counts.data(), 0, sizeof(int) * 4*columnPitchInElements * numAnchors, stream));
            // CUDACHECK(cudaMemsetAsync(d_coverages.data(), 0, sizeof(int) * columnPitchInElements * numAnchors, stream));
            // CUDACHECK(cudaMemsetAsync(d_origCoverages.data(), 0, sizeof(int) * columnPitchInElements * numAnchors, stream));
            // CUDACHECK(cudaMemsetAsync(d_weights.data(), 0, sizeof(float) * 4*columnPitchInElements * numAnchors, stream));
            // CUDACHECK(cudaMemsetAsync(d_support.data(), 0, sizeof(float) * columnPitchInElements * numAnchors, stream));
            // CUDACHECK(cudaMemsetAsync(d_origWeights.data(), 0, sizeof(float) * columnPitchInElements * numAnchors, stream));
            // CUDACHECK(cudaMemsetAsync(d_columnProperties.data(), 0, sizeof(MSAColumnProperties) * numAnchors, stream));

            multiMSA.numMSAs = numMSAs;
            multiMSA.columnPitchInElements = columnPitchInElements;
            multiMSA.counts = d_counts.data();
            multiMSA.weights = d_weights.data();
            multiMSA.coverages = d_coverages.data();
            multiMSA.consensus = d_consensusEncoded.data();
            multiMSA.support = d_support.data();
            multiMSA.origWeights = d_origWeights.data();
            multiMSA.origCoverages = d_origCoverages.data();
            multiMSA.columnProperties = d_columnProperties.data();
        }
        public:  

        int numMSAs{};
        int columnPitchInElements{};
        int* tempvalue{};

        helpers::SimpleAllocationPinnedHost<int, 0> pinnedValue{};
        CudaEvent event{cudaEventDisableTiming};

        rmm::mr::device_memory_resource* mr;

        rmm::device_uvector<std::uint8_t> d_consensusEncoded;
        rmm::device_uvector<int> d_counts;
        rmm::device_uvector<int> d_coverages;
        rmm::device_uvector<int> d_origCoverages;
        rmm::device_uvector<float> d_weights;
        rmm::device_uvector<float> d_support;
        rmm::device_uvector<float> d_origWeights;
        rmm::device_uvector<MSAColumnProperties> d_columnProperties;

        GPUMultiMSA multiMSA{};
    };


} //namespace gpu
} //namespace care

#endif