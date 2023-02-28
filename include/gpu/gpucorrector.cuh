#ifndef CARE_GPUCORRECTOR_CUH
#define CARE_GPUCORRECTOR_CUH


#include <hpc_helpers.cuh>
#include <hpc_helpers/include/nvtx_markers.cuh>

#include <gpu/gpuminhasher.cuh>
#include <gpu/kernels.hpp>
#include <gpu/gpucorrectorkernels.cuh>
#include <gpu/gpureadstorage.cuh>
#include <gpu/asyncresult.cuh>
#include <gpu/cudaerrorcheck.cuh>
#include <gpu/cubwrappers.cuh>
#include <gpu/gpumsamanaged.cuh>
#include <gpu/global_cuda_stream_pool.cuh>
#include <gpu/minhashqueryfilter.cuh>

#include <config.hpp>
#include <util.hpp>
#include <corrector_common.hpp>
#include <threadpool.hpp>

#include <options.hpp>
#include <correctedsequence.hpp>
#include <memorymanagement.hpp>
#include <msa.hpp>
#include <classification.hpp>

#include <forest.hpp>
#include <gpu/forest_gpu.cuh>

#include <algorithm>
#include <array>
#include <map>

#include <cub/cub.cuh>

#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/iterator/transform_output_iterator.h>
#include <thrust/gather.h>
#include <thrust/scan.h>
#include <thrust/unique.h>
#include <thrust/equal.h>

#include <rmm/mr/device/per_device_resource.hpp>
#include <rmm/mr/device/cuda_async_memory_resource.hpp>
#include <rmm/exec_policy.hpp>
#include <gpu/rmm_utilities.cuh>


namespace care{
namespace gpu{

    struct SequenceFlagMultiplier{
        int pitch{};
        const bool* flags{};

        __host__ __device__
        SequenceFlagMultiplier(const bool* flags_, int pitch_)
            :pitch(pitch_), flags(flags_){

        }

        __host__ __device__
        bool operator()(int i) const{
            return flags[i / pitch];
        }
    };

    template<class Iter>
    struct IteratorMultiplier{
        using value_type = typename std::iterator_traits<Iter>::value_type;

        int factor{};
        Iter data{};

        __host__ __device__
        IteratorMultiplier(Iter data_, int factor_)
            : factor(factor_), data(data_){

        }

        __host__ __device__
        value_type operator()(int i) const{
            return *(data + (i / factor));
        }
    };

    template<class Iter>
    IteratorMultiplier<Iter> make_iterator_multiplier(Iter data, int factor){
        return IteratorMultiplier<Iter>{data, factor};
    }

    struct ReplaceNumberOp{
        int doNotUseEditsValue;
        int decodedSequencePitchInBytes;

        ReplaceNumberOp(int val, int pitch) 
            : doNotUseEditsValue(val), decodedSequencePitchInBytes(pitch)
        {}

        __host__ __device__
        int operator()(const int num) const noexcept{
            return num == doNotUseEditsValue ? decodedSequencePitchInBytes : 0;
        }
    };


    class GpuErrorCorrectorInput{
    public:
        template<class T>
        using PinnedBuffer = helpers::SimpleAllocationPinnedHost<T>;


        CudaEvent event{cudaEventDisableTiming};

        rmm::mr::device_memory_resource* mr;

        PinnedBuffer<int> h_numAnchors;
        PinnedBuffer<int> h_numCandidates;
        PinnedBuffer<read_number> h_anchorReadIds;
        PinnedBuffer<read_number> h_candidate_read_ids;

        rmm::device_uvector<int> d_numAnchors;
        rmm::device_uvector<int> d_numCandidates;
        rmm::device_uvector<read_number> d_anchorReadIds;
        rmm::device_uvector<unsigned int> d_anchor_sequences_data;
        rmm::device_uvector<int> d_anchor_sequences_lengths;
        rmm::device_uvector<read_number> d_candidate_read_ids;
        rmm::device_uvector<unsigned int> d_candidate_sequences_data;
        rmm::device_uvector<int> d_candidate_sequences_lengths;
        rmm::device_uvector<int> d_candidates_per_anchor;
        rmm::device_uvector<int> d_candidates_per_anchor_prefixsum;

        GpuErrorCorrectorInput()
        : mr(rmm::mr::get_current_device_resource()),
            d_numAnchors(0, cudaStreamPerThread, mr),
            d_numCandidates(0, cudaStreamPerThread, mr),
            d_anchorReadIds(0, cudaStreamPerThread, mr),
            d_anchor_sequences_data(0, cudaStreamPerThread, mr),
            d_anchor_sequences_lengths(0, cudaStreamPerThread, mr),
            d_candidate_read_ids(0, cudaStreamPerThread, mr),
            d_candidate_sequences_data(0, cudaStreamPerThread, mr),
            d_candidate_sequences_lengths(0, cudaStreamPerThread, mr),
            d_candidates_per_anchor(0, cudaStreamPerThread, mr),
            d_candidates_per_anchor_prefixsum(0, cudaStreamPerThread, mr)
        {
            CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));
        }

        MemoryUsage getMemoryInfo() const{
            MemoryUsage info{};
            auto handleHost = [&](const auto& h){
                info.host += h.sizeInBytes();
            };
            auto handleDevice = [&](const auto& d){
                using ElementType = typename std::remove_reference<decltype(d)>::type::value_type;
                info.device[event.getDeviceId()] += d.size() * sizeof(ElementType);
            };

            handleHost(h_numAnchors);
            handleHost(h_numCandidates);
            handleHost(h_anchorReadIds);
            handleHost(h_candidate_read_ids);

            handleDevice(d_numAnchors);
            handleDevice(d_numCandidates);
            handleDevice(d_anchorReadIds);
            handleDevice(d_anchor_sequences_data);
            handleDevice(d_anchor_sequences_lengths);
            handleDevice(d_candidate_read_ids);
            handleDevice(d_candidate_sequences_data);
            handleDevice(d_candidate_sequences_lengths);
            handleDevice(d_candidates_per_anchor);
            handleDevice(d_candidates_per_anchor_prefixsum);

            return info;
        }  
    };

    class GpuErrorCorrectorRawOutput{
    public:
        template<class T>
        using PinnedBuffer = helpers::SimpleAllocationPinnedHost<T>;

        bool nothingToDo;
        int numAnchors;
        int numCandidates;
        int doNotUseEditsValue;
        std::size_t editsPitchInBytes;
        std::size_t decodedSequencePitchInBytes;
        CudaEvent event{cudaEventDisableTiming};
        PinnedBuffer<read_number> h_anchorReadIds;
        PinnedBuffer<read_number> h_candidate_read_ids;
        PinnedBuffer<bool> h_anchor_is_corrected;
        PinnedBuffer<AnchorHighQualityFlag> h_is_high_quality_anchor;
        PinnedBuffer<int> h_num_corrected_candidates_per_anchor;
        PinnedBuffer<int> h_num_corrected_candidates_per_anchor_prefixsum;
        PinnedBuffer<int> h_indices_of_corrected_candidates;

        PinnedBuffer<int> h_candidate_sequences_lengths;
        PinnedBuffer<int> h_numEditsPerCorrectedanchor;
        PinnedBuffer<EncodedCorrectionEdit> h_editsPerCorrectedanchor;
        PinnedBuffer<char> h_corrected_anchors;
        PinnedBuffer<int> h_anchor_sequences_lengths;
        PinnedBuffer<char> h_corrected_candidates;
        PinnedBuffer<int> h_alignment_shifts;
        PinnedBuffer<int> h_numEditsPerCorrectedCandidate;
        PinnedBuffer<EncodedCorrectionEdit> h_editsPerCorrectedCandidate;
        PinnedBuffer<int> h_anchorEditOffsets;
        PinnedBuffer<int> h_correctedAnchorsOffsets;
        PinnedBuffer<int> h_candidateEditOffsets;
        PinnedBuffer<int> h_correctedCandidatesOffsets;

        MemoryUsage getMemoryInfo() const{
            MemoryUsage info{};
            auto handleHost = [&](const auto& h){
                info.host += h.sizeInBytes();
            };

            handleHost(h_anchorReadIds);
            handleHost(h_candidate_read_ids);
            handleHost(h_anchor_is_corrected);
            handleHost(h_is_high_quality_anchor);
            handleHost(h_num_corrected_candidates_per_anchor);
            handleHost(h_num_corrected_candidates_per_anchor_prefixsum);
            handleHost(h_indices_of_corrected_candidates);
            handleHost(h_candidate_sequences_lengths);
            handleHost(h_numEditsPerCorrectedanchor);
            handleHost(h_editsPerCorrectedanchor);
            handleHost(h_corrected_anchors);
            handleHost(h_anchor_sequences_lengths);
            handleHost(h_corrected_candidates);
            handleHost(h_alignment_shifts);
            handleHost(h_numEditsPerCorrectedCandidate);
            handleHost(h_editsPerCorrectedCandidate);
            handleHost(h_anchorEditOffsets);
            handleHost(h_correctedAnchorsOffsets);
            handleHost(h_candidateEditOffsets);
            handleHost(h_correctedCandidatesOffsets);

            return info;
        }  
    };



    class GpuAnchorHasher{
    public:

        GpuAnchorHasher() = default;

        GpuAnchorHasher(
            const GpuReadStorage& gpuReadStorage_,
            const GpuMinhasher& gpuMinhasher_,
            ThreadPool* threadPool_,
            rmm::mr::device_memory_resource* mr_
        ) : 
            gpuReadStorage{&gpuReadStorage_},
            gpuMinhasher{&gpuMinhasher_},
            threadPool{threadPool_},
            minhashHandle{gpuMinhasher->makeMinhasherHandle()},
            readstorageHandle{gpuReadStorage->makeHandle()},
            mr{mr_}
        {
            CUDACHECK(cudaGetDevice(&deviceId));            

            maxCandidatesPerRead = gpuMinhasher->getNumResultsPerMapThreshold() * gpuMinhasher->getNumberOfMaps();

            previousBatchFinishedEvent = CudaEvent{};

            encodedSequencePitchInInts = SequenceHelpers::getEncodedNumInts2Bit(gpuReadStorage->getSequenceLengthUpperBound());
            qualityPitchInBytes = SDIV(gpuReadStorage->getSequenceLengthUpperBound(), 32) * 32;
        }

        ~GpuAnchorHasher(){
            // std::cerr << "GpuAnchorHasher::~GpuAnchorHasher(). Memory of minhash handle: ";
            // auto memoryUsage = gpuMinhasher->getMemoryInfo(minhashHandle);
            // std::cerr << memoryUsage.host;
            // for(auto pair : memoryUsage.device){
            //     std::cerr << ", [" << pair.first << "] " << pair.second;
            // }
            // std::cerr << "\n";

            gpuReadStorage->destroyHandle(readstorageHandle);
            gpuMinhasher->destroyHandle(minhashHandle);
        }

        void makeErrorCorrectorInput(
            const read_number* anchorIds,
            int numIds,
            bool useQualityScores,
            GpuErrorCorrectorInput& ecinput,
            cudaStream_t stream
        ){
            cub::SwitchDevice sd{deviceId};

            assert(cudaSuccess == ecinput.event.query());
            CUDACHECK(previousBatchFinishedEvent.synchronize());

            resizeBuffers(ecinput, numIds, stream);
    
            //copy input to pinned memory
            *ecinput.h_numAnchors.data() = numIds;
            std::copy_n(anchorIds, numIds, ecinput.h_anchorReadIds.data());

            CUDACHECK(cudaMemcpyAsync(
                ecinput.d_numAnchors.data(),
                ecinput.h_numAnchors.data(),
                sizeof(int),
                H2D,
                stream
            ));

            CUDACHECK(cudaMemcpyAsync(
                ecinput.d_anchorReadIds.data(),
                ecinput.h_anchorReadIds.data(),
                sizeof(read_number) * (*ecinput.h_numAnchors.data()),
                H2D,
                stream
            ));

            if(numIds > 0){
                nvtx::push_range("getAnchorReads", 0);
                getAnchorReads(ecinput, useQualityScores, stream);
                nvtx::pop_range();

                nvtx::push_range("getCandidateReadIdsWithMinhashing", 1);
                getCandidateReadIdsWithMinhashing(ecinput, stream);
                nvtx::pop_range();

                CUDACHECK(cudaStreamSynchronize(stream));

                getCandidateReads(ecinput, useQualityScores, stream);
            }            

            CUDACHECK(ecinput.event.record(stream));
            CUDACHECK(previousBatchFinishedEvent.record(stream));
        }

        MemoryUsage getMemoryInfo() const{
            MemoryUsage info{};
       
            info += gpuMinhasher->getMemoryInfo(minhashHandle);
          
            info += gpuReadStorage->getMemoryInfo(readstorageHandle);
            return info;
        } 

    public: //private:
        void resizeBuffers(GpuErrorCorrectorInput& ecinput, int numAnchors, cudaStream_t stream){
            ecinput.h_numAnchors.resize(1);
            ecinput.h_numCandidates.resize(1);
            ecinput.h_anchorReadIds.resize(numAnchors);

            ecinput.d_numAnchors.resize(1, stream);
            ecinput.d_numCandidates.resize(1, stream);
            ecinput.d_anchorReadIds.resize(numAnchors, stream);
            ecinput.d_anchor_sequences_data.resize(encodedSequencePitchInInts * numAnchors, stream);
            ecinput.d_anchor_sequences_lengths.resize(numAnchors, stream);
            ecinput.d_candidates_per_anchor.resize(numAnchors, stream);
            ecinput.d_candidates_per_anchor_prefixsum.resize(numAnchors + 1, stream);
        }
        
        void getAnchorReads(GpuErrorCorrectorInput& ecinput, bool /*useQualityScores*/, cudaStream_t stream){
            const int numAnchors = (*ecinput.h_numAnchors.data());

            gpuReadStorage->gatherContiguousSequences(
                readstorageHandle,
                ecinput.d_anchor_sequences_data.data(),
                encodedSequencePitchInInts,
                ecinput.h_anchorReadIds[0],
                numAnchors,
                stream,
                mr
            );

            gpuReadStorage->gatherSequenceLengths(
                readstorageHandle,
                ecinput.d_anchor_sequences_lengths.data(),
                ecinput.d_anchorReadIds.data(),
                numAnchors,
                stream
            );
        }

        void getCandidateReads(GpuErrorCorrectorInput& ecinput, bool /*useQualityScores*/, cudaStream_t stream){
            const int numCandidates = (*ecinput.h_numCandidates.data());

            ecinput.d_candidate_sequences_data.resize(encodedSequencePitchInInts * numCandidates, stream);

            gpuReadStorage->gatherSequences(
                readstorageHandle,
                ecinput.d_candidate_sequences_data.data(),
                encodedSequencePitchInInts,
                makeAsyncConstBufferWrapper(ecinput.h_candidate_read_ids.data()),
                ecinput.d_candidate_read_ids.data(),
                numCandidates,
                stream,
                mr
            );

            ecinput.d_candidate_sequences_lengths.resize(numCandidates, stream);

            gpuReadStorage->gatherSequenceLengths(
                readstorageHandle,
                ecinput.d_candidate_sequences_lengths.data(),
                ecinput.d_candidate_read_ids.data(),
                numCandidates,
                stream
            );
        }

        void getCandidateReadIdsWithMinhashing(GpuErrorCorrectorInput& ecinput, cudaStream_t stream){
            int totalNumValues = 0;

            gpuMinhasher->determineNumValues(
                minhashHandle,
                ecinput.d_anchor_sequences_data.data(),
                encodedSequencePitchInInts,
                ecinput.d_anchor_sequences_lengths.data(),
                (*ecinput.h_numAnchors.data()),
                ecinput.d_candidates_per_anchor.data(),
                totalNumValues,
                stream,
                mr
            );

            CUDACHECK(cudaStreamSynchronize(stream));

            ecinput.d_candidate_read_ids.resize(totalNumValues, stream);
            ecinput.h_candidate_read_ids.resize(totalNumValues);

            if(totalNumValues == 0){
                *ecinput.h_numCandidates = 0;
                CUDACHECK(cudaMemsetAsync(ecinput.d_numCandidates.data(), 0, sizeof(int), stream));
                CUDACHECK(cudaMemsetAsync(ecinput.d_candidates_per_anchor.data(), 0, sizeof(int) * (*ecinput.h_numAnchors), stream));
                CUDACHECK(cudaMemsetAsync(ecinput.d_candidates_per_anchor_prefixsum.data(), 0, sizeof(int) * (1 + (*ecinput.h_numAnchors)), stream));
                return;
            }

            gpuMinhasher->retrieveValues(
                minhashHandle,
                (*ecinput.h_numAnchors.data()),                
                totalNumValues,
                ecinput.d_candidate_read_ids.data(),
                ecinput.d_candidates_per_anchor.data(),
                ecinput.d_candidates_per_anchor_prefixsum.data(),
                stream,
                mr
            );

            rmm::device_uvector<read_number> d_candidate_read_ids2(totalNumValues, stream, mr);
            rmm::device_uvector<int> d_candidates_per_anchor2((*ecinput.h_numAnchors), stream, mr);
            rmm::device_uvector<int> d_candidates_per_anchor_prefixsum2(1 + (*ecinput.h_numAnchors), stream, mr);

            cub::DoubleBuffer<read_number> d_items{ecinput.d_candidate_read_ids.data(), d_candidate_read_ids2.data()};
            cub::DoubleBuffer<int> d_numItemsPerSegment{ecinput.d_candidates_per_anchor.data(), d_candidates_per_anchor2.data()};
            cub::DoubleBuffer<int> d_numItemsPerSegmentPrefixSum{ecinput.d_candidates_per_anchor_prefixsum.data(), d_candidates_per_anchor_prefixsum2.data()};

            GpuMinhashQueryFilter::keepDistinctAndNotMatching(
                ecinput.d_anchorReadIds.data(),
                d_items,
                d_numItemsPerSegment,
                d_numItemsPerSegmentPrefixSum, //numSegments + 1
                (*ecinput.h_numAnchors),
                totalNumValues,
                stream,
                mr
            );

            if(d_items.Current() != ecinput.d_candidate_read_ids.data()){
                //std::cerr << "swap d_candidate_read_ids\n";
                std::swap(ecinput.d_candidate_read_ids, d_candidate_read_ids2);
            }
            if(d_numItemsPerSegment.Current() != ecinput.d_candidates_per_anchor.data()){
                //std::cerr << "swap d_candidates_per_anchor\n";
                std::swap(ecinput.d_candidates_per_anchor, d_candidates_per_anchor2);
            }
            if(d_numItemsPerSegmentPrefixSum.Current() != ecinput.d_candidates_per_anchor_prefixsum.data()){
                //std::cerr << "swap d_candidates_per_anchor_prefixsum\n";
                std::swap(ecinput.d_candidates_per_anchor_prefixsum, d_candidates_per_anchor_prefixsum2);
            }

            gpucorrectorkernels::copyMinhashResultsKernel<<<640, 256, 0, stream>>>(
                ecinput.d_numCandidates.data(),
                ecinput.h_numCandidates.data(),
                ecinput.h_candidate_read_ids.data(),
                ecinput.d_candidates_per_anchor_prefixsum.data(),
                ecinput.d_candidate_read_ids.data(),
                *ecinput.h_numAnchors.data()
            ); CUDACHECKASYNC;

        }
    
        int deviceId;
        int maxCandidatesPerRead;
        std::size_t encodedSequencePitchInInts;
        std::size_t qualityPitchInBytes;
        CudaEvent previousBatchFinishedEvent;
        const GpuReadStorage* gpuReadStorage;
        const GpuMinhasher* gpuMinhasher;
        ThreadPool* threadPool;
        ThreadPool::ParallelForHandle pforHandle;
        MinhasherHandle minhashHandle;
        ReadStorageHandle readstorageHandle;
        rmm::mr::device_memory_resource* mr;
    };


    class OutputConstructor{
    public:

        OutputConstructor() = default;

        OutputConstructor(
            ReadCorrectionFlags& correctionFlags_,
            const ProgramOptions& programOptions_
        ) :
            correctionFlags{&correctionFlags_},
            programOptions{&programOptions_}
        {

        }

        std::vector<int> getAnchorIndicesToProcessAndUpdateCorrectionFlags(const GpuErrorCorrectorRawOutput& currentOutput) const{
            std::vector<int> anchorIndicesToProcess;
            anchorIndicesToProcess.reserve(currentOutput.numAnchors);

            nvtx::push_range("preprocess anchor results",0);

            for(int anchor_index = 0; anchor_index < currentOutput.numAnchors; anchor_index++){
                const read_number readId = currentOutput.h_anchorReadIds[anchor_index];
                const bool isCorrected = currentOutput.h_anchor_is_corrected[anchor_index];
                const bool isHQ = currentOutput.h_is_high_quality_anchor[anchor_index].hq();

                if(isHQ){
                    correctionFlags->setCorrectedAsHqAnchor(readId);
                }

                if(isCorrected){
                    anchorIndicesToProcess.emplace_back(anchor_index);
                }else{
                    correctionFlags->setCouldNotBeCorrectedAsAnchor(readId);
                }

                assert(!(isHQ && !isCorrected));
            }

            nvtx::pop_range();

            return anchorIndicesToProcess;
        }

        std::vector<std::pair<int,int>> getCandidateIndicesToProcess(const GpuErrorCorrectorRawOutput& currentOutput) const{
            std::vector<std::pair<int,int>> candidateIndicesToProcess;

            if(programOptions->correctCandidates){
                candidateIndicesToProcess.reserve(16 * currentOutput.numAnchors);
            }

            if(programOptions->correctCandidates){

                nvtx::push_range("preprocess candidate results",0);

                for(int anchor_index = 0; anchor_index < currentOutput.numAnchors; anchor_index++){

                    const int globalOffset = currentOutput.h_num_corrected_candidates_per_anchor_prefixsum[anchor_index];
                    const int n_corrected_candidates = currentOutput.h_num_corrected_candidates_per_anchor[anchor_index];

                    // const int* const my_indices_of_corrected_candidates = currentOutput.h_indices_of_corrected_candidates
                    //                                     + globalOffset;

                    for(int i = 0; i < n_corrected_candidates; ++i) {
                        //const int global_candidate_index = my_indices_of_corrected_candidates[i];
                        //const read_number candidate_read_id = currentOutput.h_candidate_read_ids[global_candidate_index];
                        const read_number candidate_read_id = currentOutput.h_candidate_read_ids[globalOffset + i];

                        if (!correctionFlags->isCorrectedAsHQAnchor(candidate_read_id)) {
                            candidateIndicesToProcess.emplace_back(std::make_pair(anchor_index, i));
                        }
                    }
                }

                nvtx::pop_range();

            }

            return candidateIndicesToProcess;
        }

        template<class ForLoop>
        CorrectionOutput constructResults(const GpuErrorCorrectorRawOutput& currentOutput, ForLoop loopExecutor) const{
            //assert(cudaSuccess == currentOutput.event.query());

            if(currentOutput.nothingToDo){
                return CorrectionOutput{};
            }

            const std::vector<int> anchorIndicesToProcess = getAnchorIndicesToProcessAndUpdateCorrectionFlags(currentOutput);
            const std::vector<std::pair<int,int>> candidateIndicesToProcess = getCandidateIndicesToProcess(currentOutput);

            const int numCorrectedAnchors = anchorIndicesToProcess.size();
            const int numCorrectedCandidates = candidateIndicesToProcess.size();

            // std::cerr << "numCorrectedAnchors: " << numCorrectedAnchors << 
            //     ", numCorrectedCandidates: " << numCorrectedCandidates << "\n";

            CorrectionOutput correctionOutput;
            correctionOutput.anchorCorrections.resize(numCorrectedAnchors);

            if(programOptions->correctCandidates){
                correctionOutput.candidateCorrections.resize(numCorrectedCandidates);
            }

            auto unpackAnchors = [&](int begin, int end){
                nvtx::push_range("Anchor unpacking " + std::to_string(end - begin), 3);

                //Edits and numEdits are stored compact, only for corrected anchors.
                //they are indexed by positionInVector instead of anchor_index
                            
                for(int positionInVector = begin; positionInVector < end; ++positionInVector) {
                    const int anchor_index = anchorIndicesToProcess[positionInVector];

                    auto& tmp = correctionOutput.anchorCorrections[positionInVector];
                    
                    const read_number readId = currentOutput.h_anchorReadIds[anchor_index];

                    tmp.hq = currentOutput.h_is_high_quality_anchor[anchor_index].hq();                    
                    tmp.type = TempCorrectedSequenceType::Anchor;
                    tmp.readId = readId;
                    tmp.edits.clear();
                    
                    const int numEdits = currentOutput.h_numEditsPerCorrectedanchor[positionInVector];
                    if(numEdits != currentOutput.doNotUseEditsValue){
                        const int editOffset = currentOutput.h_anchorEditOffsets[positionInVector];
                        const auto* myedits = currentOutput.h_editsPerCorrectedanchor + editOffset;
                        tmp.edits.insert(tmp.edits.end(), myedits, myedits + numEdits);
                        tmp.useEdits = true;
                    }else{
                        
                        tmp.useEdits = false;

                        const int sequenceOffset = currentOutput.h_correctedAnchorsOffsets[positionInVector];

                        const char* const my_corrected_anchor_data = currentOutput.h_corrected_anchors + sequenceOffset;
                        const int anchor_length = currentOutput.h_anchor_sequences_lengths[anchor_index];
                        tmp.sequence.assign(my_corrected_anchor_data, anchor_length);
                    }

                    // if(tmp.readId == 9273463){
                    //     std::cerr << tmp << "\n";
                    // }
                }

                nvtx::pop_range();
            };

            auto unpackcandidates = [&](int begin, int end){
                nvtx::push_range("candidate unpacking " + std::to_string(end - begin), 3);

                //buffers are stored compact. offsets for each anchor are given by h_num_corrected_candidates_per_anchor_prefixsum
                //Edits, numEdits, h_candidate_read_ids, h_candidate_sequences_lengths, h_alignment_shifts are stored compact, only for corrected candidates.
                //edits are only present for candidates which use edits and have numEdits > 0
                //offsets to the edits of candidates are stored in h_candidateEditOffsets
                

                for(int positionInVector = begin; positionInVector < end; ++positionInVector) {
                    

                    //TIMERSTARTCPU(setup);
                    const int anchor_index = candidateIndicesToProcess[positionInVector].first;
                    const int candidateIndex = candidateIndicesToProcess[positionInVector].second;
                    const read_number anchorReadId = currentOutput.h_anchorReadIds[anchor_index];

                    auto& tmp = correctionOutput.candidateCorrections[positionInVector];

                    const size_t offsetForCorrectedCandidateData = currentOutput.h_num_corrected_candidates_per_anchor_prefixsum[anchor_index];

                    const read_number candidate_read_id = currentOutput.h_candidate_read_ids[offsetForCorrectedCandidateData + candidateIndex];
                    const int candidate_shift = currentOutput.h_alignment_shifts[offsetForCorrectedCandidateData + candidateIndex];

                    if(programOptions->new_columns_to_correct < candidate_shift){
                        std::cerr << "readid " << anchorReadId << " candidate readid " << candidate_read_id << " : "
                        << candidate_shift << " " << programOptions->new_columns_to_correct <<"\n";

                        assert(programOptions->new_columns_to_correct >= candidate_shift);
                    }                
                    
                    tmp.type = TempCorrectedSequenceType::Candidate;
                    tmp.shift = candidate_shift;
                    tmp.readId = candidate_read_id;
                    tmp.edits.clear();
                    
                    const int numEdits = currentOutput.h_numEditsPerCorrectedCandidate[offsetForCorrectedCandidateData + candidateIndex];
                    const int editsOffset = currentOutput.h_candidateEditOffsets[offsetForCorrectedCandidateData + candidateIndex];

                    if(numEdits != currentOutput.doNotUseEditsValue){
                        const auto* myEdits = &currentOutput.h_editsPerCorrectedCandidate[editsOffset];
                        tmp.edits.insert(tmp.edits.end(), myEdits, myEdits + numEdits);
                        tmp.useEdits = true;
                    }else{
                        const int correctionOffset = currentOutput.h_correctedCandidatesOffsets[candidateIndex];
                        const int candidate_length = currentOutput.h_candidate_sequences_lengths[candidateIndex];
                        const char* const candidate_data = currentOutput.h_corrected_candidates + correctionOffset * currentOutput.decodedSequencePitchInBytes;
                        tmp.sequence.assign(candidate_data, candidate_length);
                        
                        tmp.useEdits = false;
                    }

                    // if(tmp.readId == 9273463){
                    //     std::cerr << tmp << " with anchorid " << anchorReadId << "\n";
                    // }
                }

                nvtx::pop_range();
            };


            if(!programOptions->correctCandidates){
                loopExecutor(
                    0, 
                    numCorrectedAnchors, 
                    [=](auto begin, auto end, auto /*threadId*/){
                        unpackAnchors(begin, end);
                    }
                );
            }else{
        
  
                loopExecutor(
                    0, 
                    numCorrectedAnchors, 
                    [=](auto begin, auto end, auto /*threadId*/){
                        unpackAnchors(begin, end);
                    }
                );

                loopExecutor(
                    0, 
                    numCorrectedCandidates, 
                    [=](auto begin, auto end, auto /*threadId*/){
                        unpackcandidates(begin, end);
                    }
                );      
            }

            return correctionOutput;
        }

        template<class ForLoop>
        EncodedCorrectionOutput constructEncodedResults(const GpuErrorCorrectorRawOutput& currentOutput, ForLoop loopExecutor) const{
            //assert(cudaSuccess == currentOutput.event.query());

            if(currentOutput.nothingToDo){
                return EncodedCorrectionOutput{};
            }

            const std::vector<int> anchorIndicesToProcess = getAnchorIndicesToProcessAndUpdateCorrectionFlags(currentOutput);
            const std::vector<std::pair<int,int>> candidateIndicesToProcess = getCandidateIndicesToProcess(currentOutput);

            const int numCorrectedAnchors = anchorIndicesToProcess.size();
            const int numCorrectedCandidates = candidateIndicesToProcess.size();

            // std::cerr << "numCorrectedAnchors: " << numCorrectedAnchors << 
            //     ", numCorrectedCandidates: " << numCorrectedCandidates << "\n";

            EncodedCorrectionOutput encodedCorrectionOutput;
            encodedCorrectionOutput.encodedAnchorCorrections.resize(numCorrectedAnchors);

            if(programOptions->correctCandidates){
                encodedCorrectionOutput.encodedCandidateCorrections.resize(numCorrectedCandidates);
            }

            auto unpackAnchors = [&](int begin, int end){
                nvtx::push_range("Anchor unpacking " + std::to_string(end - begin), 3);

                //Edits and numEdits are stored compact, only for corrected anchors.
                //they are indexed by positionInVector instead of anchor_index

                std::vector<CorrectionEdit> edits;
                            
                for(int positionInVector = begin; positionInVector < end; ++positionInVector) {
                    const int anchor_index = anchorIndicesToProcess[positionInVector];
                    
                    const read_number readId = currentOutput.h_anchorReadIds[anchor_index];

                    edits.clear();
                    bool useEdits = false;
                    const char* sequence = nullptr;
                    int sequenceLength = 0;
                    
                    const int numEdits = currentOutput.h_numEditsPerCorrectedanchor[positionInVector];
                    if(numEdits != currentOutput.doNotUseEditsValue){
                        const int editOffset = currentOutput.h_anchorEditOffsets[positionInVector];
                        const auto* myEdits = currentOutput.h_editsPerCorrectedanchor + editOffset;
                        edits.insert(edits.end(), myEdits, myEdits + numEdits);
                        useEdits = true;
                    }else{                        
                        const int sequenceOffset = currentOutput.h_correctedAnchorsOffsets[positionInVector];
                        const char* const my_corrected_anchor_data = currentOutput.h_corrected_anchors + sequenceOffset;
                        const int anchor_length = currentOutput.h_anchor_sequences_lengths[anchor_index];
 
                        sequenceLength = anchor_length;
                        sequence = my_corrected_anchor_data;
                    }

                    EncodedTempCorrectedSequence::encodeDataIntoEncodedCorrectedSequence(
                        encodedCorrectionOutput.encodedAnchorCorrections[positionInVector],
                        readId,
                        currentOutput.h_is_high_quality_anchor[anchor_index].hq(),
                        useEdits,
                        TempCorrectedSequenceType::Anchor,
                        0,
                        edits.size(),
                        edits.data(),
                        sequenceLength,
                        sequence
                    );
                }

                nvtx::pop_range();
            };

            auto unpackcandidates = [&](int begin, int end){
                nvtx::push_range("candidate unpacking " + std::to_string(end - begin), 3);

                //buffers are stored compact. offsets for each anchor are given by h_num_corrected_candidates_per_anchor_prefixsum
                //Edits, numEdits, h_candidate_read_ids, h_candidate_sequences_lengths, h_alignment_shifts are stored compact, only for corrected candidates.
                //edits are only present for candidates which use edits and have numEdits > 0
                //offsets to the edits of candidates are stored in h_candidateEditOffsets

                std::vector<CorrectionEdit> edits;          

                for(int positionInVector = begin; positionInVector < end; ++positionInVector) {
                    

                    const int anchor_index = candidateIndicesToProcess[positionInVector].first;
                    const int candidateIndex = candidateIndicesToProcess[positionInVector].second;
                    const read_number anchorReadId = currentOutput.h_anchorReadIds[anchor_index];

                    const size_t offsetForCorrectedCandidateData = currentOutput.h_num_corrected_candidates_per_anchor_prefixsum[anchor_index];

                    const read_number candidate_read_id = currentOutput.h_candidate_read_ids[offsetForCorrectedCandidateData + candidateIndex];
                    const int candidate_shift = currentOutput.h_alignment_shifts[offsetForCorrectedCandidateData + candidateIndex];

                    if(programOptions->new_columns_to_correct < candidate_shift){
                        std::cerr << "readid " << anchorReadId << " candidate readid " << candidate_read_id << " : "
                        << candidate_shift << " " << programOptions->new_columns_to_correct <<"\n";

                        assert(programOptions->new_columns_to_correct >= candidate_shift);
                    }

                    edits.clear();
                    
                    const int numEdits = currentOutput.h_numEditsPerCorrectedCandidate[offsetForCorrectedCandidateData + candidateIndex];
                    const int editsOffset = currentOutput.h_candidateEditOffsets[offsetForCorrectedCandidateData + candidateIndex];

                    bool useEdits = false;
                    const char* sequence = nullptr;
                    int sequenceLength = 0;

                    if(numEdits != currentOutput.doNotUseEditsValue){
                        const auto* myEdits = &currentOutput.h_editsPerCorrectedCandidate[editsOffset];
                        edits.insert(edits.end(), myEdits, myEdits + numEdits);
                        useEdits = true;
                    }else{
                        const int correctionOffset = currentOutput.h_correctedCandidatesOffsets[candidateIndex];
                        const int candidate_length = currentOutput.h_candidate_sequences_lengths[candidateIndex];
                        const char* const candidate_data = currentOutput.h_corrected_candidates + correctionOffset * currentOutput.decodedSequencePitchInBytes;

                        sequenceLength = candidate_length;
                        sequence = candidate_data;
                    }

                    EncodedTempCorrectedSequence::encodeDataIntoEncodedCorrectedSequence(
                        encodedCorrectionOutput.encodedCandidateCorrections[positionInVector],
                        candidate_read_id,
                        false,
                        useEdits,
                        TempCorrectedSequenceType::Candidate,
                        candidate_shift,
                        edits.size(),
                        edits.data(),
                        sequenceLength,
                        sequence
                    );
                }

                nvtx::pop_range();
            };


            if(!programOptions->correctCandidates){
                loopExecutor(
                    0, 
                    numCorrectedAnchors, 
                    [=](auto begin, auto end, auto /*threadId*/){
                        unpackAnchors(begin, end);
                    }
                );
            }else{
        
  
                loopExecutor(
                    0, 
                    numCorrectedAnchors, 
                    [=](auto begin, auto end, auto /*threadId*/){
                        unpackAnchors(begin, end);
                    }
                );

                loopExecutor(
                    0, 
                    numCorrectedCandidates, 
                    [=](auto begin, auto end, auto /*threadId*/){
                        unpackcandidates(begin, end);
                    }
                );      
            }

            return encodedCorrectionOutput;
        }

        MemoryUsage getMemoryInfo() const{
            MemoryUsage info{};
            return info;
        }
    public: //private:
        ReadCorrectionFlags* correctionFlags;
        const ProgramOptions* programOptions;
    };

    class GpuErrorCorrector{

    public:

        template<class T>
        using PinnedBuffer = helpers::SimpleAllocationPinnedHost<T>;

        static constexpr int getNumRefinementIterations() noexcept{
            return 5;
        }

        static constexpr bool useMsaRefinement() noexcept{
            return getNumRefinementIterations() > 0;
        }

        GpuErrorCorrector() = default;

        GpuErrorCorrector(
            const GpuReadStorage& gpuReadStorage_,
            const ReadCorrectionFlags& correctionFlags_,
            const ProgramOptions& programOptions_,
            int maxAnchorsPerCall,
            rmm::mr::device_memory_resource* mr_,
            ThreadPool* threadPool_,
            const GpuForest* gpuForestAnchor_,
            const GpuForest* gpuForestCandidate_
        ) : 
            maxAnchors{maxAnchorsPerCall},
            correctionFlags{&correctionFlags_},
            gpuReadStorage{&gpuReadStorage_},
            programOptions{&programOptions_},
            mr{mr_},
            threadPool{threadPool_},
            gpuForestAnchor{gpuForestAnchor_},
            gpuForestCandidate{gpuForestCandidate_},
            readstorageHandle{gpuReadStorage->makeHandle()},
            d_indicesForGather{0, cudaStreamPerThread, mr},
            d_anchorContainsN{0, cudaStreamPerThread, mr},
            d_candidateContainsN{0, cudaStreamPerThread, mr},
            d_candidate_sequences_lengths{0, cudaStreamPerThread, mr},
            d_candidate_sequences_data{0, cudaStreamPerThread, mr},
            d_anchorIndicesOfCandidates{0, cudaStreamPerThread, mr},
            d_alignment_overlaps{0, cudaStreamPerThread, mr},
            d_alignment_shifts{0, cudaStreamPerThread, mr},
            d_alignment_nOps{0, cudaStreamPerThread, mr},
            d_alignment_best_alignment_flags{0, cudaStreamPerThread, mr}, 
            d_indices{0, cudaStreamPerThread, mr},
            d_indices_per_anchor{0, cudaStreamPerThread, mr},
            d_indices_per_anchor_prefixsum{0, cudaStreamPerThread, mr},
            d_num_indices{0, cudaStreamPerThread, mr},
            d_corrected_anchors{0, cudaStreamPerThread, mr},
            d_corrected_candidates{0, cudaStreamPerThread, mr},
            d_num_corrected_candidates_per_anchor{0, cudaStreamPerThread, mr},
            d_num_corrected_candidates_per_anchor_prefixsum{0, cudaStreamPerThread, mr},
            d_num_total_corrected_candidates{0, cudaStreamPerThread, mr},
            d_anchor_is_corrected{0, cudaStreamPerThread, mr},
            d_is_high_quality_anchor{0, cudaStreamPerThread, mr},
            d_high_quality_anchor_indices{0, cudaStreamPerThread, mr},
            d_num_high_quality_anchor_indices{0, cudaStreamPerThread, mr}, 
            d_editsPerCorrectedanchor{0, cudaStreamPerThread, mr},
            d_numEditsPerCorrectedanchor{0, cudaStreamPerThread, mr},
            d_editsPerCorrectedCandidate{0, cudaStreamPerThread, mr},
            d_numEditsPerCorrectedCandidate{0, cudaStreamPerThread, mr},
            d_indices_of_corrected_anchors{0, cudaStreamPerThread, mr},
            d_num_indices_of_corrected_anchors{0, cudaStreamPerThread, mr},
            d_indices_of_corrected_candidates{0, cudaStreamPerThread, mr},
            d_totalNumEdits{0, cudaStreamPerThread, mr},
            d_isPairedCandidate{0, cudaStreamPerThread, mr},
            d_numAnchors{0, cudaStreamPerThread, mr},
            d_numCandidates{0, cudaStreamPerThread, mr},
            d_anchorReadIds{0, cudaStreamPerThread, mr},
            d_anchor_sequences_data{0, cudaStreamPerThread, mr},
            d_anchor_sequences_lengths{0, cudaStreamPerThread, mr},
            d_candidate_read_ids{0, cudaStreamPerThread, mr},
            d_candidates_per_anchor{0, cudaStreamPerThread, mr},
            d_candidates_per_anchor_prefixsum{0, cudaStreamPerThread, mr}
        {
            if(programOptions->correctionType != CorrectionType::Classic){
                assert(gpuForestAnchor != nullptr);
            }

            if(programOptions->correctionTypeCands != CorrectionType::Classic){
                assert(gpuForestCandidate != nullptr);
            }

            CUDACHECK(cudaGetDevice(&deviceId));

            for(auto& event: events){
                event = std::move(CudaEvent{cudaEventDisableTiming});
            }

            inputCandidateDataIsReadyEvent = CudaEvent{cudaEventDisableTiming};
            previousBatchFinishedEvent = CudaEvent{cudaEventDisableTiming};

            encodedSequencePitchInInts = SequenceHelpers::getEncodedNumInts2Bit(gpuReadStorage->getSequenceLengthUpperBound());
            decodedSequencePitchInBytes = SDIV(gpuReadStorage->getSequenceLengthUpperBound(), 4) * 4;
            qualityPitchInBytes = SDIV(gpuReadStorage->getSequenceLengthUpperBound(), 32) * 32;
            maxNumEditsPerSequence = std::max(1,gpuReadStorage->getSequenceLengthUpperBound() / 7);
            //pad to multiple of 128 bytes
            editsPitchInBytes = SDIV(maxNumEditsPerSequence * sizeof(EncodedCorrectionEdit), 128) * 128;

            const std::size_t min_overlap = std::max(
                1, 
                std::max(
                    programOptions->min_overlap, 
                    int(gpuReadStorage->getSequenceLengthUpperBound() * programOptions->min_overlap_ratio)
                )
            );
            const std::size_t msa_max_column_count = (3*gpuReadStorage->getSequenceLengthUpperBound() - 2*min_overlap);
            //round up to 32 elements
            msaColumnPitchInElements = SDIV(msa_max_column_count, 32) * 32;

            initFixedSizeBuffers(cudaStreamPerThread);

            CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));
            extraStream = streampool::get_current_device_pool()->get_stream();
        }

        ~GpuErrorCorrector(){
            gpuReadStorage->destroyHandle(readstorageHandle);
        }

        void correct(GpuErrorCorrectorInput& input, GpuErrorCorrectorRawOutput& output, cudaStream_t stream){
            CUDACHECK(previousBatchFinishedEvent.synchronize());
            CUDACHECK(cudaStreamSynchronize(stream));

            //assert(cudaSuccess == input.event.query());
            //assert(cudaSuccess == output.event.query());

            currentInput = &input;
            currentOutput = &output;

            assert(*currentInput->h_numAnchors.data() <= maxAnchors);

            currentNumAnchors = *currentInput->h_numAnchors.data();
            currentNumCandidates = *currentInput->h_numCandidates.data();

            if(gpuReadStorage->isPairedEnd()){
                assert(currentNumAnchors % 2 == 0);
            }

            currentOutput->nothingToDo = false;
            currentOutput->numAnchors = currentNumAnchors;
            currentOutput->numCandidates = currentNumCandidates;
            currentOutput->doNotUseEditsValue = getDoNotUseEditsValue();
            currentOutput->editsPitchInBytes = editsPitchInBytes;
            currentOutput->decodedSequencePitchInBytes = decodedSequencePitchInBytes;

            if(currentNumCandidates == 0 || currentNumAnchors == 0){
                currentOutput->nothingToDo = true;
                return;
            }
            
            cub::SwitchDevice sd{deviceId};

            //fixed size memory should already be allocated. However, this will also set the correct working stream for stream-ordered allocations which is important.
            initFixedSizeBuffers(stream); 

            resizeBuffers(currentNumAnchors, currentNumCandidates, stream);

            gpucorrectorkernels::copyCorrectionInputDeviceData<<<SDIV(currentNumCandidates, 256),256, 0, stream>>>(
                d_numAnchors.data(),
                d_numCandidates.data(),
                d_anchorReadIds.data(),
                d_anchor_sequences_data.data(),
                d_anchor_sequences_lengths.data(),
                d_candidate_read_ids.data(),
                d_candidates_per_anchor.data(),
                d_candidates_per_anchor_prefixsum.data(),
                encodedSequencePitchInInts,
                currentNumAnchors,
                currentNumCandidates,
                currentInput->d_anchorReadIds.data(),
                currentInput->d_anchor_sequences_data.data(),
                currentInput->d_anchor_sequences_lengths.data(),
                currentInput->d_candidate_read_ids.data(),
                currentInput->d_candidates_per_anchor.data(),
                currentInput->d_candidates_per_anchor_prefixsum.data()
            ); CUDACHECKASYNC;

            CUDACHECK(cudaMemcpyAsync(
                d_candidate_sequences_data.data(),
                currentInput->d_candidate_sequences_data.data(),
                sizeof(unsigned int) * encodedSequencePitchInInts * currentNumCandidates,
                D2D,
                stream
            ));

            CUDACHECK(cudaMemcpyAsync(
                d_candidate_sequences_lengths.data(),
                currentInput->d_candidate_sequences_lengths.data(),
                sizeof(int) * currentNumCandidates,
                D2D,
                stream
            ));

            CUDACHECK(cudaEventRecord(inputCandidateDataIsReadyEvent, stream));

            //after gpu data has been copied to local working set, the gpu data of currentInput can be reused
            CUDACHECK(currentInput->event.record(stream));

            gpucorrectorkernels::setAnchorIndicesOfCandidateskernel
                    <<<currentNumAnchors, 128, 0, stream>>>(
                d_anchorIndicesOfCandidates.data(),
                d_numAnchors.data(),
                d_candidates_per_anchor.data(),
                d_candidates_per_anchor_prefixsum.data()
            ); CUDACHECKASYNC;

            flagPairedCandidates(stream);

            getAmbiguousFlagsOfAnchors(stream);
            getAmbiguousFlagsOfCandidates(stream);


            // nvtx::push_range("getCandidateSequenceData", 3);
            // getCandidateSequenceData(stream); 
            // nvtx::pop_range();

            nvtx::push_range("getCandidateAlignments", 5);
            getCandidateAlignments(stream); 
            nvtx::pop_range();

            nvtx::push_range("buildMultipleSequenceAlignment", 6);
            buildAndRefineMultipleSequenceAlignment(stream);
            nvtx::pop_range();

            nvtx::push_range("correctanchors", 8);
            correctAnchors(stream);
            nvtx::pop_range();
            
            if(programOptions->correctCandidates) {
                CUDACHECK(cudaEventRecord(events[0], stream));
            }
            

            nvtx::push_range("copyAnchorResultsFromDeviceToHost", 3);
            copyAnchorResultsFromDeviceToHost(stream);
            nvtx::pop_range();

            if(programOptions->correctCandidates) {
                cudaStream_t candsStream = extraStream;
                CUDACHECK(cudaStreamWaitEvent(candsStream, events[0], 0)); 

                nvtx::push_range("correctCandidates", 9);
                correctCandidates(candsStream);
                nvtx::pop_range();

                nvtx::push_range("copyCandidateResultsFromDeviceToHost", 4);
                copyCandidateResultsFromDeviceToHost(candsStream);
                nvtx::pop_range(); 

                CUDACHECK(cudaEventRecord(events[0], candsStream)); 
                CUDACHECK(cudaStreamWaitEvent(stream, events[0], 0));      
            }

            managedgpumsa = nullptr;

            std::copy_n(currentInput->h_anchorReadIds.data(), currentNumAnchors, currentOutput->h_anchorReadIds.data());            

            //after the current work in stream is completed, all results in currentOutput are ready to use.
            CUDACHECK(cudaEventRecord(currentOutput->event, stream));
            CUDACHECK(cudaEventRecord(previousBatchFinishedEvent, stream));
        }

        MemoryUsage getMemoryInfo() const{
            MemoryUsage info{};
            auto handleHost = [&](const auto& h){
                info.host += h.sizeInBytes();
            };
            auto handleDevice = [&](const auto& d){
                using ElementType = typename std::remove_reference<decltype(d)>::type::value_type;
                info.device[deviceId] += d.size() * sizeof(ElementType);
            };

            info += gpuReadStorage->getMemoryInfo(readstorageHandle);
            if(managedgpumsa){
                info += managedgpumsa->getMemoryInfo();
            }

            handleHost(h_num_total_corrected_candidates);
            handleHost(h_num_indices);
            handleHost(h_numSelected);
            handleHost(h_numRemainingCandidatesAfterAlignment);
            handleHost(h_managedmsa_tmp);

            handleHost(h_indicesForGather);
            handleHost(h_isPairedCandidate);
            handleHost(h_candidates_per_anchor_prefixsum);
            handleHost(h_indices);
            handleHost(h_flagsCandidates);

            handleDevice(d_anchorContainsN);
            handleDevice(d_candidateContainsN);
            handleDevice(d_candidate_sequences_lengths);
            handleDevice(d_candidate_sequences_data);
            handleDevice(d_anchorIndicesOfCandidates);
            handleDevice(d_alignment_overlaps);
            handleDevice(d_alignment_shifts);
            handleDevice(d_alignment_nOps);
            handleDevice(d_alignment_best_alignment_flags);
            handleDevice(d_indices);
            handleDevice(d_indices_per_anchor);
            handleDevice(d_indices_per_anchor_prefixsum);
            handleDevice(d_num_indices);
            handleDevice(d_corrected_anchors);
            handleDevice(d_corrected_candidates);
            handleDevice(d_num_corrected_candidates_per_anchor);
            handleDevice(d_num_corrected_candidates_per_anchor_prefixsum);
            handleDevice(d_num_total_corrected_candidates);
            handleDevice(d_anchor_is_corrected);
            handleDevice(d_is_high_quality_anchor);
            handleDevice(d_high_quality_anchor_indices);
            handleDevice(d_num_high_quality_anchor_indices);
            handleDevice(d_editsPerCorrectedanchor);
            handleDevice(d_numEditsPerCorrectedanchor);
            handleDevice(d_editsPerCorrectedCandidate);
            handleDevice(d_numEditsPerCorrectedCandidate);
            handleDevice(d_indices_of_corrected_anchors);
            handleDevice(d_num_indices_of_corrected_anchors);
            handleDevice(d_indices_of_corrected_candidates);
            handleDevice(d_numEditsPerCorrectedanchor);
            handleDevice(d_numAnchors);
            handleDevice(d_numCandidates);
            handleDevice(d_anchorReadIds);
            handleDevice(d_anchor_sequences_data);
            handleDevice(d_anchor_sequences_lengths);
            handleDevice(d_candidate_read_ids);
            handleDevice(d_candidates_per_anchor);
            handleDevice(d_candidates_per_anchor_prefixsum);

            return info;
        } 

        void releaseMemory(cudaStream_t stream){
            auto handleDevice = [&](auto& d){
                ::destroy(d, stream);
            };

            handleDevice(d_anchorContainsN);
            handleDevice(d_candidateContainsN);
            handleDevice(d_candidate_sequences_lengths);
            handleDevice(d_candidate_sequences_data);
            handleDevice(d_anchorIndicesOfCandidates);
            handleDevice(d_alignment_overlaps);
            handleDevice(d_alignment_shifts);
            handleDevice(d_alignment_nOps);
            handleDevice(d_alignment_best_alignment_flags);
            handleDevice(d_indices);
            handleDevice(d_indices_per_anchor);
            handleDevice(d_indices_per_anchor_prefixsum);
            handleDevice(d_num_indices);
            handleDevice(d_corrected_anchors);
            handleDevice(d_corrected_candidates);
            handleDevice(d_num_corrected_candidates_per_anchor);
            handleDevice(d_num_corrected_candidates_per_anchor_prefixsum);
            handleDevice(d_num_total_corrected_candidates);
            handleDevice(d_anchor_is_corrected);
            handleDevice(d_is_high_quality_anchor);
            handleDevice(d_high_quality_anchor_indices);
            handleDevice(d_num_high_quality_anchor_indices);
            handleDevice(d_editsPerCorrectedanchor);
            handleDevice(d_numEditsPerCorrectedanchor);
            handleDevice(d_editsPerCorrectedCandidate);
            handleDevice(d_numEditsPerCorrectedCandidate);
            handleDevice(d_indices_of_corrected_anchors);
            handleDevice(d_num_indices_of_corrected_anchors);
            handleDevice(d_indices_of_corrected_candidates);
            handleDevice(d_numEditsPerCorrectedanchor);
            handleDevice(d_numAnchors);
            handleDevice(d_numCandidates);
            handleDevice(d_anchorReadIds);
            handleDevice(d_anchor_sequences_data);
            handleDevice(d_anchor_sequences_lengths);
            handleDevice(d_candidate_read_ids);
            handleDevice(d_candidates_per_anchor);
            handleDevice(d_candidates_per_anchor_prefixsum);
        } 

        void releaseCandidateMemory(cudaStream_t stream){
            auto handleDevice = [&](auto& d){
                ::destroy(d, stream);
            };

            handleDevice(d_candidateContainsN);
            handleDevice(d_candidate_sequences_lengths);
            handleDevice(d_candidate_sequences_data);
            handleDevice(d_anchorIndicesOfCandidates);
            handleDevice(d_alignment_overlaps);
            handleDevice(d_alignment_shifts);
            handleDevice(d_alignment_nOps);
            handleDevice(d_alignment_best_alignment_flags);
            handleDevice(d_indices);
            handleDevice(d_corrected_candidates);
            handleDevice(d_editsPerCorrectedCandidate);
            handleDevice(d_numEditsPerCorrectedCandidate);
            handleDevice(d_indices_of_corrected_candidates);
            handleDevice(d_candidate_read_ids);
            handleDevice(d_candidates_per_anchor);
            handleDevice(d_candidates_per_anchor_prefixsum);
        } 

        


    public: //private:

        void initFixedSizeBuffers(cudaStream_t stream){
            const std::size_t numEditsAnchors = SDIV(editsPitchInBytes * maxAnchors, sizeof(EncodedCorrectionEdit));          

            h_num_total_corrected_candidates.resize(1);
            h_num_indices.resize(1);
            h_numSelected.resize(1);
            h_numRemainingCandidatesAfterAlignment.resize(1);
            h_managedmsa_tmp.resize(1);

            d_anchorContainsN.resize(maxAnchors, stream);

            d_indices_per_anchor.resize(maxAnchors, stream);
            d_num_indices.resize(1, stream);
            d_indices_per_anchor_prefixsum.resize(maxAnchors, stream);
            d_corrected_anchors.resize(maxAnchors * decodedSequencePitchInBytes, stream);
            d_num_corrected_candidates_per_anchor.resize(maxAnchors, stream);
            d_num_corrected_candidates_per_anchor_prefixsum.resize(maxAnchors, stream);
            d_num_total_corrected_candidates.resize(1, stream);
            d_anchor_is_corrected.resize(maxAnchors, stream);
            d_is_high_quality_anchor.resize(maxAnchors, stream);
            d_high_quality_anchor_indices.resize(maxAnchors, stream);
            d_num_high_quality_anchor_indices.resize(1, stream); 
            d_editsPerCorrectedanchor.resize(numEditsAnchors, stream);
            d_numEditsPerCorrectedanchor.resize(maxAnchors, stream);
            d_indices_of_corrected_anchors.resize(maxAnchors, stream);
            d_num_indices_of_corrected_anchors.resize(1, stream);

            d_numAnchors.resize(1, stream);
            d_numCandidates.resize(1, stream);
            d_anchorReadIds.resize(maxAnchors, stream);
            d_anchor_sequences_data.resize(encodedSequencePitchInInts * maxAnchors, stream);
            d_anchor_sequences_lengths.resize(maxAnchors, stream);
            d_candidates_per_anchor.resize(maxAnchors, stream);
            h_candidates_per_anchor_prefixsum.resize(maxAnchors + 1);
            d_candidates_per_anchor_prefixsum.resize(maxAnchors + 1, stream);
            d_totalNumEdits.resize(1, stream);
        }
 
        void resizeBuffers(int /*numReads*/, int numCandidates, cudaStream_t stream){  
            //assert(numReads <= maxAnchors);

            const std::size_t numEditsAnchors = SDIV(editsPitchInBytes * maxAnchors, sizeof(EncodedCorrectionEdit));          

            //does not depend on number of candidates
            currentOutput->h_anchor_sequences_lengths.resize(maxAnchors);
            currentOutput->h_corrected_anchors.resize(maxAnchors * decodedSequencePitchInBytes);            
            currentOutput->h_anchor_is_corrected.resize(maxAnchors);
            currentOutput->h_is_high_quality_anchor.resize(maxAnchors);
            currentOutput->h_editsPerCorrectedanchor.resize(numEditsAnchors);
            currentOutput->h_numEditsPerCorrectedanchor.resize(maxAnchors);            
            currentOutput->h_num_corrected_candidates_per_anchor.resize(maxAnchors);
            currentOutput->h_num_corrected_candidates_per_anchor_prefixsum.resize(maxAnchors);
            currentOutput->h_anchorReadIds.resize(maxAnchors);
            currentOutput->h_anchorEditOffsets.resize(maxAnchors);
            currentOutput->h_correctedAnchorsOffsets.resize(maxAnchors * decodedSequencePitchInBytes);           
            
            d_anchorIndicesOfCandidates.resize(numCandidates, stream);
            d_candidateContainsN.resize(numCandidates, stream);
            d_candidate_read_ids.resize(numCandidates, stream);
            d_candidate_sequences_lengths.resize(numCandidates, stream);
            d_candidate_sequences_data.resize(numCandidates * encodedSequencePitchInInts, stream);
            d_isPairedCandidate.resize(numCandidates, stream);
            h_isPairedCandidate.resize(numCandidates);

            h_flagsCandidates.resize(numCandidates);
            
            h_indicesForGather.resize(numCandidates);
            d_indicesForGather.resize(numCandidates, stream);
            
            d_alignment_overlaps.resize(numCandidates, stream);
            d_alignment_shifts.resize(numCandidates, stream);
            d_alignment_nOps.resize(numCandidates, stream);
            d_alignment_best_alignment_flags.resize(numCandidates, stream);
            d_indices.resize(numCandidates + 1, stream);

            d_indices_of_corrected_candidates.resize(numCandidates, stream);
        }

        void flagPairedCandidates(cudaStream_t stream){

            if(gpuReadStorage->isPairedEnd()){

                assert(currentNumAnchors % 2 == 0);
                assert(currentNumAnchors != 0);

                d_isPairedCandidate.resize(currentNumCandidates, stream);

                helpers::call_fill_kernel_async(d_isPairedCandidate.data(), currentNumCandidates, false, stream);                   

                dim3 block = 128;
                dim3 grid = currentNumAnchors / 2;
                constexpr int staticSmemBytes = 4096;

                gpucorrectorkernels::flagPairedCandidatesKernel<128,staticSmemBytes>
                <<<grid, block, 0, stream>>>(
                    currentNumAnchors / 2,
                    d_candidates_per_anchor.data(),
                    d_candidates_per_anchor_prefixsum.data(),
                    d_candidate_read_ids.data(),
                    d_isPairedCandidate.data()
                ); CUDACHECKASYNC;

                #if 0
                    //remove candidates which are not paired
                    rmm::device_uvector<read_number> d_candidate_read_ids2(currentNumCandidates, stream, mr);
                    rmm::device_uvector<int> d_anchorIndicesOfCandidates2(currentNumCandidates, stream, mr);

                    CubCallWrapper(mr).cubSelectFlagged(
                        thrust::make_zip_iterator(thrust::make_tuple(
                            d_candidate_read_ids.data(),
                            d_anchorIndicesOfCandidates.data()
                        )),                        
                        thrust::make_transform_iterator(d_isPairedCandidate.data(), thrust::identity<bool>()),
                        thrust::make_zip_iterator(thrust::make_tuple(
                            d_candidate_read_ids2.data(),
                            d_anchorIndicesOfCandidates2.data()
                        )),
                        d_numCandidates.data(),
                        currentNumCandidates,
                        stream
                    );

                    CUDACHECK(cudaMemcpyAsync(
                        h_num_indices.data(),
                        d_numCandidates.data(),
                        sizeof(int),
                        D2H,
                        stream
                    ));
                    CUDACHECK(cudaStreamSynchronize(stream));

                    auto oldNumCandidates = currentNumCandidates;
                    currentNumCandidates = *h_num_indices;

                    CUDACHECK(cudaMemcpyAsync(
                        currentInput->h_candidate_read_ids.data(),
                        d_candidate_read_ids2.data(),
                        sizeof(int) * currentNumCandidates,
                        D2H,
                        stream
                    ));
                    CUDACHECK(cudaEventRecord(events[1], stream));

                    std::swap(d_candidate_read_ids, d_candidate_read_ids2);
                    std::swap(d_anchorIndicesOfCandidates, d_anchorIndicesOfCandidates2);

                    CUDACHECK(cudaMemsetAsync(
                        d_candidates_per_anchor.data(),
                        0,
                        sizeof(int) * currentNumAnchors,
                        stream
                    ));


                    if(currentNumCandidates > 0){

                        rmm::device_uvector<int> d_uniqueAnchorIndices(maxNumAnchors, stream, mr);
                        rmm::device_uvector<int> d_aggregates_out(maxNumAnchors, stream, mr);
                        rmm::device_scalar<int> d_numRuns(stream, mr);

                        CubCallWrapper(mr).cubReduceByKey(
                            d_anchorIndicesOfCandidates.data(), 
                            d_uniqueAnchorIndices.data(), 
                            thrust::make_constant_iterator(1), 
                            d_aggregates_out.data(), 
                            d_num_indices.data(), 
                            cub::Sum(), 
                            currentNumCandidates, 
                            stream
                        );

                        helpers::lambda_kernel<<<SDIV(currentNumAnchors, 256), 256, 0, stream>>>(
                            [
                                d_uniqueAnchorIndices = d_uniqueAnchorIndices.data(),
                                d_aggregates_out = d_aggregates_out.data(),
                                d_candidates_per_anchor = d_candidates_per_anchor.data(),
                                d_numRuns = d_numRuns.data()
                            ] __device__ (){
                                
                                const int tid = threadIdx.x + blockIdx.x * blockDim.x;
                                const int stride = blockDim.x * gridDim.x;

                                for(int i = tid; i < *d_numRuns; i += stride){
                                    d_candidates_per_anchor[d_uniqueAnchorIndices[i]]
                                        = d_aggregates_out[i];
                                }
                            }
                        ); CUDACHECKASYNC;

                        CubCallWrapper(mr).cubInclusiveSum(
                            d_candidates_per_anchor.data(),
                            d_candidates_per_anchor_prefixsum.data() + 1,
                            currentNumAnchors,
                            stream
                        );
                    }

                    CUDACHECK(cudaEventSynchronize(events[1])); //wait for currentInput->h_candidateReadIds
                #endif
            }else{
                CUDACHECK(cudaMemsetAsync(
                    d_isPairedCandidate.data(),
                    0,
                    sizeof(bool) * currentNumCandidates,
                    stream
                ));
            }
        }

        void copyAnchorResultsFromDeviceToHost(cudaStream_t stream){
            if(programOptions->correctionType == CorrectionType::Classic){
                copyAnchorResultsFromDeviceToHostClassic(stream);
            }else if(programOptions->correctionType == CorrectionType::Forest){
                copyAnchorResultsFromDeviceToHostForestGpu(stream);
            }else{
                throw std::runtime_error("copyAnchorResultsFromDeviceToHost not implemented for this correctionType");
            }
        }

        void copyAnchorResultsFromDeviceToHostClassic(cudaStream_t stream){

            rmm::device_uvector<int> d_tmpBuffer(2*(currentNumAnchors + 1), stream, mr);
            int* d_editsOffsetsTmp = d_tmpBuffer.data();
            int* d_correctedAnchorOffsetsTmp = d_tmpBuffer.data() + (currentNumAnchors + 1);

            const int* const d_totalNumberOfEdits = d_editsOffsetsTmp + currentNumAnchors;
            const int* const d_totalCorrectedSequencesBytes = d_correctedAnchorOffsetsTmp + currentNumAnchors;

            auto numEditsPerAnchorTransformed = thrust::make_transform_iterator(
                d_numEditsPerCorrectedanchor.data(),
                [doNotUseEditsValue = getDoNotUseEditsValue()] __device__ (const auto& num){ return num == doNotUseEditsValue ? 0 : num;}
            );

            auto correctedAnchorsPitches = thrust::make_transform_iterator(
                d_numEditsPerCorrectedanchor.data(),
                ReplaceNumberOp(getDoNotUseEditsValue(), decodedSequencePitchInBytes)
            );

            //num edits per anchor prefixsum
            //num bytes per corrected anchor sequence prefix sum
            CubCallWrapper(mr).cubInclusiveScan(
                thrust::make_zip_iterator(thrust::make_tuple(
                   numEditsPerAnchorTransformed,
                   correctedAnchorsPitches
                )), 
                thrust::make_zip_iterator(thrust::make_tuple(
                   d_editsOffsetsTmp + 1,
                   d_correctedAnchorOffsetsTmp + 1
                )),
                [] __device__ (const auto& l , const auto& r){
                    //tuple addition
                    return thrust::make_tuple(
                        thrust::get<0>(l) + thrust::get<0>(r),
                        thrust::get<1>(l) + thrust::get<1>(r)
                    );
                },
                currentNumAnchors,
                stream
            );

            h_indicesForGather.resize(2);
            helpers::lambda_kernel<<<1,1,0,stream>>>([
                output = h_indicesForGather.data(), 
                d_totalNumberOfEdits, 
                d_totalCorrectedSequencesBytes
            ]__device__(){
                output[0] = *d_totalNumberOfEdits;
                output[1] = *d_totalCorrectedSequencesBytes;
            }); CUDACHECKASYNC

            CUDACHECK(cudaStreamSynchronize(stream));

            const int& h_totalNumberOfEdits = h_indicesForGather[0];
            const int& h_totalBytesOfCorrectedAnchors = h_indicesForGather[1];

            //copy other buffers to host
            helpers::call_copy_n_kernel(
                thrust::make_zip_iterator(thrust::make_tuple(
                    d_anchor_sequences_lengths.data(), 
                    d_anchor_is_corrected.data(),
                    d_is_high_quality_anchor.data(),
                    d_numEditsPerCorrectedanchor.data()
                )), 
                currentNumAnchors, 
                thrust::make_zip_iterator(thrust::make_tuple(
                    currentOutput->h_anchor_sequences_lengths.data(), 
                    currentOutput->h_anchor_is_corrected.data(),
                    currentOutput->h_is_high_quality_anchor.data(),
                    currentOutput->h_numEditsPerCorrectedanchor.data()
                )), 
                stream
            ); CUDACHECKASYNC;

            if(h_totalNumberOfEdits > 0){
                //compact edits
                rmm::device_uvector<EncodedCorrectionEdit> d_editsPerCorrectedanchor2(h_totalNumberOfEdits + sizeof(int), stream, mr);

                gpucorrectorkernels::compactEditsKernel<<<SDIV(currentNumAnchors, 128), 128, 0, stream>>>(
                    d_editsPerCorrectedanchor.data(),
                    d_editsPerCorrectedanchor2.data(),
                    d_editsOffsetsTmp,
                    d_numAnchors.data(),
                    d_numEditsPerCorrectedanchor.data(),
                    getDoNotUseEditsValue(),
                    editsPitchInBytes
                ); CUDACHECKASYNC;

                CUDACHECK(cudaMemcpyAsync(
                    currentOutput->h_editsPerCorrectedanchor.data(),
                    d_editsPerCorrectedanchor2.data(),
                    sizeof(EncodedCorrectionEdit) * h_totalNumberOfEdits,
                    D2H,
                    stream
                ));
                CUDACHECK(cudaMemcpyAsync(
                    currentOutput->h_anchorEditOffsets.data() + 1,
                    d_editsOffsetsTmp + 1,
                    sizeof(int) * (currentNumAnchors-1),
                    D2H,
                    stream
                ));
                currentOutput->h_anchorEditOffsets[0] = 0;
            }else{
                std::fill_n(currentOutput->h_anchorEditOffsets.data(), currentNumAnchors, 0);
            }

            if(h_totalBytesOfCorrectedAnchors > 0){
                //compact corrected anchor sequences with numEdits == getDoNotUseEditsValue          
                rmm::device_uvector<char> d_corrected_anchors2(SDIV(h_totalBytesOfCorrectedAnchors, sizeof(int)) * sizeof(int), stream, mr);

                gpucorrectorkernels::compactCorrectedSequencesKernel<32><<<SDIV(currentNumAnchors, 128), 128, 0, stream>>>(
                    d_corrected_anchors.data(),
                    d_corrected_anchors2.data(),
                    this->decodedSequencePitchInBytes,
                    d_num_indices_of_corrected_anchors.data(),
                    d_numEditsPerCorrectedanchor.data(),
                    getDoNotUseEditsValue(),
                    d_correctedAnchorOffsetsTmp,
                    d_indices_of_corrected_anchors.data()
                ); CUDACHECKASYNC;

                CUDACHECK(cudaMemcpyAsync(
                    currentOutput->h_corrected_anchors.data(),
                    d_corrected_anchors2.data(),
                    h_totalBytesOfCorrectedAnchors,
                    D2H,
                    stream
                ));

                CUDACHECK(cudaMemcpyAsync(
                    currentOutput->h_correctedAnchorsOffsets.data() + 1,
                    d_correctedAnchorOffsetsTmp + 1,
                    sizeof(int) * (currentNumAnchors-1),
                    D2H,
                    stream
                ));
                currentOutput->h_correctedAnchorsOffsets[0] = 0;
            }else{
                std::fill_n(currentOutput->h_correctedAnchorsOffsets.data(), currentNumAnchors, 0);
            }

        }

        void copyAnchorResultsFromDeviceToHostForestGpu(cudaStream_t stream){
            copyAnchorResultsFromDeviceToHostClassic(stream);
        }

        void copyCandidateResultsFromDeviceToHost(cudaStream_t stream){
            if(programOptions->correctionTypeCands == CorrectionType::Classic){
                copyCandidateResultsFromDeviceToHostClassic(stream);
            }else if(programOptions->correctionTypeCands == CorrectionType::Forest){
                copyCandidateResultsFromDeviceToHostForestGpu(stream);
            }else{
                throw std::runtime_error("copyCandidateResultsFromDeviceToHost not implemented for this correctionTypeCands");
            }
        }

        void copyCandidateResultsFromDeviceToHostClassic(cudaStream_t stream){
            const int numTotalCorrectedCandidates = *h_num_total_corrected_candidates;
            rmm::device_uvector<int> d_tmpBuffer(2*(numTotalCorrectedCandidates + 1), stream, mr);

            int* d_editsOffsetsTmp = d_tmpBuffer.data();
            int* d_correctedCandidatesOffsetsTmp = d_tmpBuffer.data() + numTotalCorrectedCandidates + 1;

            const int* const d_totalNumberOfEdits = d_editsOffsetsTmp + numTotalCorrectedCandidates;
            const int* const d_totalCorrectedSequencesBytes = d_correctedCandidatesOffsetsTmp + numTotalCorrectedCandidates;

            auto numEditsPerCandidateTransformed = thrust::make_transform_iterator(
                d_numEditsPerCorrectedCandidate.data(),
                [doNotUseEditsValue = getDoNotUseEditsValue()] __device__ (const auto& num){ 
                    return num == doNotUseEditsValue ? 0 : num;
                }
            );

            auto correctedCandidatesPitches = thrust::make_transform_iterator(
                d_numEditsPerCorrectedCandidate.data(),
                ReplaceNumberOp(getDoNotUseEditsValue(), decodedSequencePitchInBytes)
            );

            //compute edit offsets for compacted edits
            //compact corrected candidate sequences
            CubCallWrapper(mr).cubInclusiveScan(
                thrust::make_zip_iterator(thrust::make_tuple(
                   numEditsPerCandidateTransformed,
                   correctedCandidatesPitches
                )), 
                thrust::make_zip_iterator(thrust::make_tuple(
                   d_editsOffsetsTmp + 1,
                   d_correctedCandidatesOffsetsTmp + 1
                )),
                [] __device__ (const auto& l , const auto& r){
                    //tuple addition
                    return thrust::make_tuple(
                        thrust::get<0>(l) + thrust::get<0>(r),
                        thrust::get<1>(l) + thrust::get<1>(r)
                    );
                },
                numTotalCorrectedCandidates,
                stream
            );

            h_indicesForGather.resize(2);
            helpers::lambda_kernel<<<1,1,0,stream>>>([
                output = h_indicesForGather.data(), 
                d_totalNumberOfEdits, 
                d_totalCorrectedSequencesBytes
            ]__device__(){
                output[0] = *d_totalNumberOfEdits;
                output[1] = *d_totalCorrectedSequencesBytes;
            }); CUDACHECKASYNC

            CUDACHECK(cudaStreamSynchronize(stream));

            const int& h_totalNumberOfEdits = h_indicesForGather[0];
            const int& h_totalBytesOfCorrectedCandidates = h_indicesForGather[1];

            CubCallWrapper(mr).cubExclusiveSum(
                d_num_corrected_candidates_per_anchor.data(), 
                d_num_corrected_candidates_per_anchor_prefixsum.data(), 
                currentNumAnchors, 
                stream
            );

            helpers::call_copy_n_kernel(
                thrust::make_zip_iterator(thrust::make_tuple(
                    d_num_corrected_candidates_per_anchor_prefixsum.data(), 
                    d_num_corrected_candidates_per_anchor.data()
                )),
                currentNumAnchors,
                thrust::make_zip_iterator(thrust::make_tuple(
                    currentOutput->h_num_corrected_candidates_per_anchor_prefixsum.data(), 
                    currentOutput->h_num_corrected_candidates_per_anchor.data()
                )),
                stream
            );

            if(numTotalCorrectedCandidates > 0){

                {

                    rmm::device_uvector<int> d_alignment_shifts2(numTotalCorrectedCandidates, stream, mr);
                    rmm::device_uvector<read_number> d_candidate_read_ids2(numTotalCorrectedCandidates, stream, mr);
                    rmm::device_uvector<int> d_candidate_sequences_lengths2(numTotalCorrectedCandidates, stream, mr);

                    helpers::call_compact_kernel_async(
                        thrust::make_zip_iterator(thrust::make_tuple(
                            d_alignment_shifts2.data(), 
                            d_candidate_read_ids2.data(),
                            d_candidate_sequences_lengths2.data()
                        )),
                        thrust::make_zip_iterator(thrust::make_tuple(
                            d_alignment_shifts.data(), 
                            d_candidate_read_ids.data(),
                            d_candidate_sequences_lengths.data()
                        )),
                        d_indices_of_corrected_candidates.data(), 
                        numTotalCorrectedCandidates,
                        stream
                    );

                    helpers::call_copy_n_kernel(
                        thrust::make_zip_iterator(thrust::make_tuple(
                            d_alignment_shifts2.data(), 
                            d_candidate_read_ids2.data(),
                            d_candidate_sequences_lengths2.data(),
                            d_indices_of_corrected_candidates.data(),
                            d_numEditsPerCorrectedCandidate.data()
                        )), 
                        numTotalCorrectedCandidates, 
                        thrust::make_zip_iterator(thrust::make_tuple(
                            currentOutput->h_alignment_shifts.data(), 
                            currentOutput->h_candidate_read_ids.data(),
                            currentOutput->h_candidate_sequences_lengths.data(),
                            currentOutput->h_indices_of_corrected_candidates.data(),
                            currentOutput->h_numEditsPerCorrectedCandidate.data()
                        )), 
                        stream
                    );

                }

                if(h_totalNumberOfEdits > 0){
                    //compact edits
                    rmm::device_uvector<EncodedCorrectionEdit> d_editsPerCorrectedCandidate2(h_totalNumberOfEdits + sizeof(int), stream, mr);

                    gpucorrectorkernels::compactEditsKernel<<<SDIV(numTotalCorrectedCandidates, 128), 128, 0, stream>>>(
                        d_editsPerCorrectedCandidate.data(),
                        d_editsPerCorrectedCandidate2.data(),
                        d_editsOffsetsTmp,
                        d_num_total_corrected_candidates.data(),
                        d_numEditsPerCorrectedCandidate.data(),
                        getDoNotUseEditsValue(),
                        editsPitchInBytes
                    ); CUDACHECKASYNC;

                    CUDACHECK(cudaMemcpyAsync(
                        currentOutput->h_editsPerCorrectedCandidate.data(),
                        d_editsPerCorrectedCandidate2.data(),
                        sizeof(EncodedCorrectionEdit) * h_totalNumberOfEdits,
                        D2H,
                        stream
                    ));

                    CUDACHECK(cudaMemcpyAsync(
                        currentOutput->h_candidateEditOffsets.data() + 1,
                        d_editsOffsetsTmp + 1,
                        sizeof(int) * (numTotalCorrectedCandidates-1),
                        D2H,
                        stream
                    ));

                    currentOutput->h_candidateEditOffsets[0] = 0;
                }else{
                    std::fill_n(currentOutput->h_candidateEditOffsets.data(), numTotalCorrectedCandidates, 0);
                }
                
                if(h_totalBytesOfCorrectedCandidates > 0){
                    //compact sequences
                    rmm::device_uvector<char> d_corrected_candidates2(h_totalBytesOfCorrectedCandidates, stream, mr);

                    gpucorrectorkernels::compactCorrectedSequencesKernel<32><<<SDIV(numTotalCorrectedCandidates, 128), 128, 0, stream>>>(
                        d_corrected_candidates.data(),
                        d_corrected_candidates2.data(),
                        this->decodedSequencePitchInBytes,
                        d_num_total_corrected_candidates.data(),
                        d_numEditsPerCorrectedCandidate.data(),
                        getDoNotUseEditsValue(),
                        d_correctedCandidatesOffsetsTmp,
                        thrust::make_counting_iterator(0)
                    ); CUDACHECKASYNC;

                    CUDACHECK(cudaMemcpyAsync(
                        currentOutput->h_corrected_candidates.data(),
                        d_corrected_candidates2.data(),
                        h_totalBytesOfCorrectedCandidates,
                        D2H,
                        stream
                    ));

                    CUDACHECK(cudaMemcpyAsync(
                        currentOutput->h_correctedCandidatesOffsets.data() + 1,
                        d_correctedCandidatesOffsetsTmp + 1,
                        sizeof(int) * (numTotalCorrectedCandidates - 1),
                        D2H,
                        stream
                    ));

                    currentOutput->h_correctedCandidatesOffsets[0] = 0;
                }else{
                    std::fill_n(currentOutput->h_correctedCandidatesOffsets.data(), numTotalCorrectedCandidates, 0);
                }                
            }

        }

        void copyCandidateResultsFromDeviceToHostForestGpu(cudaStream_t stream){
            copyCandidateResultsFromDeviceToHostClassic(stream);
        }

        void getAmbiguousFlagsOfAnchors(cudaStream_t stream){

            gpuReadStorage->areSequencesAmbiguous(
                readstorageHandle,
                d_anchorContainsN.data(), 
                d_anchorReadIds.data(), 
                currentNumAnchors,
                stream
            );
        }

        void getAmbiguousFlagsOfCandidates(cudaStream_t stream){
            gpuReadStorage->areSequencesAmbiguous(
                readstorageHandle,
                d_candidateContainsN.data(), 
                d_candidate_read_ids.data(), 
                currentNumCandidates,
                stream
            ); 
        }

        void getCandidateAlignments(cudaStream_t stream){

            const bool removeAmbiguousAnchors = programOptions->excludeAmbiguousReads;
            const bool removeAmbiguousCandidates = programOptions->excludeAmbiguousReads;
   
            callShiftedHammingDistanceKernel(
                d_alignment_overlaps.data(),
                d_alignment_shifts.data(),
                d_alignment_nOps.data(),
                d_alignment_best_alignment_flags.data(),
                d_anchor_sequences_data.data(),
                d_candidate_sequences_data.data(),
                d_anchor_sequences_lengths.data(),
                d_candidate_sequences_lengths.data(),
                d_anchorIndicesOfCandidates.data(),
                currentNumAnchors,
                currentNumCandidates,
                d_anchorContainsN.data(),
                removeAmbiguousAnchors,
                d_candidateContainsN.data(),
                removeAmbiguousCandidates,
                gpuReadStorage->getSequenceLengthUpperBound(),
                gpuReadStorage->getSequenceLengthUpperBound(),
                encodedSequencePitchInInts,
                encodedSequencePitchInInts,
                programOptions->min_overlap,
                programOptions->maxErrorRate,
                programOptions->min_overlap_ratio,
                programOptions->estimatedErrorrate,
                stream
            );

            #if 1
            if(!gpuReadStorage->isPairedEnd()){
                //default kernel
                call_cuda_filter_alignments_by_mismatchratio_kernel_async(
                    d_alignment_best_alignment_flags.data(),
                    d_alignment_nOps.data(),
                    d_alignment_overlaps.data(),
                    d_candidates_per_anchor_prefixsum.data(),
                    currentNumAnchors,
                    currentNumCandidates,
                    programOptions->estimatedErrorrate,
                    programOptions->estimatedCoverage * programOptions->m_coverage,
                    stream
                );
            }else{
                helpers::lambda_kernel<<<SDIV(currentNumCandidates, 128), 128, 0, stream>>>(
                    [
                        bestAlignmentFlags = d_alignment_best_alignment_flags.data(),
                        nOps = d_alignment_nOps.data(),
                        overlaps = d_alignment_overlaps.data(),
                        currentNumCandidates = currentNumCandidates,
                        d_isPairedCandidate = d_isPairedCandidate.data(),
                        pairedFilterThreshold = programOptions->pairedFilterThreshold
                    ] __device__(){
                        const int tid = threadIdx.x + blockIdx.x * blockDim.x;
                        const int stride = blockDim.x * gridDim.x;

                        for(int candidate_index = tid; candidate_index < currentNumCandidates; candidate_index += stride){
                            if(!d_isPairedCandidate[candidate_index]){
                                if(bestAlignmentFlags[candidate_index] != AlignmentOrientation::None) {

                                    const int alignment_overlap = overlaps[candidate_index];
                                    const int alignment_nops = nOps[candidate_index];

                                    const float mismatchratio = float(alignment_nops) / alignment_overlap;

                                    if(mismatchratio >= pairedFilterThreshold) {
                                        bestAlignmentFlags[candidate_index] = AlignmentOrientation::None;
                                    }
                                }
                            }
                        }
                    }
                );
            }
            #else
                //default kernel
                call_cuda_filter_alignments_by_mismatchratio_kernel_async(
                    d_alignment_best_alignment_flags.data(),
                    d_alignment_nOps.data(),
                    d_alignment_overlaps.data(),
                    d_candidates_per_anchor_prefixsum.data(),
                    d_numAnchors.data(),
                    d_numCandidates.data(),
                    maxAnchors,
                    currentNumCandidates,
                    programOptions->estimatedErrorrate,
                    programOptions->estimatedCoverage * programOptions->m_coverage,
                    stream
                );
            #endif

            callSelectIndicesOfGoodCandidatesKernelAsync(
                d_indices.data(),
                d_indices_per_anchor.data(),
                d_num_indices.data(),
                d_alignment_best_alignment_flags.data(),
                d_candidates_per_anchor.data(),
                d_candidates_per_anchor_prefixsum.data(),
                d_anchorIndicesOfCandidates.data(),
                currentNumAnchors,
                currentNumCandidates,
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                h_numRemainingCandidatesAfterAlignment.data(),
                d_num_indices.data(),
                sizeof(int),
                D2H,
                stream
            ));

            CUDACHECK(cudaEventRecord(events[1], stream));

            //compactCandidatesByAlignmentFlag2(stream);
        }

        void compactCandidatesByAlignmentFlag(cudaStream_t stream){
            auto policy = rmm::exec_policy_nosync(stream, mr);

            rmm::device_uvector<int> d_num_indices_new(1, stream, mr);        
            rmm::device_uvector<int> d_numCandidates_new(1, stream, mr);            

            rmm::device_uvector<int> d_inputPositions(currentNumCandidates, stream, mr);
            auto newNumCandidates = thrust::distance(
                d_inputPositions.begin(),
                thrust::copy_if(
                    policy,
                    thrust::make_counting_iterator(0),
                    thrust::make_counting_iterator(0) + currentNumCandidates,
                    d_alignment_best_alignment_flags.begin(),
                    d_inputPositions.begin(),              
                    []__host__ __device__ (const AlignmentOrientation& o){
                        return o != AlignmentOrientation::None;
                    }
                )
            );

            CUDACHECK(cudaMemcpyAsync(
                d_numCandidates_new.data(),
                &newNumCandidates,
                sizeof(int),
                H2D,
                stream
            ));

            CUDACHECK(cudaMemcpyAsync(
                d_num_indices_new.data(),
                d_numCandidates_new.data(),
                sizeof(int),
                D2D,
                stream
            ));


            auto newNumCandidatesRounded = newNumCandidates;

            rmm::device_uvector<bool> d_candidateContainsN_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_candidate_sequences_lengths_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_alignment_overlaps_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_alignment_shifts_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_alignment_nOps_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<AlignmentOrientation> d_alignment_best_alignment_flags_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<bool> d_isPairedCandidate_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<read_number> d_candidate_read_ids_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_anchorIndicesOfCandidates_new(newNumCandidatesRounded, stream, mr);

            thrust::gather(
                policy,
                d_inputPositions.begin(),
                d_inputPositions.begin() + newNumCandidates,
                thrust::make_zip_iterator(thrust::make_tuple(
                    d_candidateContainsN.begin(),
                    d_candidate_sequences_lengths.begin(),
                    d_alignment_overlaps.begin(),
                    d_alignment_shifts.begin(),
                    d_alignment_nOps.begin(),
                    d_alignment_best_alignment_flags.begin(),
                    d_isPairedCandidate.begin(),
                    d_candidate_read_ids.begin(),
                    d_anchorIndicesOfCandidates.begin()
                )),
                thrust::make_zip_iterator(thrust::make_tuple(
                    d_candidateContainsN_new.begin(),
                    d_candidate_sequences_lengths_new.begin(),
                    d_alignment_overlaps_new.begin(),
                    d_alignment_shifts_new.begin(),
                    d_alignment_nOps_new.begin(),
                    d_alignment_best_alignment_flags_new.begin(),
                    d_isPairedCandidate_new.begin(),
                    d_candidate_read_ids_new.begin(),
                    d_anchorIndicesOfCandidates_new.begin()
                ))
            );


            rmm::device_uvector<unsigned int> d_candidate_sequences_data_new(newNumCandidatesRounded * encodedSequencePitchInInts, stream, mr);           

            helpers::lambda_kernel<<<SDIV(newNumCandidates, 128 / 8), 128, 0, stream>>>(
                [
                    d_inputPositions = d_inputPositions.data(),
                    d_candidate_sequences_data_new = d_candidate_sequences_data_new.data(),
                    d_candidate_sequences_data = d_candidate_sequences_data.data(),
                    newNumCandidates = newNumCandidates,
                    encodedSequencePitchInInts = encodedSequencePitchInInts
                ] __device__ (){
                    constexpr int groupsize = 8;
                    auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());

                    const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
                    const int numGroups = (blockDim.x * gridDim.x) / groupsize;

                    for(int c = groupId; c < newNumCandidates; c += numGroups){
                        const std::size_t srcindex = *(d_inputPositions + c);
                        const std::size_t destindex = c;

                        for(int k = group.thread_rank(); k < encodedSequencePitchInInts; k += group.size()){
                            d_candidate_sequences_data_new[destindex * encodedSequencePitchInInts + k]
                                = d_candidate_sequences_data[srcindex * encodedSequencePitchInInts + k];
                        }
                    }
                }
            ); CUDACHECKASYNC

            rmm::device_uvector<int> d_candidates_per_anchor_new(currentNumAnchors, stream, mr);
            rmm::device_uvector<int> d_candidates_per_anchor_prefixsum_new(currentNumAnchors+1, stream, mr);
            rmm::device_uvector<int> d_indices_per_anchor_new(currentNumAnchors, stream, mr);

            rmm::device_uvector<int> d_tmp1(currentNumAnchors, stream, mr);
            rmm::device_uvector<int> d_tmp2(currentNumAnchors, stream, mr);
            auto nonEmpty = thrust::distance(d_tmp1.begin(),
                thrust::reduce_by_key(
                    policy,
                    d_anchorIndicesOfCandidates_new.begin(),
                    d_anchorIndicesOfCandidates_new.begin() + newNumCandidates,
                    thrust::make_constant_iterator(1),
                    d_tmp1.begin(),
                    d_tmp2.begin()
                ).first
            );
            thrust::fill(policy, d_candidates_per_anchor_new.begin(), d_candidates_per_anchor_new.end(), 0);
            thrust::scatter(
                policy,
                d_tmp2.begin(),
                d_tmp2.begin() + nonEmpty,
                d_tmp1.begin(),
                d_candidates_per_anchor_new.begin()
            );

            thrust::inclusive_scan(
                policy,
                d_candidates_per_anchor_new.begin(),
                d_candidates_per_anchor_new.begin() + currentNumAnchors,
                d_candidates_per_anchor_prefixsum_new.begin() + 1
            );

            CUDACHECK(cudaMemsetAsync(d_candidates_per_anchor_prefixsum_new.data(), 0, sizeof(int), stream));
            CUDACHECK(cudaMemcpyAsync(
                d_indices_per_anchor_new.data(),
                d_candidates_per_anchor_new.data(),
                sizeof(int) * currentNumAnchors,
                D2D,
                stream
            ));

            rmm::device_uvector<int> d_indices_new(newNumCandidatesRounded, stream, mr);

            helpers::lambda_kernel<<<SDIV(currentNumAnchors, 128 / 32), 128, 0, stream>>>(
                [
                    d_candidates_per_anchor_new = d_candidates_per_anchor_new.data(),
                    d_candidates_per_anchor_prefixsum_new = d_candidates_per_anchor_prefixsum_new.data(),
                    d_indices_new = d_indices_new.data(),
                    currentNumAnchors = currentNumAnchors,
                    newNumCandidatesRounded = int(newNumCandidatesRounded)
                ] __device__ (){
                    constexpr int groupsize = 32;
                    auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());

                    const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
                    const int numGroups = (blockDim.x * gridDim.x) / groupsize;

                    for(int a = groupId; a < currentNumAnchors; a += numGroups){
                        const int numCands = d_candidates_per_anchor_new[a];
                        const int offset = d_candidates_per_anchor_prefixsum_new[a];

                        for(int i = group.thread_rank(); i < numCands; i += group.size()){
                            d_indices_new[offset + i] = i;
                        }
                    }
                }
            ); CUDACHECKASYNC

            #if 1

            std::swap(d_num_indices, d_num_indices_new);
            std::swap(d_numCandidates, d_numCandidates_new);
            std::swap(d_candidateContainsN, d_candidateContainsN_new);
            std::swap(d_candidate_sequences_lengths, d_candidate_sequences_lengths_new);
            std::swap(d_alignment_overlaps, d_alignment_overlaps_new);
            std::swap(d_alignment_shifts, d_alignment_shifts_new);
            std::swap(d_alignment_nOps, d_alignment_nOps_new);
            std::swap(d_alignment_best_alignment_flags, d_alignment_best_alignment_flags_new);
            std::swap(d_isPairedCandidate, d_isPairedCandidate_new);
            std::swap(d_candidate_read_ids, d_candidate_read_ids_new);
            std::swap(d_anchorIndicesOfCandidates, d_anchorIndicesOfCandidates_new);
            std::swap(d_candidate_sequences_data, d_candidate_sequences_data_new);
            std::swap(d_candidates_per_anchor, d_candidates_per_anchor_new);
            std::swap(d_candidates_per_anchor_prefixsum, d_candidates_per_anchor_prefixsum_new);
            std::swap(d_indices_per_anchor, d_indices_per_anchor_new);
            std::swap(d_indices, d_indices_new);

            currentNumCandidates = newNumCandidates;
            *h_num_indices = newNumCandidates;

            #endif
        }

        void compactCandidatesByAlignmentFlag2(cudaStream_t stream){
            auto policy = rmm::exec_policy_nosync(stream, mr);

            rmm::device_uvector<int> d_num_indices_new(1, stream, mr);        
            rmm::device_uvector<int> d_numCandidates_new(1, stream, mr);            

            rmm::device_uvector<int> d_inputPositions(currentNumCandidates, stream, mr);
            CubCallWrapper(mr).cubSelectFlagged(
                thrust::make_counting_iterator(0),
                thrust::make_transform_iterator(
                    d_alignment_best_alignment_flags.begin(),            
                    []__host__ __device__ (const AlignmentOrientation& o){
                        return o != AlignmentOrientation::None;
                    }
                ),
                d_inputPositions.begin(),
                d_numCandidates_new.data(),
                currentNumCandidates,
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                h_num_indices,
                d_numCandidates_new.data(),
                sizeof(int),
                D2H,
                stream
            ));

            CUDACHECK(cudaMemcpyAsync(
                d_num_indices_new.data(),
                d_numCandidates_new.data(),
                sizeof(int),
                D2D,
                stream
            ));

            CUDACHECK(cudaStreamSynchronize(stream));
            const int newNumCandidates = *h_num_indices;

            auto newNumCandidatesRounded = newNumCandidates;

            rmm::device_uvector<bool> d_candidateContainsN_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_candidate_sequences_lengths_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_alignment_overlaps_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_alignment_shifts_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_alignment_nOps_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<AlignmentOrientation> d_alignment_best_alignment_flags_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<bool> d_isPairedCandidate_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<read_number> d_candidate_read_ids_new(newNumCandidatesRounded, stream, mr);
            rmm::device_uvector<int> d_anchorIndicesOfCandidates_new(newNumCandidatesRounded, stream, mr);

            helpers::call_compact_kernel_async(
                thrust::make_zip_iterator(thrust::make_tuple(
                    d_candidateContainsN_new.begin(),
                    d_candidate_sequences_lengths_new.begin(),
                    d_alignment_overlaps_new.begin(),
                    d_alignment_shifts_new.begin(),
                    d_alignment_nOps_new.begin(),
                    d_alignment_best_alignment_flags_new.begin(),
                    d_isPairedCandidate_new.begin(),
                    d_candidate_read_ids_new.begin(),
                    d_anchorIndicesOfCandidates_new.begin()
                )),
                thrust::make_zip_iterator(thrust::make_tuple(
                    d_candidateContainsN.begin(),
                    d_candidate_sequences_lengths.begin(),
                    d_alignment_overlaps.begin(),
                    d_alignment_shifts.begin(),
                    d_alignment_nOps.begin(),
                    d_alignment_best_alignment_flags.begin(),
                    d_isPairedCandidate.begin(),
                    d_candidate_read_ids.begin(),
                    d_anchorIndicesOfCandidates.begin()
                )),
                d_inputPositions.begin(),
                d_numCandidates_new.data(),
                currentNumCandidates,
                stream
            );

            rmm::device_uvector<unsigned int> d_candidate_sequences_data_new(newNumCandidatesRounded * encodedSequencePitchInInts, stream, mr);           

            helpers::lambda_kernel<<<SDIV(newNumCandidates, 128 / 8), 128, 0, stream>>>(
                [
                    d_inputPositions = d_inputPositions.data(),
                    d_candidate_sequences_data_new = d_candidate_sequences_data_new.data(),
                    d_candidate_sequences_data = d_candidate_sequences_data.data(),
                    newNumCandidates = newNumCandidates,
                    encodedSequencePitchInInts = encodedSequencePitchInInts
                ] __device__ (){
                    constexpr int groupsize = 8;
                    auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());

                    const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
                    const int numGroups = (blockDim.x * gridDim.x) / groupsize;

                    for(int c = groupId; c < newNumCandidates; c += numGroups){
                        const std::size_t srcindex = *(d_inputPositions + c);
                        const std::size_t destindex = c;

                        for(int k = group.thread_rank(); k < encodedSequencePitchInInts; k += group.size()){
                            d_candidate_sequences_data_new[destindex * encodedSequencePitchInInts + k]
                                = d_candidate_sequences_data[srcindex * encodedSequencePitchInInts + k];
                        }
                    }
                }
            ); CUDACHECKASYNC

            rmm::device_uvector<int> d_candidates_per_anchor_new(currentNumAnchors, stream, mr);
            rmm::device_uvector<int> d_candidates_per_anchor_prefixsum_new(currentNumAnchors+1, stream, mr);
            rmm::device_uvector<int> d_indices_per_anchor_new(currentNumAnchors, stream, mr);

            rmm::device_uvector<int> d_tmp1(currentNumAnchors, stream, mr);
            rmm::device_uvector<int> d_tmp2(currentNumAnchors, stream, mr);
            rmm::device_scalar<int> d_numNonEmptySegments(stream, mr);

            #if 1

            CubCallWrapper(mr).cubReduceByKey(
                d_anchorIndicesOfCandidates_new.data(),
                d_tmp1.data(),
                thrust::make_constant_iterator(1),
                d_tmp2.data(),
                d_numNonEmptySegments.data(),
                thrust::plus<int>{},
                newNumCandidates,
                stream
            );
            #else
            int nonEmpty = thrust::distance(d_tmp1.begin(),
                thrust::reduce_by_key(
                    policy,
                    d_anchorIndicesOfCandidates_new.begin(),
                    d_anchorIndicesOfCandidates_new.begin() + newNumCandidates,
                    thrust::make_constant_iterator(1),
                    d_tmp1.begin(),
                    d_tmp2.begin()
                ).first
            );
            #endif

            helpers::call_fill_kernel_async(d_candidates_per_anchor_new.begin(), currentNumAnchors, 0, stream);

            #if 1
            helpers::lambda_kernel<<<SDIV(newNumCandidates, 128), 128, 0, stream>>>(
                [
                    outputpositions = d_tmp1.data(),
                    input = d_tmp2.data(),
                    output = d_candidates_per_anchor_new.data(),
                    d_numNonEmptySegments = d_numNonEmptySegments.data()
                ] __device__ (){
                    for(int i = threadIdx.x + blockIdx.x * blockDim.x; i < *d_numNonEmptySegments; i += gridDim.x * blockDim.x){
                        output[outputpositions[i]] = input[i];
                    }
                }
            ); CUDACHECKASYNC

            #else
            helpers::lambda_kernel<<<SDIV(nonEmpty, 128), 128, 0, stream>>>(
                [
                    outputpositions = d_tmp1.data(),
                    input = d_tmp2.data(),
                    output = d_candidates_per_anchor_new.data(),
                    nonEmpty
                ] __device__ (){
                    for(int i = threadIdx.x + blockIdx.x * blockDim.x; i < nonEmpty; i += gridDim.x * blockDim.x){
                        output[outputpositions[i]] = input[i];
                    }
                }
            ); CUDACHECKASYNC
            #endif

            CubCallWrapper(mr).cubInclusiveSum(
                d_candidates_per_anchor_new.begin(),
                d_candidates_per_anchor_prefixsum_new.begin() + 1,
                currentNumAnchors,
                stream
            );

            CUDACHECK(cudaMemsetAsync(d_candidates_per_anchor_prefixsum_new.data(), 0, sizeof(int), stream));
            CUDACHECK(cudaMemcpyAsync(
                d_indices_per_anchor_new.data(),
                d_candidates_per_anchor_new.data(),
                sizeof(int) * currentNumAnchors,
                D2D,
                stream
            ));

            rmm::device_uvector<int> d_indices_new(newNumCandidatesRounded, stream, mr);

            helpers::lambda_kernel<<<SDIV(currentNumAnchors, 128 / 32), 128, 0, stream>>>(
                [
                    d_candidates_per_anchor_new = d_candidates_per_anchor_new.data(),
                    d_candidates_per_anchor_prefixsum_new = d_candidates_per_anchor_prefixsum_new.data(),
                    d_indices_new = d_indices_new.data(),
                    currentNumAnchors = currentNumAnchors,
                    newNumCandidatesRounded = int(newNumCandidatesRounded)
                ] __device__ (){
                    constexpr int groupsize = 32;
                    auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());

                    const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
                    const int numGroups = (blockDim.x * gridDim.x) / groupsize;

                    for(int a = groupId; a < currentNumAnchors; a += numGroups){
                        const int numCands = d_candidates_per_anchor_new[a];
                        const int offset = d_candidates_per_anchor_prefixsum_new[a];

                        for(int i = group.thread_rank(); i < numCands; i += group.size()){
                            d_indices_new[offset + i] = i;
                        }
                    }
                }
            ); CUDACHECKASYNC

            #if 1

            std::swap(d_num_indices, d_num_indices_new);
            std::swap(d_numCandidates, d_numCandidates_new);
            std::swap(d_candidateContainsN, d_candidateContainsN_new);
            std::swap(d_candidate_sequences_lengths, d_candidate_sequences_lengths_new);
            std::swap(d_alignment_overlaps, d_alignment_overlaps_new);
            std::swap(d_alignment_shifts, d_alignment_shifts_new);
            std::swap(d_alignment_nOps, d_alignment_nOps_new);
            std::swap(d_alignment_best_alignment_flags, d_alignment_best_alignment_flags_new);
            std::swap(d_isPairedCandidate, d_isPairedCandidate_new);
            std::swap(d_candidate_read_ids, d_candidate_read_ids_new);
            std::swap(d_anchorIndicesOfCandidates, d_anchorIndicesOfCandidates_new);
            std::swap(d_candidate_sequences_data, d_candidate_sequences_data_new);
            std::swap(d_candidates_per_anchor, d_candidates_per_anchor_new);
            std::swap(d_candidates_per_anchor_prefixsum, d_candidates_per_anchor_prefixsum_new);
            std::swap(d_indices_per_anchor, d_indices_per_anchor_new);
            std::swap(d_indices, d_indices_new);

            currentNumCandidates = newNumCandidates;
            *h_num_indices = newNumCandidates;

            #endif
        }

        void buildAndRefineMultipleSequenceAlignment(cudaStream_t stream){

            char* d_anchor_qual = nullptr;
            char* d_cand_qual = nullptr;

            if(programOptions->useQualityScores){
                #if 1

                cudaStream_t qualityStream = extraStream;

                CUDACHECK(cudaStreamWaitEvent(qualityStream, inputCandidateDataIsReadyEvent, 0));

                d_anchor_qual = reinterpret_cast<char*>(mr->allocate(currentNumAnchors * qualityPitchInBytes, qualityStream));

                gpuReadStorage->gatherContiguousQualities(
                    readstorageHandle,
                    d_anchor_qual,
                    qualityPitchInBytes,
                    currentInput->h_anchorReadIds[0],
                    currentNumAnchors,
                    qualityStream,
                    mr
                );

                d_cand_qual = reinterpret_cast<char*>(mr->allocate(currentNumCandidates * qualityPitchInBytes, qualityStream));

                gpuReadStorage->gatherQualities(
                    readstorageHandle,
                    d_cand_qual,
                    qualityPitchInBytes,
                    makeAsyncConstBufferWrapper(currentInput->h_candidate_read_ids.data()),
                    d_candidate_read_ids.data(),
                    currentNumCandidates,
                    qualityStream,
                    mr
                );
                
                CUDACHECK(cudaEventRecord(events[0], qualityStream));
                CUDACHECK(cudaStreamWaitEvent(stream, events[0], 0));

                #else 

                rmm::device_uvector<int> d_prefixsum(maxAnchors + 1, stream, mr);

                CubCallWrapper(mr).cubExclusiveSum(
                    d_indices_per_anchor.data(),
                    d_prefixsum.data(),
                    maxAnchors,
                    stream
                );

                auto zippedValid = thrust::make_zip_iterator(thrust::make_tuple(
                    h_indicesForGather.data(),
                    d_indicesForGather.data()
                ));

                auto duplicateInput = [] __host__ __device__ (read_number id){
                    return thrust::make_tuple(id, id);
                };
                auto duplicatedIds = thrust::make_transform_iterator(
                    d_candidate_read_ids.data(),
                    duplicateInput
                );
                auto copyifend = thrust::copy_if(
                    rmm::exec_policy_nosync(stream, mr),
                    duplicatedIds,
                    duplicatedIds + currentNumCandidates,
                    d_alignment_best_alignment_flags.data(),
                    zippedValid,
                    [] __host__ __device__ (const AlignmentOrientation& o){
                        return o != AlignmentOrientation::None;
                    }
                );

                const int hNumIndices = thrust::distance(zippedValid, copyifend);

                d_anchor_qual = reinterpret_cast<char*>(mr->allocate(currentNumAnchors * qualityPitchInBytes, stream));

                gpuReadStorage->gatherContiguousQualities(
                    readstorageHandle,
                    d_anchor_qual,
                    qualityPitchInBytes,
                    currentInput->h_anchorReadIds[0],
                    currentNumAnchors,
                    qualityStream,
                    mr
                );               

                rmm::device_uvector<char> d_candidate_qualities_compact(hNumIndices * qualityPitchInBytes, stream, mr);

                nvtx::push_range("get compact qscores " + std::to_string(hNumIndices) + " " + std::to_string(currentNumCandidates), 6);
                gpuReadStorage->gatherQualities(
                    readstorageHandle,
                    d_candidate_qualities_compact.data(),
                    qualityPitchInBytes,
                    makeAsyncConstBufferWrapper(h_indicesForGather.data(), events[1]),
                    d_indicesForGather.data(),
                    hNumIndices,
                    stream,
                    mr
                );
                nvtx::pop_range();

                d_cand_qual = reinterpret_cast<char*>(mr->allocate(currentNumCandidates * qualityPitchInBytes, stream));

                #if 0
                //scatter compact quality scores to correct positions
                helpers::lambda_kernel<<<SDIV(hNumIndices, 256 / 8), 256, 0, stream>>>(
                    [
                        d_candidate_qualities_compact = d_candidate_qualities_compact.data(),
                        d_candidate_qualities = d_cand_qual,
                        d_candidate_sequences_lengths = d_candidate_sequences_lengths.data(),
                        qualityPitchInBytes = qualityPitchInBytes,
                        d_indices = d_indices.data(),
                        d_indices_per_anchor = d_indices_per_anchor.data(),
                        d_indices_per_anchor_prefixsum = d_prefixsum.data(),
                        d_num_indices = d_num_indices.data(),
                        d_candidates_per_anchor_prefixsum = d_candidates_per_anchor_prefixsum.data(),
                        currentNumAnchors = currentNumAnchors,
                        hNumIndices = hNumIndices
                    ] __device__ (){
                        constexpr int groupsize = 8;
                        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());

                        const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
                        const int numGroups = (blockDim.x * gridDim.x) / groupsize;

                        assert(qualityPitchInBytes % sizeof(int) == 0);

                        for(int c = groupId; c < hNumIndices; c += numGroups){
                            const int anchorIndex = thrust::distance(
                                d_indices_per_anchor_prefixsum,
                                thrust::lower_bound(
                                    thrust::seq,
                                    d_indices_per_anchor_prefixsum,
                                    d_indices_per_anchor_prefixsum + currentNumAnchors,
                                    c + 1
                                )
                            )-1;

                            const int segmentOffset = d_candidates_per_anchor_prefixsum[anchorIndex];
                            const int* const myIndices = d_indices + segmentOffset;
                            const int localCandidatePositionInAnchor = groupId - d_indices_per_anchor_prefixsum[anchorIndex];
                            const int outputCandidateIndex = segmentOffset + myIndices[localCandidatePositionInAnchor];

                            const int candidateLength = d_candidate_sequences_lengths[outputCandidateIndex];
                            const int iters = SDIV(candidateLength, sizeof(int));

                            const int* const input = (const int*)(d_candidate_qualities_compact + size_t(c) * qualityPitchInBytes);
                            int* const output = (int*)(d_candidate_qualities + size_t(outputCandidateIndex) * qualityPitchInBytes);

                            for(int k = group.thread_rank(); k < iters; k += group.size()){
                                output[k] = input[k];
                            }
                        }
                    }
                ); CUDACHECKASYNC;

                #else

                //scatter compact quality scores to correct positions
                helpers::lambda_kernel<<<maxAnchors, 256, 0, stream>>>(
                    [
                        d_candidate_qualities_compact = d_candidate_qualities_compact.data(),
                        d_candidate_qualities = d_cand_qual,
                        d_candidate_sequences_lengths = d_candidate_sequences_lengths.data(),
                        qualityPitchInBytes = qualityPitchInBytes,
                        d_indices = d_indices.data(),
                        d_indices_per_anchor = d_indices_per_anchor.data(),
                        d_indices_per_anchor_prefixsum = d_prefixsum.data(),
                        d_num_indices = d_num_indices.data(),
                        d_candidates_per_anchor_prefixsum = d_candidates_per_anchor_prefixsum.data(),
                        currentNumAnchors = currentNumAnchors
                    ] __device__ (){
                        constexpr int groupsize = 8;
                        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());

                        const int groupId = threadIdx.x / groupsize;
                        const int numgroups = blockDim.x / groupsize;

                        assert(qualityPitchInBytes % sizeof(int) == 0);

                        for(int anchor = blockIdx.x; anchor < currentNumAnchors; anchor += gridDim.x){

                            const int globalCandidateOffset = d_candidates_per_anchor_prefixsum[anchor];
                            const int* const myIndices = d_indices + globalCandidateOffset;
                            const int numIndices = d_indices_per_anchor[anchor];
                            const int offset = d_indices_per_anchor_prefixsum[anchor];

                            for(int c = groupId; c < numIndices; c += numgroups){
                                const int outputpos = globalCandidateOffset + myIndices[c];
                                const int inputpos = offset + c;
                                const int length = d_candidate_sequences_lengths[outputpos];

                                const int iters = SDIV(length, sizeof(int));

                                const int* const input = (const int*)(d_candidate_qualities_compact + size_t(inputpos) * qualityPitchInBytes);
                                int* const output = (int*)(d_candidate_qualities + size_t(outputpos) * qualityPitchInBytes);

                                for(int k = group.thread_rank(); k < iters; k += group.size()){
                                    output[k] = input[k];
                                }
                            }
                        }
                    }
                ); CUDACHECKASYNC;
                #endif

                #endif

            }

            managedgpumsa = std::make_unique<ManagedGPUMultiMSA>(stream, mr, h_managedmsa_tmp.data());

            #if 0

            if(useMsaRefinement()){
                rmm::device_uvector<int> d_indices_tmp(currentNumCandidates+1, stream, mr);
                rmm::device_uvector<int> d_indices_per_anchor_tmp(maxAnchors+1, stream, mr);
                rmm::device_uvector<int> d_num_indices_tmp(1, stream, mr);
     
                managedgpumsa->constructAndRefine(
                    d_indices_tmp.data(),
                    d_indices_per_anchor_tmp.data(),
                    d_num_indices_tmp.data(),
                    d_alignment_overlaps.data(),
                    d_alignment_shifts.data(),
                    d_alignment_nOps.data(),
                    d_alignment_best_alignment_flags.data(),
                    d_anchor_sequences_lengths.data(),
                    d_candidate_sequences_lengths.data(),
                    d_indices.data(),
                    d_indices_per_anchor.data(),
                    d_candidates_per_anchor_prefixsum.data(),
                    d_anchor_sequences_data.data(),
                    d_candidate_sequences_data.data(),
                    d_isPairedCandidate.data(),
                    d_anchor_qual,
                    d_cand_qual,
                    currentNumAnchors,
                    currentNumCandidates,
                    programOptions->maxErrorRate,
                    programOptions->useQualityScores,
                    encodedSequencePitchInInts,
                    qualityPitchInBytes,
                    programOptions->estimatedCoverage,
                    getNumRefinementIterations(),
                    MSAColumnCount{static_cast<int>(msaColumnPitchInElements)},
                    stream
                );
                std::swap(d_indices_tmp, d_indices);
                std::swap(d_indices_per_anchor_tmp, d_indices_per_anchor);
                std::swap(d_num_indices_tmp, d_num_indices);
            }else{
                managedgpumsa->construct(
                    d_alignment_overlaps.data(),
                    d_alignment_shifts.data(),
                    d_alignment_nOps.data(),
                    d_alignment_best_alignment_flags.data(),
                    d_indices.data(),
                    d_indices_per_anchor.data(),
                    d_candidates_per_anchor_prefixsum.data(),
                    d_anchor_sequences_lengths.data(),
                    d_anchor_sequences_data.data(),
                    d_anchor_qual,
                    currentNumAnchors,
                    d_candidate_sequences_lengths.data(),
                    d_candidate_sequences_data.data(),
                    d_cand_qual,
                    d_isPairedCandidate.data(),
                    encodedSequencePitchInInts,
                    qualityPitchInBytes,
                    programOptions->useQualityScores,
                    programOptions->maxErrorRate,
                    MSAColumnCount{static_cast<int>(msaColumnPitchInElements)},
                    stream
                );
            }

            #else

            managedgpumsa->construct(
                d_alignment_overlaps.data(),
                d_alignment_shifts.data(),
                d_alignment_nOps.data(),
                d_alignment_best_alignment_flags.data(),
                d_indices.data(),
                d_indices_per_anchor.data(),
                d_candidates_per_anchor_prefixsum.data(),
                d_anchor_sequences_lengths.data(),
                d_anchor_sequences_data.data(),
                d_anchor_qual,
                currentNumAnchors,
                d_candidate_sequences_lengths.data(),
                d_candidate_sequences_data.data(),
                d_cand_qual,
                d_isPairedCandidate.data(),
                encodedSequencePitchInInts,
                qualityPitchInBytes,
                programOptions->useQualityScores,
                programOptions->maxErrorRate,
                MSAColumnCount{static_cast<int>(msaColumnPitchInElements)},
                stream
            );

            if(useMsaRefinement()){
                
                rmm::device_uvector<int> d_indices_tmp(currentNumCandidates+1, stream, mr);
                rmm::device_uvector<int> d_indices_per_anchor_tmp(maxAnchors+1, stream, mr);
                rmm::device_uvector<int> d_num_indices_tmp(1, stream, mr);

                managedgpumsa->refine(
                    d_indices_tmp.data(),
                    d_indices_per_anchor_tmp.data(),
                    d_num_indices_tmp.data(),
                    d_alignment_overlaps.data(),
                    d_alignment_shifts.data(),
                    d_alignment_nOps.data(),
                    d_alignment_best_alignment_flags.data(),
                    d_indices.data(),
                    d_indices_per_anchor.data(),
                    d_candidates_per_anchor_prefixsum.data(),
                    d_anchor_sequences_lengths.data(),
                    d_anchor_sequences_data.data(),
                    d_anchor_qual,
                    currentNumAnchors,
                    d_candidate_sequences_lengths.data(),
                    d_candidate_sequences_data.data(),
                    d_cand_qual,
                    d_isPairedCandidate.data(),
                    currentNumCandidates,
                    encodedSequencePitchInInts,
                    qualityPitchInBytes,
                    programOptions->useQualityScores,
                    programOptions->maxErrorRate,
                    programOptions->estimatedCoverage,
                    getNumRefinementIterations(),
                    stream
                );

                std::swap(d_indices_tmp, d_indices);
                std::swap(d_indices_per_anchor_tmp, d_indices_per_anchor);
                std::swap(d_num_indices_tmp, d_num_indices);

            }

            #endif

            if(programOptions->useQualityScores){
                mr->deallocate(d_anchor_qual, currentNumAnchors * qualityPitchInBytes, stream);
                mr->deallocate(d_cand_qual, currentNumCandidates * qualityPitchInBytes, stream);
            }
        }


        void correctAnchors(cudaStream_t stream){
            if(programOptions->correctionType == CorrectionType::Classic){
                correctAnchorsClassic(stream);
            }else if(programOptions->correctionType == CorrectionType::Forest){
                correctAnchorsForestGpu(stream);
            }else{
                throw std::runtime_error("correctAnchors not implemented for this correctionType");
            }
        }

        void correctAnchorsClassic(cudaStream_t stream){

            const float avg_support_threshold = 1.0f - 1.0f * programOptions->estimatedErrorrate;
            const float min_support_threshold = 1.0f - 3.0f * programOptions->estimatedErrorrate;
            // coverage is always >= 1
            const float min_coverage_threshold = std::max(1.0f,
                programOptions->m_coverage / 6.0f * programOptions->estimatedCoverage);
            //const float max_coverage_threshold = 0.5 * programOptions->estimatedCoverage;

            // correct anchors

            call_msaCorrectAnchorsKernel_async(
                d_corrected_anchors.data(),
                d_anchor_is_corrected.data(),
                d_is_high_quality_anchor.data(),
                managedgpumsa->multiMSAView(),
                d_anchor_sequences_data.data(),
                d_indices_per_anchor.data(),
                d_numAnchors.data(),
                maxAnchors,
                encodedSequencePitchInInts,
                decodedSequencePitchInBytes,
                programOptions->estimatedErrorrate,
                avg_support_threshold,
                min_support_threshold,
                min_coverage_threshold,
                gpuReadStorage->getSequenceLengthUpperBound(),
                stream
            );

            gpucorrectorkernels::selectIndicesOfFlagsOneBlock<256><<<1,256,0, stream>>>(
                d_indices_of_corrected_anchors.data(),
                d_num_indices_of_corrected_anchors.data(),
                d_anchor_is_corrected.data(),
                d_numAnchors.data()
            ); CUDACHECKASYNC;

            helpers::call_fill_kernel_async(d_numEditsPerCorrectedanchor.data(), currentNumAnchors, 0, stream);

            callConstructSequenceCorrectionResultsKernel(
                d_editsPerCorrectedanchor.data(),
                d_numEditsPerCorrectedanchor.data(),
                getDoNotUseEditsValue(),
                d_indices_of_corrected_anchors.data(),
                d_num_indices_of_corrected_anchors.data(),
                d_anchorContainsN.data(),
                d_anchor_sequences_data.data(),
                d_anchor_sequences_lengths.data(),
                d_corrected_anchors.data(),
                currentNumAnchors,
                false,
                maxNumEditsPerSequence,
                encodedSequencePitchInInts,
                decodedSequencePitchInBytes,
                editsPitchInBytes,      
                stream
            );
            
        }

        void correctAnchorsForestGpu(cudaStream_t stream){

            const float avg_support_threshold = 1.0f - 1.0f * programOptions->estimatedErrorrate;
            const float min_support_threshold = 1.0f - 3.0f * programOptions->estimatedErrorrate;
            // coverage is always >= 1
            const float min_coverage_threshold = std::max(1.0f,
                programOptions->m_coverage / 6.0f * programOptions->estimatedCoverage);
            const float max_coverage_threshold = 0.5 * programOptions->estimatedCoverage;

            // correct anchors

            callMsaCorrectAnchorsWithForestKernel(
                d_corrected_anchors.data(),
                d_anchor_is_corrected.data(),
                d_is_high_quality_anchor.data(),
                managedgpumsa->multiMSAView(),
                *gpuForestAnchor,
                programOptions->thresholdAnchor,
                d_anchor_sequences_data.data(),
                d_indices_per_anchor.data(),
                currentNumAnchors,
                encodedSequencePitchInInts,
                decodedSequencePitchInBytes,
                gpuReadStorage->getSequenceLengthUpperBound(),
                programOptions->estimatedErrorrate,
                programOptions->estimatedCoverage,
                avg_support_threshold,
                min_support_threshold,
                min_coverage_threshold,
                stream
            );

            gpucorrectorkernels::selectIndicesOfFlagsOneBlock<256><<<1,256,0, stream>>>(
                d_indices_of_corrected_anchors.data(),
                d_num_indices_of_corrected_anchors.data(),
                d_anchor_is_corrected.data(),
                d_numAnchors.data()
            ); CUDACHECKASYNC;

            helpers::call_fill_kernel_async(d_numEditsPerCorrectedanchor.data(), currentNumAnchors, 0, stream);

            callConstructSequenceCorrectionResultsKernel(
                d_editsPerCorrectedanchor.data(),
                d_numEditsPerCorrectedanchor.data(),
                getDoNotUseEditsValue(),
                d_indices_of_corrected_anchors.data(),
                d_num_indices_of_corrected_anchors.data(),
                d_anchorContainsN.data(),
                d_anchor_sequences_data.data(),
                d_anchor_sequences_lengths.data(),
                d_corrected_anchors.data(),
                currentNumAnchors,
                false,
                maxNumEditsPerSequence,
                encodedSequencePitchInInts,
                decodedSequencePitchInBytes,
                editsPitchInBytes,      
                stream
            );

        }
        
        void correctCandidates(cudaStream_t stream){
            if(programOptions->correctionTypeCands == CorrectionType::Classic){
                correctCandidatesClassic(stream);
            }else if(programOptions->correctionTypeCands == CorrectionType::Forest){
                correctCandidatesForestGpu(stream);
            }else{
                throw std::runtime_error("correctCandidates not implemented for this correctionTypeCands");
            }
        }

        void correctCandidatesClassic(cudaStream_t stream){

            const float min_support_threshold = 1.0f-3.0f*programOptions->estimatedErrorrate;
            // coverage is always >= 1
            const float min_coverage_threshold = std::max(1.0f,
                programOptions->m_coverage / 6.0f * programOptions->estimatedCoverage);
            const int new_columns_to_correct = programOptions->new_columns_to_correct;

            rmm::device_uvector<bool> d_candidateCanBeCorrected(currentNumCandidates, stream, mr);

            cub::TransformInputIterator<bool, IsHqAnchor, AnchorHighQualityFlag*>
                d_isHqanchor(d_is_high_quality_anchor.data(), IsHqAnchor{});

            gpucorrectorkernels::selectIndicesOfFlagsOneBlock<256><<<1,256,0, stream>>>(
                d_high_quality_anchor_indices.data(),
                d_num_high_quality_anchor_indices.data(),
                d_isHqanchor,
                d_numAnchors.data()
            ); CUDACHECKASYNC;

            gpucorrectorkernels::initArraysBeforeCandidateCorrectionKernel<<<SDIV(currentNumCandidates, 128), 128, 0, stream>>>(
                currentNumCandidates,
                d_numAnchors.data(),
                d_num_corrected_candidates_per_anchor.data(),
                d_candidateCanBeCorrected.data()
            ); CUDACHECKASYNC;

            #if 0
                rmm::device_uvector<bool> d_flagsCandidates(currentNumCandidates, stream, mr);
                bool* d_excludeFlags = d_flagsCandidates.data();
                bool* h_excludeFlags = h_flagsCandidates.data();

                //corrections of candidates for which a high quality anchor correction exists will not be used
                //-> don't compute them
                for(int i = 0; i < currentNumCandidates; i++){
                    const read_number candidateReadId = currentInput->h_candidate_read_ids[i];
                    h_excludeFlags[i] = correctionFlags->isCorrectedAsHQAnchor(candidateReadId);
                }

                helpers::call_copy_n_kernel(
                    (const int*)h_excludeFlags,
                    SDIV(currentNumCandidates, sizeof(int)),
                    (int*)d_excludeFlags,
                    stream
                );

                callFlagCandidatesToBeCorrectedWithExcludeFlagsKernel(
                    d_candidateCanBeCorrected,
                    d_num_corrected_candidates_per_anchor.data(),
                    managedgpumsa->multiMSAView(),
                    d_excludeFlags,
                    d_alignment_shifts.data(),
                    d_candidate_sequences_lengths.data(),
                    d_anchorIndicesOfCandidates.data(),
                    d_is_high_quality_anchor.data(),
                    d_candidates_per_anchor_prefixsum,
                    d_indices.data(),
                    d_indices_per_anchor.data(),
                    d_numAnchors,
                    d_numCandidates,
                    min_support_threshold,
                    min_coverage_threshold,
                    new_columns_to_correct,
                    stream
                );
            #else
            callFlagCandidatesToBeCorrectedKernel_async(
                d_candidateCanBeCorrected.data(),
                d_num_corrected_candidates_per_anchor.data(),
                managedgpumsa->multiMSAView(),
                d_alignment_shifts.data(),
                d_candidate_sequences_lengths.data(),
                d_anchorIndicesOfCandidates.data(),
                d_is_high_quality_anchor.data(),
                d_candidates_per_anchor_prefixsum.data(),
                d_indices.data(),
                d_indices_per_anchor.data(),
                d_numAnchors.data(),
                d_numCandidates.data(),
                min_support_threshold,
                min_coverage_threshold,
                new_columns_to_correct,
                stream
            );
            #endif

            CubCallWrapper(mr).cubSelectFlagged(
                cub::CountingInputIterator<int>(0),
                d_candidateCanBeCorrected.data(),
                d_indices_of_corrected_candidates.data(),
                d_num_total_corrected_candidates.data(),
                currentNumCandidates,
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                h_num_total_corrected_candidates.data(),
                d_num_total_corrected_candidates.data(),
                sizeof(int),
                D2H,
                stream
            ));
            CUDACHECK(cudaStreamSynchronize(stream));

            if((*h_num_total_corrected_candidates) > 0){

                d_corrected_candidates.resize(decodedSequencePitchInBytes * (*h_num_total_corrected_candidates), stream);
                d_numEditsPerCorrectedCandidate.resize((*h_num_total_corrected_candidates), stream);
                std::size_t numEditsCandidates = SDIV(editsPitchInBytes * (*h_num_total_corrected_candidates), sizeof(EncodedCorrectionEdit));
                d_editsPerCorrectedCandidate.resize(numEditsCandidates, stream);

                currentOutput->h_candidate_sequences_lengths.resize((*h_num_total_corrected_candidates));
                currentOutput->h_corrected_candidates.resize((*h_num_total_corrected_candidates) * decodedSequencePitchInBytes);
                currentOutput->h_editsPerCorrectedCandidate.resize(numEditsCandidates);
                currentOutput->h_numEditsPerCorrectedCandidate.resize((*h_num_total_corrected_candidates));
                currentOutput->h_candidateEditOffsets.resize((*h_num_total_corrected_candidates));
                currentOutput->h_indices_of_corrected_candidates.resize((*h_num_total_corrected_candidates));
                currentOutput->h_alignment_shifts.resize((*h_num_total_corrected_candidates));
                currentOutput->h_candidate_read_ids.resize((*h_num_total_corrected_candidates));
                currentOutput->h_correctedCandidatesOffsets.resize((*h_num_total_corrected_candidates) * decodedSequencePitchInBytes);
                
                callCorrectCandidatesKernel(
                    d_corrected_candidates.data(),            
                    managedgpumsa->multiMSAView(),
                    d_alignment_shifts.data(),
                    d_alignment_best_alignment_flags.data(),
                    d_candidate_sequences_data.data(),
                    d_candidate_sequences_lengths.data(),
                    d_candidateContainsN.data(),
                    d_indices_of_corrected_candidates.data(),
                    d_num_total_corrected_candidates.data(),
                    d_anchorIndicesOfCandidates.data(),
                    *h_num_total_corrected_candidates,
                    encodedSequencePitchInInts,
                    decodedSequencePitchInBytes,
                    gpuReadStorage->getSequenceLengthUpperBound(),
                    stream
                );            

                callConstructSequenceCorrectionResultsKernel(
                    d_editsPerCorrectedCandidate.data(),
                    d_numEditsPerCorrectedCandidate.data(),
                    getDoNotUseEditsValue(),
                    d_indices_of_corrected_candidates.data(),
                    d_num_total_corrected_candidates.data(),
                    d_candidateContainsN.data(),
                    d_candidate_sequences_data.data(),
                    d_candidate_sequences_lengths.data(),
                    d_corrected_candidates.data(),
                    (*h_num_total_corrected_candidates),
                    true,
                    maxNumEditsPerSequence,
                    encodedSequencePitchInInts,
                    decodedSequencePitchInBytes,
                    editsPitchInBytes,      
                    stream
                );
            }  
        }

        void correctCandidatesForestGpu(cudaStream_t stream){

            const float min_support_threshold = 1.0f-3.0f*programOptions->estimatedErrorrate;
            // coverage is always >= 1
            const float min_coverage_threshold = std::max(1.0f,
                programOptions->m_coverage / 6.0f * programOptions->estimatedCoverage);
            const int new_columns_to_correct = programOptions->new_columns_to_correct;

            rmm::device_uvector<bool> d_candidateCanBeCorrected(currentNumCandidates, stream, mr);

            cub::TransformInputIterator<bool, IsHqAnchor, AnchorHighQualityFlag*>
                d_isHqanchor(d_is_high_quality_anchor.data(), IsHqAnchor{});

            gpucorrectorkernels::selectIndicesOfFlagsOneBlock<256><<<1,256,0, stream>>>(
                d_high_quality_anchor_indices.data(),
                d_num_high_quality_anchor_indices.data(),
                d_isHqanchor,
                d_numAnchors.data()
            ); CUDACHECKASYNC;

            gpucorrectorkernels::initArraysBeforeCandidateCorrectionKernel<<<SDIV(currentNumCandidates, 128), 128, 0, stream>>>(
                currentNumCandidates,
                d_numAnchors.data(),
                d_num_corrected_candidates_per_anchor.data(),
                d_candidateCanBeCorrected.data()
            ); CUDACHECKASYNC;

   
            #if 1
                rmm::device_uvector<bool> d_flagsCandidates(currentNumCandidates, stream, mr);
                bool* d_excludeFlags = d_flagsCandidates.data();
                bool* h_excludeFlags = h_flagsCandidates.data();

                //corrections of candidates for which a high quality anchor correction exists will not be used
                //-> don't compute them
                for(int i = 0; i < currentNumCandidates; i++){
                    const read_number candidateReadId = currentInput->h_candidate_read_ids[i];
                    h_excludeFlags[i] = correctionFlags->isCorrectedAsHQAnchor(candidateReadId);
                }

                cudaMemcpyAsync(
                    d_excludeFlags,
                    h_excludeFlags,
                    sizeof(bool) * currentNumCandidates,
                    H2D,
                    stream
                );

                callFlagCandidatesToBeCorrectedWithExcludeFlagsKernel(
                    d_candidateCanBeCorrected.data(),
                    d_num_corrected_candidates_per_anchor.data(),
                    managedgpumsa->multiMSAView(),
                    d_excludeFlags,
                    d_alignment_shifts.data(),
                    d_candidate_sequences_lengths.data(),
                    d_anchorIndicesOfCandidates.data(),
                    d_is_high_quality_anchor.data(),
                    d_candidates_per_anchor_prefixsum.data(),
                    d_indices.data(),
                    d_indices_per_anchor.data(),
                    d_numAnchors.data(),
                    d_numCandidates.data(),
                    min_support_threshold,
                    min_coverage_threshold,
                    new_columns_to_correct,
                    stream
                );
            #else
            callFlagCandidatesToBeCorrectedKernel_async(
                d_candidateCanBeCorrected.data(),
                d_num_corrected_candidates_per_anchor.data(),
                managedgpumsa->multiMSAView(),
                d_alignment_shifts.data(),
                d_candidate_sequences_lengths.data(),
                d_anchorIndicesOfCandidates.data(),
                d_is_high_quality_anchor.data(),
                d_candidates_per_anchor_prefixsum,
                d_indices.data(),
                d_indices_per_anchor.data(),
                d_numAnchors,
                d_numCandidates,
                min_support_threshold,
                min_coverage_threshold,
                new_columns_to_correct,
                stream
            );
            #endif

            CubCallWrapper(mr).cubSelectFlagged(
                cub::CountingInputIterator<int>(0),
                d_candidateCanBeCorrected.data(),
                d_indices_of_corrected_candidates.data(),
                d_num_total_corrected_candidates.data(),
                currentNumCandidates,
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                h_num_total_corrected_candidates.data(),
                d_num_total_corrected_candidates.data(),
                sizeof(int),
                D2H,
                stream
            ));
            CUDACHECK(cudaStreamSynchronize(stream));

            d_corrected_candidates.resize(decodedSequencePitchInBytes * (*h_num_total_corrected_candidates), stream);
            d_numEditsPerCorrectedCandidate.resize((*h_num_total_corrected_candidates), stream);
            std::size_t numEditsCandidates = SDIV(editsPitchInBytes * (*h_num_total_corrected_candidates), sizeof(EncodedCorrectionEdit));
            d_editsPerCorrectedCandidate.resize(numEditsCandidates, stream);

            currentOutput->h_candidate_sequences_lengths.resize((*h_num_total_corrected_candidates));
            currentOutput->h_corrected_candidates.resize((*h_num_total_corrected_candidates) * decodedSequencePitchInBytes);
            currentOutput->h_editsPerCorrectedCandidate.resize(numEditsCandidates);
            currentOutput->h_numEditsPerCorrectedCandidate.resize((*h_num_total_corrected_candidates));
            currentOutput->h_candidateEditOffsets.resize((*h_num_total_corrected_candidates));
            currentOutput->h_indices_of_corrected_candidates.resize((*h_num_total_corrected_candidates));
            currentOutput->h_alignment_shifts.resize((*h_num_total_corrected_candidates));
            currentOutput->h_candidate_read_ids.resize((*h_num_total_corrected_candidates));
            currentOutput->h_correctedCandidatesOffsets.resize((*h_num_total_corrected_candidates) * decodedSequencePitchInBytes);
            
            callMsaCorrectCandidatesWithForestKernel(
                d_corrected_candidates.data(),            
                managedgpumsa->multiMSAView(),
                *gpuForestCandidate,
                programOptions->thresholdCands,
                programOptions->estimatedCoverage,
                d_alignment_shifts.data(),
                d_alignment_best_alignment_flags.data(),
                d_candidate_sequences_data.data(),
                d_candidate_sequences_lengths.data(),
                d_indices_of_corrected_candidates.data(),
                d_anchorIndicesOfCandidates.data(),
                *h_num_total_corrected_candidates,
                encodedSequencePitchInInts,
                decodedSequencePitchInBytes,
                gpuReadStorage->getSequenceLengthUpperBound(),
                stream
            );  

            callConstructSequenceCorrectionResultsKernel(
                d_editsPerCorrectedCandidate.data(),
                d_numEditsPerCorrectedCandidate.data(),
                getDoNotUseEditsValue(),
                d_indices_of_corrected_candidates.data(),
                d_num_total_corrected_candidates.data(),
                d_candidateContainsN.data(),
                d_candidate_sequences_data.data(),
                d_candidate_sequences_lengths.data(),
                d_corrected_candidates.data(),
                *h_num_total_corrected_candidates,
                true,
                maxNumEditsPerSequence,
                encodedSequencePitchInInts,
                decodedSequencePitchInBytes,
                editsPitchInBytes,      
                stream
            );
            
        }

        static constexpr int getDoNotUseEditsValue() noexcept{
            return -1;
        }

    private:

        int deviceId;
        std::array<CudaEvent, 2> events;
        cudaStream_t extraStream;

        CudaEvent previousBatchFinishedEvent;
        CudaEvent inputCandidateDataIsReadyEvent;

        std::size_t msaColumnPitchInElements;
        std::size_t encodedSequencePitchInInts;
        std::size_t decodedSequencePitchInBytes;
        std::size_t qualityPitchInBytes;
        std::size_t editsPitchInBytes;

        int maxAnchors;
        int maxNumEditsPerSequence;
        int currentNumAnchors;
        int currentNumCandidates;

        std::map<int, int> numCandidatesPerReadMap{};

        const ReadCorrectionFlags* correctionFlags;

        const GpuReadStorage* gpuReadStorage;

        const ProgramOptions* programOptions;

        GpuErrorCorrectorInput* currentInput;
        GpuErrorCorrectorRawOutput* currentOutput;
        GpuErrorCorrectorRawOutput currentOutputData;

        rmm::mr::device_memory_resource* mr;
        ThreadPool* threadPool;
        ThreadPool::ParallelForHandle pforHandle;

        const GpuForest* gpuForestAnchor{};
        const GpuForest* gpuForestCandidate{};

        ReadStorageHandle readstorageHandle;

        PinnedBuffer<int> h_num_total_corrected_candidates;
        PinnedBuffer<int> h_num_indices;
        PinnedBuffer<int> h_numSelected;
        PinnedBuffer<int> h_numRemainingCandidatesAfterAlignment;
        PinnedBuffer<int> h_managedmsa_tmp;

        PinnedBuffer<read_number> h_indicesForGather;
        rmm::device_uvector<read_number> d_indicesForGather;

        rmm::device_uvector<bool> d_anchorContainsN;
        rmm::device_uvector<bool> d_candidateContainsN;
        rmm::device_uvector<int> d_candidate_sequences_lengths;
        rmm::device_uvector<unsigned int> d_candidate_sequences_data;
        rmm::device_uvector<int> d_anchorIndicesOfCandidates;
        rmm::device_uvector<int> d_alignment_overlaps;
        rmm::device_uvector<int> d_alignment_shifts;
        rmm::device_uvector<int> d_alignment_nOps;
        rmm::device_uvector<AlignmentOrientation> d_alignment_best_alignment_flags; 
        rmm::device_uvector<int> d_indices;
        rmm::device_uvector<int> d_indices_per_anchor;
        rmm::device_uvector<int> d_indices_per_anchor_prefixsum;
        rmm::device_uvector<int> d_num_indices;
        rmm::device_uvector<char> d_corrected_anchors;
        rmm::device_uvector<char> d_corrected_candidates;
        rmm::device_uvector<int> d_num_corrected_candidates_per_anchor;
        rmm::device_uvector<int> d_num_corrected_candidates_per_anchor_prefixsum;
        rmm::device_uvector<int> d_num_total_corrected_candidates;
        rmm::device_uvector<bool> d_anchor_is_corrected;
        rmm::device_uvector<AnchorHighQualityFlag> d_is_high_quality_anchor;
        rmm::device_uvector<int> d_high_quality_anchor_indices;
        rmm::device_uvector<int> d_num_high_quality_anchor_indices; 
        rmm::device_uvector<EncodedCorrectionEdit> d_editsPerCorrectedanchor;
        rmm::device_uvector<int> d_numEditsPerCorrectedanchor;
        rmm::device_uvector<EncodedCorrectionEdit> d_editsPerCorrectedCandidate;
        
        rmm::device_uvector<int> d_numEditsPerCorrectedCandidate;
        rmm::device_uvector<int> d_indices_of_corrected_anchors;
        rmm::device_uvector<int> d_num_indices_of_corrected_anchors;
        rmm::device_uvector<int> d_indices_of_corrected_candidates;
        rmm::device_uvector<int> d_totalNumEdits;
        rmm::device_uvector<bool> d_isPairedCandidate;
        PinnedBuffer<bool> h_isPairedCandidate;

        rmm::device_uvector<int> d_numAnchors;
        rmm::device_uvector<int> d_numCandidates;
        rmm::device_uvector<read_number> d_anchorReadIds;
        rmm::device_uvector<unsigned int> d_anchor_sequences_data;
        rmm::device_uvector<int> d_anchor_sequences_lengths;
        rmm::device_uvector<read_number> d_candidate_read_ids;
        rmm::device_uvector<int> d_candidates_per_anchor;
        rmm::device_uvector<int> d_candidates_per_anchor_prefixsum; 

        PinnedBuffer<int> h_candidates_per_anchor_prefixsum; 
        PinnedBuffer<int> h_indices;

        PinnedBuffer<bool> h_flagsCandidates;

        std::unique_ptr<ManagedGPUMultiMSA> managedgpumsa;
    };



}
}






#endif
