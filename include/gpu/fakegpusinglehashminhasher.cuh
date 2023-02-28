#if 0

#ifndef CARE_FAKEGPU_SINGLE_HASH_MINHASHER_CUH
#define CARE_FAKEGPU_SINGLE_HASH_MINHASHER_CUH

#include <config.hpp>

#include <gpu/gpureadstorage.cuh>
#include <gpu/cuda_unique.cuh>
#include <cpuhashtable.hpp>
#include <gpu/gpuminhasher.cuh>
#include <groupbykey.hpp>
#include <gpu/cudaerrorcheck.cuh>
#include <gpu/gpusequencehasher.cuh>
#include <gpu/cubwrappers.cuh>

#include <cpusequencehasher.hpp>

#include <options.hpp>
#include <util.hpp>
#include <hpc_helpers.cuh>
#include <filehelpers.hpp>

#include <sequencehelpers.hpp>
#include <memorymanagement.hpp>
#include <threadpool.hpp>

#include <sharedmutex.hpp>


#include <vector>
#include <memory>
#include <limits>
#include <string>
#include <fstream>
#include <algorithm>

#include <cub/cub.cuh>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/device_vector.h>

#include <rmm/mr/device/device_memory_resource.hpp>
#include <rmm/mr/device/per_device_resource.hpp>
#include <rmm/device_uvector.hpp>
#include <gpu/rmm_utilities.cuh>

namespace care{
namespace gpu{

namespace fakegpusinglehashminhasherkernels{

    template<class T, class InputLengthIter, class InputOffsetIter, class OutputOffsetIter, int groupsize = 32>
    __global__
    void copyRangesKernel(
        const T* __restrict__ inputData,
        InputOffsetIter inputOffsets,
        InputLengthIter inputLengths,
        int numRanges,
        T* __restrict__ outputData,
        OutputOffsetIter outputOffsets
    ){
        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        const int numGroups = (blockDim.x * gridDim.x) / groupsize;

        for(int s = groupId; s < numRanges; s += numGroups){

            const int inputOffset = inputOffsets[s];
            const int outputOffset = outputOffsets[s];
            const int numElementsToCopy = inputLengths[s];

            care::gpu::memcpy<int>(
                group, 
                outputData + outputOffset, 
                inputData + inputOffset, 
                sizeof(int) * numElementsToCopy
            );
        }
    }

}

    /*
        Minhasher which can store query results in gpu memory and uses the gpu to parallelize some portions of the code
        However, hash tables reside on the host
    */
    class FakeGpuSingleHashMinhasher : public GpuMinhasher{
    public:
        using Key_t = GpuMinhasher::Key;
        using Value_t = read_number;
    private:
        using HashTable = CpuReadOnlyMultiValueHashTable<kmer_type, read_number>;

        using Range_t = std::pair<const Value_t*, const Value_t*>;
        //using SingleHashTableConstIter = std::unordered_multimap<kmer_type, read_number>::const_iterator;
        //using Range_t = std::pair<SingleHashTableConstIter, SingleHashTableConstIter>;

        struct QueryData{
            static constexpr int overprovisioningPercent = 0;
            
            template<class T>
            using PinnedBuffer = helpers::SimpleAllocationPinnedHost<T, overprovisioningPercent>;

            enum class Stage{
                None,
                NumValues,
                Retrieve
            };


            bool isInitialized = false;
            int deviceId;
            Stage previousStage = Stage::None;

            SetUnionHandle suHandle{};

            PinnedBuffer<kmer_type> h_minhashSignatures{};
            PinnedBuffer<read_number> h_candidate_read_ids_tmp{};
            PinnedBuffer<int> h_begin_offsets{};
            PinnedBuffer<int> h_end_offsets{};
            PinnedBuffer<int> h_numValuesPerSequence{};

            std::vector<Range_t> allRanges{};
            thrust::device_vector<int> d_hashFunctionNumbers{};

            MemoryUsage getMemoryInfo() const{
                MemoryUsage info;
                info.host = 0;
                info.device[deviceId] = 0;
    
                auto handlehost = [&](const auto& buff){
                    info.host += buff.capacityInBytes();
                };
    
                auto handledevice = [&](const auto& buff){
                    using ElementType = typename std::remove_reference<decltype(buff)>::type::value_type;
                    info.device[deviceId] += buff.size() * sizeof(ElementType);
                };

                auto handlevector = [&](const auto& buff){
                    info.host += 
                        sizeof(typename std::remove_reference<decltype(buff)>::type::value_type) * buff.capacity();
                };
    
                handlehost(h_minhashSignatures);
                handlehost(h_candidate_read_ids_tmp);
                handlehost(h_begin_offsets);
                handlehost(h_end_offsets);
                handlehost(h_numValuesPerSequence);
    
                handlevector(allRanges);
                handledevice(d_hashFunctionNumbers);
    
                return info;
            }

            void destroy(){
                cub::SwitchDevice sd{deviceId};

                h_minhashSignatures.destroy();
                h_candidate_read_ids_tmp.destroy();
                h_begin_offsets.destroy();
                h_end_offsets.destroy();
                h_numValuesPerSequence.destroy();

                allRanges.clear();
                allRanges.shrink_to_fit();

                d_hashFunctionNumbers.clear();
                d_hashFunctionNumbers.shrink_to_fit();

                isInitialized = false;
            }
        };

        
    public:

        FakeGpuSingleHashMinhasher() : FakeGpuSingleHashMinhasher(0, 50, 16, 0.8f){

        }

        FakeGpuSingleHashMinhasher(int maxNumKeys_, int maxValuesPerKey, int k, float loadfactor_)
            : loadfactor(loadfactor_), maxNumKeys(maxNumKeys_), kmerSize(k), resultsPerMapThreshold(maxValuesPerKey){

        }

        FakeGpuSingleHashMinhasher(const FakeGpuSingleHashMinhasher&) = delete;
        FakeGpuSingleHashMinhasher(FakeGpuSingleHashMinhasher&&) = default;
        FakeGpuSingleHashMinhasher& operator=(const FakeGpuSingleHashMinhasher&) = delete;
        FakeGpuSingleHashMinhasher& operator=(FakeGpuSingleHashMinhasher&&) = default;


        void constructFromReadStorage(
            const ProgramOptions& programOptions,
            std::uint64_t /*nReads*/,
            const GpuReadStorage& gpuReadStorage
        ){
            
            std::size_t memoryBeforeSinglehashtable = getAvailableMemoryInKB() * 1024;

            {
                auto& readStorage = gpuReadStorage;
                const auto& deviceIds = programOptions.deviceIds;

                int deviceId = deviceIds[0];

                cub::SwitchDevice sd{deviceId};

                const int requestedNumberOfMaps = programOptions.numHashFunctions;
                numSmallest = requestedNumberOfMaps;

                const read_number numReads = readStorage.getNumberOfReads();
                const int maximumSequenceLength = readStorage.getSequenceLengthUpperBound();

                auto sequencehandle = gpuReadStorage.makeHandle();
                const std::size_t encodedSequencePitchInInts = SequenceHelpers::getEncodedNumInts2Bit(maximumSequenceLength);

                rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource();

                constexpr read_number parallelReads = 1000000;
                //constexpr read_number parallelReads = 20;
                const int numIters = SDIV(numReads, parallelReads);

                const MemoryUsage memoryUsageOfReadStorage = readStorage.getMemoryInfo();
                std::size_t totalLimit = programOptions.memoryTotalLimit;
                if(totalLimit > memoryUsageOfReadStorage.host){
                    totalLimit -= memoryUsageOfReadStorage.host;
                }else{
                    totalLimit = 0;
                }
                if(totalLimit == 0){
                    throw std::runtime_error("Not enough memory available for hash tables. Abort!");
                }
                std::size_t maxMemoryForTables = getAvailableMemoryInKB() * 1024;
                // std::cerr << "available: " << maxMemoryForTables 
                //         << ",memoryForHashtables: " << programOptions.memoryForHashtables
                //         << ", memoryTotalLimit: " << programOptions.memoryTotalLimit
                //         << ", rsHostUsage: " << memoryUsageOfReadStorage.host << "\n";

                std::size_t maxNumPairs = std::size_t(numReads) * std::size_t(numSmallest);

                kvtable = std::make_unique<DoublePassMultiValueHashTable<kmer_type, read_number>>(maxNumPairs / 8, loadfactor);

                cudaStream_t stream = cudaStreamPerThread;
                
                rmm::device_uvector<unsigned int> d_sequenceData(encodedSequencePitchInInts * parallelReads, stream, mr);
                rmm::device_uvector<int> d_lengths(parallelReads, stream, mr);
                rmm::device_uvector<read_number> d_indices(parallelReads, stream, mr);
                
                helpers::SimpleAllocationPinnedHost<read_number, 0> h_indices(parallelReads);
                setMemoryLimitForConstruction(maxMemoryForTables);
                
                std::cout << "Constructing the single hash table\n";

                std::size_t numKeys = 0;

                helpers::CpuTimer firstpasstimer("firstpass");

                for (int iter = 0; iter < numIters; iter++){
                    read_number readIdBegin = iter * parallelReads;
                    read_number readIdEnd = std::min((iter + 1) * parallelReads, numReads);

                    const std::size_t curBatchsize = readIdEnd - readIdBegin;

                    std::iota(h_indices.get(), h_indices.get() + curBatchsize, readIdBegin);

                    CUDACHECK(cudaMemcpyAsync(d_indices.data(), h_indices, sizeof(read_number) * curBatchsize, H2D, stream));

                    gpuReadStorage.gatherSequences(
                        sequencehandle,
                        d_sequenceData.data(),
                        encodedSequencePitchInInts,
                        makeAsyncConstBufferWrapper(h_indices.data()),
                        d_indices.data(),
                        curBatchsize,
                        stream,
                        mr
                    );
                
                    gpuReadStorage.gatherSequenceLengths(
                        sequencehandle,
                        d_lengths.data(),
                        d_indices.data(),
                        curBatchsize,
                        stream
                    );

                    std::vector<kmer_type> tmpkeys(curBatchsize * numSmallest);
                    std::vector<read_number> tmpids(curBatchsize * numSmallest);

                    auto numNewKeys = computeHashesAndReadIds(
                        tmpkeys.data(),
                        tmpids.data(),
                        d_sequenceData.data(),
                        curBatchsize,
                        d_lengths.data(),
                        encodedSequencePitchInInts,
                        d_indices.data(),
                        h_indices,
                        stream
                    );

                    numKeys += numNewKeys;
                    tmpkeys.erase(tmpkeys.begin() + numNewKeys, tmpkeys.end());
                    tmpids.erase(tmpids.begin() + numNewKeys, tmpids.end());

                    CUDACHECK(cudaStreamSynchronize(stream));

                    kvtable->firstPassInsert(tmpkeys.data(), tmpids.data(), tmpkeys.size());
                }

                CUDACHECK(cudaStreamSynchronize(stream));

                firstpasstimer.print();

                //kvtable->firstPassDone(2, 75);
                kvtable->firstPassDone(2, resultsPerMapThreshold);

                helpers::CpuTimer secondPassTimer("secondpass");

                for (int iter = 0; iter < numIters; iter++){
                    read_number readIdBegin = iter * parallelReads;
                    read_number readIdEnd = std::min((iter + 1) * parallelReads, numReads);

                    const std::size_t curBatchsize = readIdEnd - readIdBegin;

                    std::iota(h_indices.get(), h_indices.get() + curBatchsize, readIdBegin);

                    CUDACHECK(cudaMemcpyAsync(d_indices.data(), h_indices, sizeof(read_number) * curBatchsize, H2D, stream));

                    gpuReadStorage.gatherSequences(
                        sequencehandle,
                        d_sequenceData.data(),
                        encodedSequencePitchInInts,
                        makeAsyncConstBufferWrapper(h_indices.data()),
                        d_indices.data(),
                        curBatchsize,
                        stream,
                        mr
                    );
                
                    gpuReadStorage.gatherSequenceLengths(
                        sequencehandle,
                        d_lengths.data(),
                        d_indices.data(),
                        curBatchsize,
                        stream
                    );

                    std::vector<kmer_type> tmpkeys(curBatchsize * numSmallest);
                    std::vector<read_number> tmpids(curBatchsize * numSmallest);

                    auto numNewKeys = computeHashesAndReadIds(
                        tmpkeys.data(),
                        tmpids.data(),
                        d_sequenceData.data(),
                        curBatchsize,
                        d_lengths.data(),
                        encodedSequencePitchInInts,
                        d_indices.data(),
                        h_indices,
                        stream
                    );

                    numKeys += numNewKeys;
                    tmpkeys.erase(tmpkeys.begin() + numNewKeys, tmpkeys.end());
                    tmpids.erase(tmpids.begin() + numNewKeys, tmpids.end());

                    CUDACHECK(cudaStreamSynchronize(stream));

                    kvtable->secondPassInsert(tmpkeys.data(), tmpids.data(), tmpkeys.size());
                }


                secondPassTimer.print();

                kvtable->secondPassDone();

                std::cerr << "numKeys = " << numKeys << "\n";

                setThreadPool(nullptr); 
                
                gpuReadStorage.destroyHandle(sequencehandle);
            }
        }
    

        MinhasherHandle makeMinhasherHandle() const override {
            auto data = std::make_unique<QueryData>();
            data->d_hashFunctionNumbers.resize(getNumberOfMaps());

            thrust::sequence(thrust::device, data->d_hashFunctionNumbers.begin(), data->d_hashFunctionNumbers.end(), 0);

            CUDACHECK(cudaGetDevice(&data->deviceId));
            data->isInitialized = true;

            //std::unique_lock<std::shared_mutex> lock(sharedmutex);
            std::unique_lock<SharedMutex> lock(sharedmutex);
            const int handleid = counter++;
            MinhasherHandle h = constructHandle(handleid);

            tempdataVector.emplace_back(std::move(data));
            return h;
        }

        void destroyHandle(MinhasherHandle& handle) const override{

            std::unique_lock<SharedMutex> lock(sharedmutex);

            const int id = handle.getId();
            assert(id < int(tempdataVector.size()));
            
            tempdataVector[id] = nullptr;
            handle = constructHandle(std::numeric_limits<int>::max());
        }

        void determineNumValues(
            MinhasherHandle& queryHandle,
            const unsigned int* d_sequenceData2Bit,
            std::size_t encodedSequencePitchInInts,
            const int* d_sequenceLengths,
            int numSequences,
            int* d_numValuesPerSequence,
            int& totalNumValues,
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr
        ) const override {
            QueryData* const queryData = getQueryDataFromHandle(queryHandle);

            DEBUGSTREAMSYNC(stream);

            assert(queryData->isInitialized);
            if(numSequences == 0) return;

            // queryData->h_minhashSignatures.resize(numSmallest * numSequences);
            queryData->h_numValuesPerSequence.resize(numSequences);
            std::vector<Range_t>& allRanges = queryData->allRanges;

            allRanges.resize(numSmallest * numSequences);

            GPUSequenceHasher<kmer_type> hasher;

            auto hashResult = hasher.getTopSmallestKmerHashes(
                d_sequenceData2Bit,
                encodedSequencePitchInInts,
                numSequences,
                d_sequenceLengths,
                getKmerSize(),
                numSmallest,
                stream,
                mr
            );

            std::vector<kmer_type> h_hashes(numSequences * numSmallest);
            std::vector<int> h_numHashesPerSequence(numSequences);

            CUDACHECK(cudaMemcpyAsync(
                h_hashes.data(), 
                hashResult.d_hashvalues.data(), 
                sizeof(kmer_type) * numSmallest * numSequences, 
                D2H, 
                stream
            ));

            CUDACHECK(cudaMemcpyAsync(
                h_numHashesPerSequence.data(), 
                hashResult.d_numPerSequences.data(), 
                sizeof(int) * numSequences, 
                D2H, 
                stream
            ));

            CUDACHECK(cudaStreamSynchronize(stream));

            auto hashesPtr = h_hashes.data();
            for(int s = 0; s < numSequences; s++){
                const int numHashes = h_numHashesPerSequence[s];
                assert(numHashes <= numSmallest);

                int numValues = 0;

                for(int i = 0; i < numHashes; i++){
                    const kmer_type key = *(hashesPtr++);
                    const auto qr = kvtable->query(key);
                    allRanges[s * numSmallest + i] = {qr.valuesBegin, qr.valuesBegin + qr.numValues};
                    numValues += qr.numValues;
                }

                for(int i = numHashes; i < numSmallest; i++){
                    allRanges[s * numSmallest + i] = {nullptr, nullptr};
                }

                queryData->h_numValuesPerSequence[s] = numValues;

                totalNumValues += numValues;
            }

            DEBUGSTREAMSYNC(stream);

            CUDACHECK(cudaMemcpyAsync(d_numValuesPerSequence, queryData->h_numValuesPerSequence.data(), sizeof(int) * numSequences, H2D, stream));

            DEBUGSTREAMSYNC(stream);

            queryData->previousStage = QueryData::Stage::NumValues;
        }

        void retrieveValues(
            MinhasherHandle& queryHandle,
            const read_number* d_readIds,
            int numSequences,
            int totalNumValues,
            read_number* d_values,
            int* d_numValuesPerSequence,
            int* d_offsets, //numSequences + 1
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr
        ) const override {
            QueryData* const queryData = getQueryDataFromHandle(queryHandle);

            DEBUGSTREAMSYNC(stream);

            assert(queryData->isInitialized);
            if(numSequences == 0) return;

            assert(queryData->previousStage == QueryData::Stage::NumValues);

            //std::cerr << "totalNumValues: " << totalNumValues << "\n";

            if(totalNumValues == 0){
                cudaMemsetAsync(d_numValuesPerSequence, 0, sizeof(int) * numSequences, stream);
                cudaMemsetAsync(d_offsets, 0, sizeof(int) * (numSequences + 1), stream);

                queryData->previousStage = QueryData::Stage::Retrieve;

                DEBUGSTREAMSYNC(stream);

                return;
            }

            constexpr int roundUpTo = 10000;
            const int roundedTotalNum = SDIV(totalNumValues, roundUpTo) * roundUpTo;
            queryData->h_candidate_read_ids_tmp.resize(roundedTotalNum);

            rmm::device_uvector<read_number> d_candidate_read_ids_tmp(roundedTotalNum, stream, mr);

            queryData->h_begin_offsets.resize(numSequences+1);
            queryData->h_end_offsets.resize(numSequences+1);


            //results will be in Current() buffer
            cub::DoubleBuffer<read_number> d_values_dblbuf(d_values, d_candidate_read_ids_tmp.data());

            auto copyHitsToPinnedMemory = [queryData, numSequences, minhasher = this](){
                int* h_numValuesPerSequence = queryData->h_numValuesPerSequence.data();

                queryData->h_begin_offsets[0] = 0;
                std::partial_sum(
                    h_numValuesPerSequence, 
                    h_numValuesPerSequence + numSequences - 1, 
                    queryData->h_begin_offsets.begin() + 1
                );
                std::partial_sum(
                    h_numValuesPerSequence, 
                    h_numValuesPerSequence + numSequences, 
                    queryData->h_end_offsets.begin()
                );

                for(int s = 0; s < numSequences; s++){                    

                    //prefetch first element of next range if the next range is not empty
                    // {
                    //     constexpr int nextprefetch = 2;
                    //     int prefetchSequence = 0;
                    //     int prefetchMap = 0;
                    //     if(i + nextprefetch < numSequences){
                    //         prefetchMap = map;
                    //         prefetchSequence = i + nextprefetch;
                    //     }else{
                    //         prefetchMap = map + 1;
                    //         prefetchSequence = i + nextprefetch - numSequences;
                    //     }

                    //     if(prefetchMap < minhasher->getNumberOfMaps()){
                    //         const auto& range = queryData->allRanges[prefetchMap * numSequences + prefetchSequence];
                    //         if(range.first != range.second){
                    //             __builtin_prefetch(range.first, 0, 0);
                    //         }
                    //     }
                    // }

                    auto destination = queryData->h_candidate_read_ids_tmp.data() + queryData->h_begin_offsets[s];

                    for(int i = 0; i < minhasher->numSmallest; i++){

                        const auto& range = queryData->allRanges[s * minhasher->numSmallest + i];

                        destination = std::copy(
                            range.first, 
                            range.second, 
                            destination
                        );
                    }

                    assert(destination == queryData->h_candidate_read_ids_tmp.data() + queryData->h_end_offsets[s]);
                }
            };

            

            copyHitsToPinnedMemory();

            DEBUGSTREAMSYNC(stream);

            CUDACHECK(cudaMemcpyAsync(
                d_values_dblbuf.Current(),
                queryData->h_candidate_read_ids_tmp.data(),
                sizeof(read_number) * totalNumValues,
                H2D,
                stream
            ));

            DEBUGSTREAMSYNC(stream);

            rmm::device_uvector<int> d_begin_offsets((numSequences + 1), stream, mr);
            rmm::device_uvector<int> d_end_offsets((numSequences + 1), stream, mr);
            rmm::device_uvector<int> d_global_begin_offsets((numSequences), stream, mr);

            DEBUGSTREAMSYNC(stream);

            //copy h_endoffsets to d_endoffsets.
            //Then copy d_endoffsets to d_begin_offsets shifted to the right by 1.
            //Then copy d_begin_offsets to d_global_begin_offsets
            helpers::lambda_kernel<<<SDIV(numSequences, 1024), 1024, 0, stream>>>(
                [
                    numSequences,
                    h_end_offsets = queryData->h_end_offsets.data(),
                    d_end_offsets = d_end_offsets.data(),
                    d_begin_offsets = d_begin_offsets.data(),
                    d_global_begin_offsets = d_global_begin_offsets.data()
                ] __device__ (){
                    const int tid = threadIdx.x + blockIdx.x * blockDim.x;

                    if(tid < numSequences){
                        const int data = h_end_offsets[tid];
                        d_end_offsets[tid] = data;

                        if(tid < numSequences - 1){
                            d_begin_offsets[tid + 1] = data;
                            d_global_begin_offsets[tid + 1] = data;
                        }
                        
                        if(tid == 0){
                            d_begin_offsets[0] = 0;
                            d_global_begin_offsets[0] = 0;
                        }
                    }
                }
            );

            DEBUGSTREAMSYNC(stream);

            GpuSegmentedUnique::unique(
                d_values_dblbuf.Current(), //input
                totalNumValues,
                d_values_dblbuf.Alternate(), //output
                d_numValuesPerSequence,
                numSequences,
                d_begin_offsets.data(), //device accessible
                d_end_offsets.data(), //device accessible
                queryData->h_begin_offsets,
                queryData->h_end_offsets,
                0,
                sizeof(read_number) * 8,
                stream,
                mr
            );

            DEBUGSTREAMSYNC(stream);

            ::destroy(d_begin_offsets, stream);
            ::destroy(d_end_offsets, stream);

            if(d_readIds != nullptr){

                //remove self read ids (inplace)
                //--------------------------------------------------------------------
                callFindAndRemoveFromSegmentKernel<read_number,128,4>(
                    d_readIds,
                    d_values_dblbuf.Alternate(),
                    numSequences,
                    d_numValuesPerSequence,
                    d_global_begin_offsets.data(),
                    stream
                );

                DEBUGSTREAMSYNC(stream);

            }

            int* d_newOffsets = d_offsets;

            CUDACHECK(cudaMemsetAsync(d_newOffsets, 0, sizeof(int), stream));

            DEBUGSTREAMSYNC(stream);

            std::size_t cubtempbytes = 0;
            CUDACHECK(cub::DeviceScan::InclusiveSum(
                nullptr,
                cubtempbytes,
                d_numValuesPerSequence,
                d_newOffsets + 1,
                numSequences,
                stream
            ));

            rmm::device_uvector<char> d_cubTemp(cubtempbytes, stream);

            CUDACHECK(cub::DeviceScan::InclusiveSum(
                d_cubTemp.data(),
                cubtempbytes,
                d_numValuesPerSequence,
                d_newOffsets + 1,
                numSequences,
                stream
            ));

            DEBUGSTREAMSYNC(stream);

            //copy final remaining values into contiguous range
            helpers::lambda_kernel<<<numSequences, 128, 0, stream>>>(
                [
                    d_values_in = d_values_dblbuf.Alternate(),
                    d_values_out = d_values_dblbuf.Current(),
                    numSequences,
                    d_numValuesPerSequence,
                    d_offsets = d_global_begin_offsets.data(),
                    d_newOffsets
                ] __device__ (){

                    for(int s = blockIdx.x; s < numSequences; s += gridDim.x){
                        const int numValues = d_numValuesPerSequence[s];
                        const int inOffset = d_offsets[s];
                        const int outOffset = d_newOffsets[s];

                        for(int c = threadIdx.x; c < numValues; c += blockDim.x){
                            d_values_out[outOffset + c] = d_values_in[inOffset + c];    
                        }
                    }
                }
            ); CUDACHECKASYNC;

            DEBUGSTREAMSYNC(stream);

            queryData->previousStage = QueryData::Stage::Retrieve;
        }

        MemoryUsage getMemoryInfo() const noexcept override{
            MemoryUsage result;

            result = kvtable->getMemoryInfo();

            return result;
        }

        MemoryUsage getMemoryInfo(const MinhasherHandle& handle) const noexcept override{
            return getQueryDataFromHandle(handle)->getMemoryInfo();
        }

        int getNumResultsPerMapThreshold() const noexcept override{
            return resultsPerMapThreshold;
        }
        
        int getNumberOfMaps() const noexcept override{
            return 1;
        }

        void destroy() {
            kvtable = nullptr;
        }

        int getKmerSize() const noexcept override{
            return kmerSize;
        }

        bool hasGpuTables() const noexcept override {
            return false;
        }

        std::uint64_t getKmerMask() const{
            constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;

            return std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - getKmerSize()) * 2);
        }

        void writeToStream(std::ostream& os) const{

            os.write(reinterpret_cast<const char*>(&kmerSize), sizeof(int));
            os.write(reinterpret_cast<const char*>(&numSmallest), sizeof(int));
            os.write(reinterpret_cast<const char*>(&resultsPerMapThreshold), sizeof(int));

            os.write(reinterpret_cast<const char*>(&loadfactor), sizeof(float));

            kvtable->writeToStream(os);
        }

        int loadFromStream(std::ifstream& is, int /*numMapsUpperLimit = std::numeric_limits<int>::max()*/){
            destroy();

            is.read(reinterpret_cast<char*>(&kmerSize), sizeof(int));
            is.read(reinterpret_cast<char*>(&numSmallest), sizeof(int));
            is.read(reinterpret_cast<char*>(&resultsPerMapThreshold), sizeof(int));

            is.read(reinterpret_cast<char*>(&loadfactor), sizeof(float));

            kvtable = std::make_unique<DoublePassMultiValueHashTable<kmer_type, read_number>>(1, loadfactor);
            kvtable->loadFromStream(is);

            return 0;
        }

        std::size_t computeHashesAndReadIds(
            kmer_type* h_keyoutput,
            read_number* h_readIdsOutput,
            const unsigned int* d_sequenceData2Bit,
            int numSequences,
            const int* d_sequenceLengths,
            std::size_t encodedSequencePitchInInts,
            const read_number* /*d_readIds*/,
            const read_number* h_readIds,
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
        ){

            CubCallWrapper cub(mr);
            GPUSequenceHasher<kmer_type> hasher;

            const bool debug = false;

            auto hashResult = hasher.getTopSmallestKmerHashes(
                d_sequenceData2Bit,
                encodedSequencePitchInInts,
                numSequences,
                d_sequenceLengths,
                getKmerSize(),
                numSmallest,
                stream,
                mr,
                debug
            );

            std::vector<int> h_numHashesPerSequence(numSequences);

            CUDACHECK(cudaMemcpyAsync(
                h_keyoutput, 
                hashResult.d_hashvalues.data(), 
                sizeof(kmer_type) * numSmallest * numSequences, 
                D2H, 
                stream
            ));

            CUDACHECK(cudaMemcpyAsync(
                h_numHashesPerSequence.data(), 
                hashResult.d_numPerSequences.data(), 
                sizeof(int) * numSequences, 
                D2H, 
                stream
            ));

            CUDACHECK(cudaStreamSynchronize(stream));

            std::size_t numKeys = 0;

            auto iterator = h_readIdsOutput;
            for(int s = 0; s < numSequences; s++){
                const int num = h_numHashesPerSequence[s];
                std::fill(iterator, iterator + num, h_readIds[s]);
                iterator += num;
                numKeys += num;
            }

            return numKeys;
        }   

        void setThreadPool(ThreadPool* tp){
            threadPool = tp;
        }

        void setMemoryLimitForConstruction(std::size_t limit){
            memoryLimit = limit;
        }

    private:

        
        QueryData* getQueryDataFromHandle(const MinhasherHandle& queryHandle) const{
            std::shared_lock<SharedMutex> lock(sharedmutex);

            return tempdataVector[queryHandle.getId()].get();
        }
        

        mutable int counter = 0;
        mutable SharedMutex sharedmutex{};

        float loadfactor = 0.8f;
        int numSmallest{};
        int maxNumKeys{};
        int kmerSize{};
        int resultsPerMapThreshold{};
        ThreadPool* threadPool;
        std::size_t memoryLimit;
        std::size_t singlehashtablememoryusage{};
        std::vector<std::unique_ptr<HashTable>> minhashTables{};

        mutable std::vector<std::unique_ptr<QueryData>> tempdataVector{};

        std::unique_ptr<DoublePassMultiValueHashTable<kmer_type, read_number>> kvtable{};
    };






    
}
}



#endif

#endif
