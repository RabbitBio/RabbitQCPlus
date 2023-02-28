#ifdef CARE_HAS_WARPCORE

#ifndef CARE_SINGLE_GPU_MINHASHER_CUH
#define CARE_SINGLE_GPU_MINHASHER_CUH

#include <config.hpp>

#include <gpu/gpureadstorage.cuh>
#include <gpu/cuda_unique.cuh>
#include <cpuhashtable.hpp>
#include <gpu/gpuhashtable.cuh>
#include <gpu/gpuminhasher.cuh>
#include <gpu/cudaerrorcheck.cuh>
#include <gpu/cubwrappers.cuh>
#include <gpu/gpusequencehasher.cuh>

#include <options.hpp>
#include <util.hpp>
#include <hpc_helpers.cuh>
#include <filehelpers.hpp>

#include <sequencehelpers.hpp>
#include <memorymanagement.hpp>
#include <threadpool.hpp>
#include <sharedmutex.hpp>

#include <cub/cub.cuh>

#include <vector>
#include <memory>
#include <limits>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <mutex>

#include <rmm/mr/device/device_memory_resource.hpp>
#include <rmm/mr/device/per_device_resource.hpp>
#include <rmm/device_uvector.hpp>
#include <rmm/device_scalar.hpp>
#include <gpu/rmm_utilities.cuh>

namespace care{
namespace gpu{

    namespace sgpuminhasherkernels{

        /*
            Insert kv-pairs [i * numKeysPerTable, (i+1) * numKeysPerTable] into table i .
            There are numKeysPerTable * numTables keys, and numKeysPerTable values.
            Values are the same for each table
        */
 
        template<class DeviceTableInsertView, class Key, class Value>
        __global__
        void insertIntoTablesKernel(
            DeviceTableInsertView* __restrict__ tables,
            const int numTables,
            const Key* __restrict__ keys,
            const bool* __restrict__ isValid,
            const int numKeysPerTable,
            const Value* __restrict__ values
        ){
            const int tid = threadIdx.x + blockDim.x * blockIdx.x;
            const int stride = blockDim.x * gridDim.x;

            constexpr int tilesize = DeviceTableInsertView::cg_size();
            assert(stride % tilesize == 0);

            constexpr int warpsize = 32;
            constexpr int tilesPerWarp = warpsize / tilesize;

            auto warp = cg::tiled_partition<warpsize>(cg::this_thread_block());
            auto tile = cg::tiled_partition<tilesize>(warp);

            const int warpId = tid / warpsize;
            const int numWarpsInGrid = stride / warpsize;
            const int warpsPerTable = SDIV(numKeysPerTable, tilesPerWarp);

            for(int w = warpId; w < warpsPerTable * numTables; w += numWarpsInGrid){
                const int tableIndex = w / warpsPerTable;
                const int valueIndex = tilesPerWarp * (w % warpsPerTable) + tile.meta_group_rank();
                const int keyIndex = numKeysPerTable * tableIndex + valueIndex;
                if(keyIndex < numKeysPerTable * numTables){
                    if(isValid[keyIndex]){
                        DeviceTableInsertView table = tables[tableIndex];
                        table.insert(keys[keyIndex], values[valueIndex], tile);
                    }
                }
                //ensure that different groups in the same warp do not operate on different hashtables (warpcore issue)
                warp.sync();
            }
        }

        template<class DeviceTableInsertView, class Key>
        __global__
        void hashAndFixAndInsertIntoTablesKernel(
            DeviceTableInsertView* __restrict__ tables,
            const int numTables,
            const unsigned int* __restrict__ sequences2Bit,
            std::size_t encodedSequencePitchInInts,
            const int numSequences,
            const int* __restrict__ sequenceLengths,
            const read_number* __restrict__ readIds,
            int k,
            const int* __restrict__ hashFunctionNumbers
        ){
            const int tid = threadIdx.x + blockDim.x * blockIdx.x;
            const int stride = blockDim.x * gridDim.x;

            constexpr int tilesize = DeviceTableInsertView::cg_size();
            constexpr int warpsize = 32;
            constexpr int tilesPerWarp = warpsize / tilesize;

            auto warp = cg::tiled_partition<warpsize>(cg::this_thread_block());
            auto tile = cg::tiled_partition<tilesize>(warp);
            const int numWarpsInGrid = stride / warpsize;
            const int warpId = tid / warpsize;

            constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;
            const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);

            static_assert(tilesize == 8, "tilesize not 8. only tested for 8");

            //use 1 warp per tilesPerWarp sequences
            const int maxNumSequencesPerWarp = tilesPerWarp;
            const int numWarpIterations = SDIV(numSequences, maxNumSequencesPerWarp);


            for(int w = warpId; w < numWarpIterations; w += numWarpsInGrid){
                const int sequenceOffset = maxNumSequencesPerWarp * w;
                const int numWarpSequences = std::min(maxNumSequencesPerWarp, numSequences - sequenceOffset);
                //process sequences [sequenceOffset, sequenceOffset + numWarpSequences)

                constexpr int regKeys = 1;
                Key keys[regKeys];
                bool keyIsValid[regKeys];

                constexpr int xStride = 8 * regKeys;

                /*
                thread mapping for regKeys == 1
                    h0 h1 h2 h3 h4 h5 h6 h7
                s0  0  1  2  3  4  5  6  7
                s1  8  9 10 11 12 13 14 15
                s2 16 17 18 19 20 21 22 23
                s3 24 25 26 27 28 19 30 31
                */

                /*
                thread mapping for regKeys == 2
                    h0 h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13 h14 h15
                s0   0  0  1  1  2  2  3  3  4  4  5   5   6   6   7   7
                s1   8  8  9  9  10 10 11 11 12 12 13  13  14  14  15  15
                s2   16 16 17 17 18 18 19 19 20 20 21  21  22  22  23  23
                s3   24 24 25 25 26 26 27 27 28 28 29  29  30  30  31  31
                */


                for(int x = 0; x < numTables; x += xStride){
                    const int numHashesToComputePerSequence = std::min(xStride, numTables - x);
                    //compute hashes [x, x+numHashesToComputePerSequence) for sequences [sequenceOffset, sequenceOffset + numWarpSequences)

                    const int s = warp.thread_rank() / 8;
                    const unsigned int* const sequenceData = sequences2Bit + encodedSequencePitchInInts * (sequenceOffset + s);
                    const int sequenceLength = sequenceLengths[(sequenceOffset + s)];

                    if(sequenceLength >= k){
                        std::uint64_t minHashValue[regKeys];
                        int hashFuncId[regKeys];
                        #pragma unroll
                        for(int c = 0; c < regKeys; c++){
                            minHashValue[c] = std::numeric_limits<std::uint64_t>::max();

                            const int h = regKeys*(warp.thread_rank() % 8) + c;
                            if(h < numTables){
                                hashFuncId[c] = hashFunctionNumbers[x+h];
                            }else{
                                hashFuncId[c] = 0;
                            }
                        }

                        SequenceHelpers::forEachEncodedCanonicalKmerFromEncodedSequence(
                            sequenceData,
                            sequenceLength,
                            k,
                            [&](std::uint64_t kmer, int /*pos*/){
                                using hasher = hashers::MurmurHash<std::uint64_t>;

                                #pragma unroll
                                for(int c = 0; c < regKeys; c++){
                                    const auto hashvalue = hasher::hash(kmer + hashFuncId[c]);
                                    minHashValue[c] = min(minHashValue[c], hashvalue);
                                }
                            }
                        );

                        #pragma unroll
                        for(int c = 0; c < regKeys; c++){
                            keys[c] = Key(minHashValue[c] & kmer_mask);
                            GpuHashtableKeyCheck<Key, read_number> isValidKey;

                            while(!isValidKey(keys[c])){
                                keys[c]++;
                            }
                            keyIsValid[c] = true;
                        }

                    }else{
                        for(int i = 0; i < regKeys; i++){
                            keyIsValid[i] = false;
                        }
                    }

                    //insert keys
                    for(int t = 0; t < numHashesToComputePerSequence; t += regKeys){
                        const int shflSrcIndex = tile.meta_group_rank() * 8 + t / regKeys;
                        const int tableOffset = x + t;

                        for(int c = 0; c < regKeys; c++){                            
                            const Key keyToInsert = warp.shfl(keys[c], shflSrcIndex);
                            const bool keyIsValidToInsert = warp.shfl(keyIsValid[c], shflSrcIndex);
                            const read_number valueToInsert = readIds[(sequenceOffset + s)];
                            if(keyIsValidToInsert){
                                DeviceTableInsertView table = tables[tableOffset + c];
                                table.insert(keyToInsert, valueToInsert, tile);
                            }
                        }
                    }
                }
            }
        }


    }

    class MultiGpuMinhasher; //forward declaration

    class SingleGpuMinhasher : public GpuMinhasher{
        friend class MultiGpuMinhasher;
    public:
        using Key = kmer_type;
        using Value = read_number;
    private:
        using GpuTable = GpuHashtable<Key, Value>;

        using DeviceSwitcher = cub::SwitchDevice;

        template<class T>
        using HostBuffer = helpers::SimpleAllocationPinnedHost<T, 5>;
        struct QueryData{
            enum class Stage{
                None,
                NumValues,
                Retrieve
            };

            int deviceId{};
            Stage previousStage = Stage::None;
            int* d_numValuesPerSequence{};
            rmm::device_uvector<char> d_singlepersistentbuffer;

            QueryData(rmm::mr::device_memory_resource* mr)
            : d_singlepersistentbuffer(0, cudaStreamPerThread, mr)
            {
                CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));
            }

            MemoryUsage getMemoryInfo() const{
                MemoryUsage mem{};

                auto handledevice = [&](const auto& buff){
                    using ElementType = typename std::remove_reference<decltype(buff)>::type::value_type;
                    mem.device[deviceId] += buff.size() * sizeof(ElementType);
                };

                handledevice(d_singlepersistentbuffer);

                return mem;
            }
        };

    public:

        SingleGpuMinhasher(int maxNumKeys_, int maxValuesPerKey, int k, float loadfactor_)
            : maxNumKeys(maxNumKeys_), kmerSize(k), resultsPerMapThreshold(maxValuesPerKey), loadfactor(loadfactor_)
        {
            CUDACHECK(cudaGetDevice(&deviceId));
        }

        std::unique_ptr<SingleGpuMinhasher> makeCopy(int targetDeviceId) const{
            DeviceSwitcher ds(targetDeviceId);

            auto result = std::make_unique<SingleGpuMinhasher>(0,0,0, loadfactor);
            if(!result) return nullptr;
            
            result->maxNumKeys = maxNumKeys;
            result->kmerSize = kmerSize;
            result->resultsPerMapThreshold = resultsPerMapThreshold;
            result->h_currentHashFunctionNumbers.resize(h_currentHashFunctionNumbers.size());
            std::copy(h_currentHashFunctionNumbers.begin(), h_currentHashFunctionNumbers.end(), result->h_currentHashFunctionNumbers.begin());

            std::size_t requiredTempBytes = 0;
            for(const auto& ptr : gpuHashTables){
                std::size_t bytes = ptr->getMakeCopyTempBytes();
                requiredTempBytes = std::max(requiredTempBytes, bytes);
            }

            thrust::device_vector<char> d_copytemp(requiredTempBytes);

            for(const auto& ptr : gpuHashTables){
                auto newtableptr = ptr->makeCopy(thrust::raw_pointer_cast(d_copytemp.data()), targetDeviceId);
                if(newtableptr){
                    result->gpuHashTables.push_back(std::move(newtableptr));
                }else{
                    cudaGetLastError();
                    return nullptr;
                }
            }

            std::vector<GpuTable::DeviceTableView> views;
            for(const auto& ptr : result->gpuHashTables){
                views.emplace_back(ptr->makeDeviceView());
            }

            result->d_deviceAccessibleTableViews.resize(views.size());
            CUDACHECK(cudaMemcpyAsync(
                result->d_deviceAccessibleTableViews.data(),
                views.data(),
                sizeof(GpuTable::DeviceTableView) * views.size(),
                H2D,
                cudaStreamPerThread
            ));

            CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));

            return result;
        }


        int addHashTables(int numAdditionalTables, const int* hashFunctionIds, cudaStream_t stream) override {
            
            DeviceSwitcher ds(deviceId);

            int added = 0;
            int cur = gpuHashTables.size();

            assert(!(numAdditionalTables + cur > 64));
            std::vector<int> tmpNumbers(h_currentHashFunctionNumbers.begin(), h_currentHashFunctionNumbers.end());

            for(int i = 0; i < numAdditionalTables; i++){
                auto ptr = std::make_unique<GpuTable>(std::size_t(maxNumKeys / getLoad()),
                    getLoad(),
                    resultsPerMapThreshold,
                    stream
                );

                auto status = ptr->pop_status(stream);
                CUDACHECK(cudaStreamSynchronize(stream));
                if(status.has_any_errors()){
                    if(!status.has_out_of_memory()){
                        std::cerr << "observed error when initializing hash function " << (gpuHashTables.size() + 1) << " : " << i << ", " << status << "\n";
                    }
                    break;
                }else{

                    assert(!status.has_any_errors()); 

                    gpuHashTables.emplace_back(std::move(ptr));

                    added++;
                    tmpNumbers.push_back(hashFunctionIds[i]);
                }
            }

            h_currentHashFunctionNumbers.resize(tmpNumbers.size());
            std::copy(tmpNumbers.begin(), tmpNumbers.end(), h_currentHashFunctionNumbers.begin());

            return added;
        }

        void insert(
            const unsigned int* d_sequenceData2Bit,
            int numSequences,
            const int* d_sequenceLengths,
            std::size_t encodedSequencePitchInInts,
            const read_number* d_readIds,
            const read_number* /*h_readIds*/,
            int firstHashfunction,
            int numHashFunctions,
            const int* h_hashFunctionNumbers,
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr
        ) override {

            const std::size_t signaturesRowPitchElements = numHashFunctions;

            assert(firstHashfunction + numHashFunctions <= int(gpuHashTables.size()));

            DeviceSwitcher ds(deviceId);

            rmm::device_uvector<int> d_hashFunctionNumbers(numHashFunctions, stream, mr);
            
            CUDACHECK(cudaMemcpyAsync(
                d_hashFunctionNumbers.data(), 
                h_hashFunctionNumbers, 
                sizeof(int) * numHashFunctions, 
                H2D, 
                stream
            ));

            #if 0

            h_insertTemp.resize(numHashFunctions);
            d_insertTemp.resize(numHashFunctions);

            //create device views
            for(int i = 0; i < numHashFunctions; i++){
                new (&h_insertTemp[i]) GpuTable::DeviceTableInsertView(gpuHashTables[firstHashfunction + i]->makeDeviceInsertView());
            }
            //transfer to device
            CUDACHECK(cudaMemcpyAsync(
                d_insertTemp.data(),
                h_insertTemp.data(),
                sizeof(GpuTable::DeviceTableInsertView) * numHashFunctions,
                H2D,
                stream
            ));

            callHashAndInsertKernel(
                d_insertTemp.data(),
                numHashFunctions,
                d_sequenceData2Bit,
                encodedSequencePitchInInts,
                numSequences,
                d_sequenceLengths,
                d_readIds,
                getKmerSize(),
                d_hashFunctionNumbers.data(),
                stream
            );

            #else

            GPUSequenceHasher<kmer_type> hasher;

            auto hashResult = hasher.hash(
                d_sequenceData2Bit,
                encodedSequencePitchInInts,
                numSequences,
                d_sequenceLengths,
                getKmerSize(),
                numHashFunctions,
                d_hashFunctionNumbers.data(),
                stream,
                mr
            );

            rmm::device_uvector<kmer_type> d_signatures_transposed(signaturesRowPitchElements * numSequences, stream, mr);
            helpers::call_transpose_kernel(
                d_signatures_transposed.data(),
                hashResult.d_hashvalues.data(),
                numSequences, 
                signaturesRowPitchElements, 
                signaturesRowPitchElements,
                stream
            );

            fixKeysForGpuHashTable<Key, Value>(d_signatures_transposed.data(), numSequences * numHashFunctions, stream);

            rmm::device_uvector<bool> d_isValid_transposed(
                numHashFunctions * numSequences,
                stream,
                mr
            );

            helpers::call_transpose_kernel(
                d_isValid_transposed.data(), 
                hashResult.d_isValid.data(),
                numSequences, 
                numHashFunctions, 
                numHashFunctions,
                stream
            );

            #if 1
            auto& h_insertTemp = h_insertTempMap[stream];
            auto& d_insertTemp = d_insertTempMap[stream];
            h_insertTemp.resize(numHashFunctions);
            d_insertTemp.resize(numHashFunctions);

            //create device views
            for(int i = 0; i < numHashFunctions; i++){
                new (&h_insertTemp[i]) GpuTable::DeviceTableInsertView(gpuHashTables[firstHashfunction + i]->makeDeviceInsertView());
            }
            //transfer to device
            CUDACHECK(cudaMemcpyAsync(
                d_insertTemp.data(),
                h_insertTemp.data(),
                sizeof(GpuTable::DeviceTableInsertView) * numHashFunctions,
                H2D,
                stream
            ));

            callInsertKernel(
                d_insertTemp.data(),
                numHashFunctions,
                d_signatures_transposed.data(),
                d_isValid_transposed.data(),
                numSequences,
                d_readIds,
                stream
            );

            #else

            for(int i = 0; i < numHashFunctions; i++){
                gpuHashTables[firstHashfunction + i]->insert(
                    d_signatures_transposed.data() + i * numSequences,
                    d_readIds,
                    numSequences,
                    stream
                );
            }

            #endif

            #endif
        }

        int checkInsertionErrors(
            int firstHashfunction,
            int numHashFunctions,
            cudaStream_t stream        
        ) override{
            int count = 0;
            for(int i = 0; i < numHashFunctions; i++){
                auto status = gpuHashTables[firstHashfunction + i]->pop_status(stream);
                CUDACHECK(cudaStreamSynchronize(stream));

                if(status.has_any_errors()){
                    count++;
                    std::cerr << "Error table " << (firstHashfunction + i) << " after insertion: " << status << "\n";
                }
            }
            return count;
        }

        MinhasherHandle makeMinhasherHandle() const override {
            auto data = std::make_unique<QueryData>(rmm::mr::get_current_device_resource());
            CUDACHECK(cudaGetDevice(&data->deviceId));

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
            
            {
                cub::SwitchDevice sd(tempdataVector[id]->deviceId);
                tempdataVector[id] = nullptr;
            }
            handle = constructHandle(std::numeric_limits<int>::max());
        }

        void compact(cudaStream_t stream) override {
            DeviceSwitcher ds(deviceId);

            std::size_t required_temp_bytes = 0;

            for(auto& table : gpuHashTables){
                std::size_t temp_bytes2 = 0;
                table->compact(nullptr, temp_bytes2, stream);
                required_temp_bytes = std::max(required_temp_bytes, temp_bytes2);
            }

            std::size_t freeMem, totalMem; 
            CUDACHECK(cudaMemGetInfo(&freeMem, &totalMem));

            void* temp = nullptr;
            if(required_temp_bytes < freeMem){
                CUDACHECK(cudaMalloc(&temp, required_temp_bytes));
            }else{
                CUDACHECK(cudaMallocManaged(&temp, required_temp_bytes));
                int deviceId = 0;
                CUDACHECK(cudaGetDevice(&deviceId));
                CUDACHECK(cudaMemAdvise(temp, required_temp_bytes, cudaMemAdviseSetAccessedBy, deviceId));
            }

            for(auto& table : gpuHashTables){
                table->compact(temp, required_temp_bytes, stream);
            }

            CUDACHECK(cudaFree(temp));
        }

        MemoryUsage getMemoryInfo() const noexcept override{
            MemoryUsage mem{};

            for(const auto& table : gpuHashTables){
                mem += table->getMemoryInfo();
            }

            return mem;
        }

        MemoryUsage getMemoryInfo(const MinhasherHandle& handle) const noexcept override{
            return getQueryDataFromHandle(handle)->getMemoryInfo();
        }

        int getNumResultsPerMapThreshold() const noexcept override{
            return resultsPerMapThreshold;
        }
        
        int getNumberOfMaps() const noexcept override{
            return gpuHashTables.size();
        }

        void destroy(){
            DeviceSwitcher sd(getDeviceId());
            gpuHashTables.clear();
        }

        bool hasGpuTables() const noexcept override {
            return true;
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

            std::size_t persistent_storage_bytes = 0;

            determineNumValues(
                nullptr,            
                persistent_storage_bytes,
                d_sequenceData2Bit,
                encodedSequencePitchInInts,
                d_sequenceLengths,
                numSequences,
                d_numValuesPerSequence,
                totalNumValues,
                stream,
                mr
            );

            queryData->d_singlepersistentbuffer.resize(persistent_storage_bytes, stream);

            determineNumValues(
                queryData->d_singlepersistentbuffer.data(),
                persistent_storage_bytes,
                d_sequenceData2Bit,
                encodedSequencePitchInInts,
                d_sequenceLengths,
                numSequences,
                d_numValuesPerSequence,
                totalNumValues,
                stream,
                mr
            );

            queryData->previousStage = QueryData::Stage::NumValues;
            queryData->d_numValuesPerSequence = d_numValuesPerSequence;
        }

        void retrieveValues(
            MinhasherHandle& queryHandle,
            int numSequences,
            int totalNumValues,
            read_number* d_values,
            const int* d_numValuesPerSequence,
            int* d_offsets, //numSequences + 1
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr
        ) const override {
            QueryData* const queryData = getQueryDataFromHandle(queryHandle);

            assert(queryData->previousStage == QueryData::Stage::NumValues);
            //STUPID INTERFACE EXPECTS d_numValuesPerSequence TO CONTAIN THE SAME VALUES AS RETURNED BY determineNumValues. This needs to be refactored. 
            //assert(queryData->d_numValuesPerSequence = d_numValuesPerSequence);

            std::size_t persistent_storage_bytes = queryData->d_singlepersistentbuffer.size();

            retrieveValues(
                queryData->d_singlepersistentbuffer.data(),
                persistent_storage_bytes,
                numSequences,
                totalNumValues,
                d_values,
                d_numValuesPerSequence,
                d_offsets, //numSequences + 1
                stream,
                mr
            );

            queryData->previousStage = QueryData::Stage::Retrieve;

        }


        //to make the following to functions private, their lambda kernels have to be replaced by normal kernels

        void determineNumValues(
            void* persistent_storage,            
            std::size_t& persistent_storage_bytes,
            const unsigned int* d_sequenceData2Bit,
            std::size_t encodedSequencePitchInInts,
            const int* d_sequenceLengths,
            int numSequences,
            int* d_numValuesPerSequence,
            int& totalNumValues,
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr
        ) const {

            const int numHashFunctions = gpuHashTables.size();
            const std::size_t signaturesRowPitchElements = numHashFunctions;

            void* persistent_allocations[3]{};
            std::size_t persistent_allocation_sizes[3]{};

            persistent_allocation_sizes[0] = sizeof(kmer_type) * numHashFunctions * numSequences; // d_sig_trans
            persistent_allocation_sizes[1] = sizeof(int) * numSequences * numHashFunctions; // d_numValuesPerSequencePerHash
            persistent_allocation_sizes[2] = sizeof(int) * numSequences * numHashFunctions; // d_numValuesPerSequencePerHashExclPSVert

            CUDACHECK(cub::AliasTemporaries(
                persistent_storage,
                persistent_storage_bytes,
                persistent_allocations,
                persistent_allocation_sizes
            ));

            if(persistent_storage == nullptr){
                return;
            }

            kmer_type* const d_signatures_transposed = static_cast<kmer_type*>(persistent_allocations[0]);
            int* const d_numValuesPerSequencePerHash = static_cast<int*>(persistent_allocations[1]);
            int* const d_numValuesPerSequencePerHashExclPSVert = static_cast<int*>(persistent_allocations[2]);

            rmm::device_uvector<int> d_hashFunctionNumbers(numHashFunctions, stream, mr);

            DeviceSwitcher ds(deviceId);

            CUDACHECK(cudaMemcpyAsync(
                d_hashFunctionNumbers.data(),
                h_currentHashFunctionNumbers.data(), 
                sizeof(int) * numHashFunctions, 
                H2D, 
                stream
            ));           

            GPUSequenceHasher<kmer_type> hasher;

            auto hashResult = hasher.hash(
                d_sequenceData2Bit,
                encodedSequencePitchInInts,
                numSequences,
                d_sequenceLengths,
                getKmerSize(),
                getNumberOfMaps(),
                d_hashFunctionNumbers.data(),
                stream,
                mr
            );

            helpers::call_transpose_kernel(
                d_signatures_transposed, 
                hashResult.d_hashvalues.data(),
                numSequences, 
                signaturesRowPitchElements, 
                signaturesRowPitchElements,
                stream
            );

            fixKeysForGpuHashTable<Key, Value>(d_signatures_transposed, numSequences * numHashFunctions, stream);

            //determine number of values per hashfunction per sequence
            #if 1

            {

            const int signaturesPitchInElements = numSequences;
            const int numValuesPerKeyPitchInElements = numSequences;
            constexpr int cgsize = GpuTable::DeviceTableView::cg_size();

            dim3 block(256, 1, 1);
            const int numBlocksPerTable = SDIV(numSequences, (block.x / cgsize));
            dim3 grid(numBlocksPerTable, std::min(65535, numHashFunctions), 1);

            gpuhashtablekernels::numValuesPerKeyCompactMultiTableKernel<<<grid, block, 0, stream>>>(
                d_deviceAccessibleTableViews.data(),
                numHashFunctions,
                resultsPerMapThreshold,
                d_signatures_transposed,
                signaturesPitchInElements,
                numSequences,
                d_numValuesPerSequencePerHash,
                numValuesPerKeyPitchInElements
            );

            }
            #else

            
            for(int i = 0; i < numHashFunctions; i++){
                gpuHashTables[i]->numValuesPerKeyCompact(
                    d_signatures_transposed + i * numSequences,
                    numSequences,
                    d_numValuesPerSequencePerHash + i * numSequences,
                    stream
                );
            }

            #endif

            // accumulate number of values per sequence in d_numValuesPerSequence
            // calculate vertical exclusive prefix sum
            helpers::lambda_kernel<<<SDIV(numSequences, 256), 256, 0, stream>>>(
                [=] __device__ (){
                    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
                    const int stride = blockDim.x * gridDim.x;

                    for(int i = tid; i < numSequences; i += stride){
                        d_numValuesPerSequencePerHashExclPSVert[0 * numSequences + i] = 0;
                    }

                    for(int i = tid; i < numSequences; i += stride){
                        int vertPS = 0;
                        for(int k = 0; k < numHashFunctions; k++){
                            const int num = d_numValuesPerSequencePerHash[k * numSequences + i];

                            vertPS += num;
                            if(k < numHashFunctions - 1){
                                d_numValuesPerSequencePerHashExclPSVert[(k+1) * numSequences + i] = vertPS;
                            }else{
                                d_numValuesPerSequence[i] = vertPS;
                            }
                        }
                    }
                }
            );

            rmm::device_scalar<int> d_totalNumValues(stream, mr);

            CubCallWrapper(mr).cubReduceSum(
                d_numValuesPerSequence, 
                d_totalNumValues.data(), 
                numSequences, 
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                &totalNumValues,
                d_totalNumValues.data(),
                sizeof(int),
                D2H,
                stream
            ));
        }

        void retrieveValues(
            void* persistentbufferFromNumValues,            
            std::size_t persistent_storage_bytes,
            int numSequences,
            int totalNumValues,
            read_number* d_values,
            const int* d_numValuesPerSequence,
            int* d_offsets, //numSequences + 1
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr
        ) const {
            if(totalNumValues == 0){
                CUDACHECK(cudaMemsetAsync(d_offsets, 0, sizeof(int) * (numSequences + 1), stream));
                return;
            }

            assert(persistentbufferFromNumValues != nullptr);

            const int numHashFunctions = gpuHashTables.size();

            void* persistent_allocations[3]{};
            std::size_t persistent_allocation_sizes[3]{};

            persistent_allocation_sizes[0] = sizeof(kmer_type) * numHashFunctions * numSequences; // d_sig_trans
            persistent_allocation_sizes[1] = sizeof(int) * numSequences * numHashFunctions; // d_numValuesPerSequencePerHash
            persistent_allocation_sizes[2] = sizeof(int) * numSequences * numHashFunctions; // d_numValuesPerSequencePerHashExclPSVert

            CUDACHECK(cub::AliasTemporaries(
                persistentbufferFromNumValues,
                persistent_storage_bytes,
                persistent_allocations,
                persistent_allocation_sizes
            ));

            DeviceSwitcher ds(deviceId);

            kmer_type* const d_signatures_transposed = static_cast<kmer_type*>(persistent_allocations[0]);
            int* const d_numValuesPerSequencePerHash = static_cast<int*>(persistent_allocations[1]);
            int* const d_numValuesPerSequencePerHashExclPSVert = static_cast<int*>(persistent_allocations[2]);

            rmm::device_uvector<int> d_queryOffsetsPerSequencePerHash(numSequences * numHashFunctions, stream, mr);


            //calculate global offsets for each sequence in output array
            CUDACHECK(cudaMemsetAsync(d_offsets, 0, sizeof(int), stream));

            CubCallWrapper(mr).cubInclusiveSum(
                d_numValuesPerSequence,
                d_offsets + 1,
                numSequences,
                stream
            );

            // compute destination offsets for each hashtable such that values of different tables 
            // for the same sequence are stored contiguous in the result array

            helpers::lambda_kernel<<<SDIV(numSequences, 256), 256, 0, stream>>>(
                [
                    d_queryOffsetsPerSequencePerHash = d_queryOffsetsPerSequencePerHash.data(),
                    d_numValuesPerSequencePerHashExclPSVert,
                    numSequences,
                    numHashFunctions,
                    d_offsets
                ] __device__ (){
                    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
                    const int stride = blockDim.x * gridDim.x;

                    for(int i = tid; i < numSequences; i += stride){
                        
                        const int base = d_offsets[i];

                        //k == 0 is a copy from d_offsets
                        d_queryOffsetsPerSequencePerHash[0 * numSequences + i] = base;

                        for(int k = 1; k < numHashFunctions; k++){
                            d_queryOffsetsPerSequencePerHash[k * numSequences + i] = base + d_numValuesPerSequencePerHashExclPSVert[k * numSequences + i];
                        }
                    }
                }
            );

            //retrieve values

            #if 1
            {
            const int signaturesPitchInElements = numSequences;
            const int numValuesPerKeyPitchInElements = numSequences;
            const int beginOffsetsPitchInElements = numSequences;
            constexpr int cgsize = GpuTable::DeviceTableView::cg_size();

            dim3 block(256, 1, 1);
            const int numBlocksPerTable = SDIV(numSequences, (block.x / cgsize));
            dim3 grid(numBlocksPerTable, std::min(65535, numHashFunctions), 1);

            gpuhashtablekernels::retrieveCompactKernel<<<grid, block, 0, stream>>>(
                d_deviceAccessibleTableViews.data(),
                numHashFunctions,
                d_signatures_transposed,
                signaturesPitchInElements,
                d_queryOffsetsPerSequencePerHash.data(),
                beginOffsetsPitchInElements,
                d_numValuesPerSequencePerHash,
                numValuesPerKeyPitchInElements,
                resultsPerMapThreshold,
                numSequences,
                d_values
            );
            }
            #else

            for(int i = 0; i < numHashFunctions; i++){
                gpuHashTables[i]->retrieveCompact(
                    d_signatures_transposed + i * numSequences,
                    d_queryOffsetsPerSequencePerHash.data()  + i * numSequences,
                    d_numValuesPerSequencePerHash + i * numSequences,
                    numSequences,
                    d_values,
                    stream
                );
            }

            #endif
        }

        int getKmerSize() const noexcept override{
            return kmerSize;
        }


        constexpr int getDeviceId() const noexcept{
            return deviceId;
        }

        void setThreadPool(ThreadPool* /*tp*/) override {}

        void setHostMemoryLimitForConstruction(std::size_t /*bytes*/) override{

        }

        void setDeviceMemoryLimitsForConstruction(const std::vector<std::size_t>&) override {

        }

        void constructionIsFinished(cudaStream_t stream) override {
            auto numberOfAvailableHashFunctions = h_currentHashFunctionNumbers.size();
            std::vector<GpuTable::DeviceTableView> views;
            for(const auto& ptr : gpuHashTables){
                views.emplace_back(ptr->makeDeviceView());
            }

            d_deviceAccessibleTableViews.resize(numberOfAvailableHashFunctions);
            CUDACHECK(cudaMemcpyAsync(
                d_deviceAccessibleTableViews.data(),
                views.data(),
                sizeof(GpuTable::DeviceTableView) * numberOfAvailableHashFunctions,
                H2D,
                stream
            ));

            CUDACHECK(cudaStreamSynchronize(stream));

            h_insertTempMap.clear();
            d_insertTempMap.clear();
        }

        void writeToStream(std::ostream& /*os*/) const override{
            std::cerr << "SingleGpuMinhasher::writeToStream not supported\n";
        }

        int loadFromStream(std::ifstream& /*is*/, int /*numMapsUpperLimit*/) override{
            std::cerr << "SingleGpuMinhasher::loadFromStream not supported\n";
            return 0;
        } 

        bool canWriteToStream() const noexcept override { return false; };
        bool canLoadFromStream() const noexcept override { return false; };

private:

        void finalize(cudaStream_t stream = 0){
            compact(stream);
        }

        std::uint64_t getKmerMask() const{
            constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;

            return std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - getKmerSize()) * 2);
        }

        constexpr float getLoad() const noexcept{
            return loadfactor;
        }

        QueryData* getQueryDataFromHandle(const MinhasherHandle& queryHandle) const{
            std::shared_lock<SharedMutex> lock(sharedmutex);

            return tempdataVector[queryHandle.getId()].get();
        }

        void callInsertKernel(
            GpuTable::DeviceTableInsertView* d_insertViews,
            int numHashFunctions,
            const kmer_type* d_signatures_transposed,
            const bool* d_isValid_transposed,
            int numSequences,
            const read_number* d_readIds,
            cudaStream_t stream
        ){
            auto insertkernel = sgpuminhasherkernels::insertIntoTablesKernel<GpuTable::DeviceTableInsertView, Key, Value>;

            constexpr int blocksize = 512;
            int deviceId = 0;
            int numSMs = 0;
            int maxBlocksPerSM = 0;
            CUDACHECK(cudaGetDevice(&deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
            CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                &maxBlocksPerSM,
                insertkernel,
                blocksize, 
                0
            ));
        
            const int maxBlocks = maxBlocksPerSM * numSMs;  
        
            dim3 block(blocksize, 1, 1);

            constexpr int groupsize = GpuTable::DeviceTableInsertView::cg_size();
            constexpr int warpsize = 32;
            constexpr int groupsPerWarp = warpsize / groupsize;
            const int warpsPerTable = SDIV(numSequences, groupsPerWarp);

            const int numBlocks = SDIV(numHashFunctions * warpsPerTable, blocksize / warpsize);
            dim3 grid(std::min(numBlocks, maxBlocks), 1, 1);

            insertkernel<<<grid, block, 0, stream>>>(
                d_insertViews,
                numHashFunctions,
                d_signatures_transposed,
                d_isValid_transposed,
                numSequences,
                d_readIds
            );
            CUDACHECKASYNC;
        }

        void callHashAndInsertKernel(
            GpuTable::DeviceTableInsertView* d_insertViews,
            const int numTables,
            const unsigned int* d_sequences2Bit,
            std::size_t encodedSequencePitchInInts,
            const int numSequences,
            const int* d_sequenceLengths,
            const read_number* d_readIds,
            int k,
            const int* d_hashFunctionNumbers,
            cudaStream_t stream
        ){
            auto hashAndInsertkernel = sgpuminhasherkernels::hashAndFixAndInsertIntoTablesKernel<GpuTable::DeviceTableInsertView, Key>;

            constexpr int blocksize = 512;
            int deviceId = 0;
            int numSMs = 0;
            int maxBlocksPerSM = 0;
            CUDACHECK(cudaGetDevice(&deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
            CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                &maxBlocksPerSM,
                hashAndInsertkernel,
                blocksize, 
                0
            ));
        
            const int maxBlocks = maxBlocksPerSM * numSMs;  
        
            constexpr int groupsize = GpuTable::DeviceTableInsertView::cg_size();
            constexpr int warpsize = 32;
            constexpr int groupsPerWarp = warpsize / groupsize;
            const int warpsPerTable = SDIV(numSequences, groupsPerWarp);

            const int numBlocks = SDIV(numTables * warpsPerTable, blocksize / warpsize);
            dim3 block(blocksize, 1, 1);
            dim3 grid(std::min(numBlocks, maxBlocks), 1, 1);
            // dim3 block(32,1,1);
            // dim3 grid(1, 1, 1);

            hashAndInsertkernel<<<grid, block, 0, stream>>>(
                d_insertViews,
                numTables,
                d_sequences2Bit,
                encodedSequencePitchInInts,
                numSequences,
                d_sequenceLengths,
                d_readIds,
                k,
                d_hashFunctionNumbers
            );
            CUDACHECKASYNC;
        }

        mutable int counter = 0;
        mutable SharedMutex sharedmutex{};
        //mutable std::shared_mutex sharedmutex{};

        int deviceId{};
        int maxNumKeys{};
        int kmerSize{};
        int resultsPerMapThreshold{};
        float loadfactor{};
        HostBuffer<int> h_currentHashFunctionNumbers{};
        std::vector<std::unique_ptr<GpuTable>> gpuHashTables{};
        std::map<cudaStream_t, helpers::SimpleAllocationPinnedHost<GpuTable::DeviceTableInsertView>> h_insertTempMap{};
        std::map<cudaStream_t, helpers::SimpleAllocationDevice<GpuTable::DeviceTableInsertView>> d_insertTempMap{};
        helpers::SimpleAllocationDevice<GpuTable::DeviceTableView, 0> d_deviceAccessibleTableViews{};
        mutable std::vector<std::unique_ptr<QueryData>> tempdataVector{};
    };


}
}




#endif

#endif //#ifdef CARE_HAS_WARPCORE