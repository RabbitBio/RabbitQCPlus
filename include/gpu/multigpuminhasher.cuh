#ifdef CARE_HAS_WARPCORE


#ifndef CARE_MULTI_GPU_MINHASHER_CUH
#define CARE_MULTI_GPU_MINHASHER_CUH

#include <config.hpp>


#include <gpu/gpureadstorage.cuh>
#include <gpu/cuda_unique.cuh>
#include <gpu/singlegpuminhasher.cuh>
#include <gpu/gpuminhasher.cuh>
#include <gpu/cudaerrorcheck.cuh>
#include <gpu/cubwrappers.cuh>

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
#include <cassert>

#include <rmm/mr/device/device_memory_resource.hpp>
#include <rmm/mr/device/per_device_resource.hpp>
#include <rmm/device_uvector.hpp>
#include <rmm/device_scalar.hpp>
#include <gpu/rmm_utilities.cuh>
namespace care{
namespace gpu{

    namespace multigpuminhasherkernels{
        template<int blocksize, int itemsPerThread>
        __global__
        void aggregatePartitionResultsSingleBlockKernel(
            const int* __restrict__ numResultsPerSequencePerPartition,
            int* __restrict__ numResultsPerSequence,
            int numSequences,
            int numPartitions,
            int* __restrict__ maxNumResultsPerSequence,
            int* __restrict__ offsets
        ){
            assert(gridDim.x * gridDim.y * gridDim.z == 1);

            struct BlockPrefixCallbackOp{
                // Running prefix
                int running_total;
                // Constructor
                __device__ BlockPrefixCallbackOp(int running_total) : running_total(running_total) {}
                // Callback operator to be entered by the first warp of threads in the block.
                // Thread-0 is responsible for returning a value for seeding the block-wide scan.
                __device__ int operator()(int block_aggregate)
                {
                    int old_prefix = running_total;
                    running_total += block_aggregate;
                    return old_prefix;
                }
            };

            using BlockReduce = cub::BlockReduce<int, blocksize>;
            using BlockScan = cub::BlockScan<int, blocksize>;

            __shared__ typename BlockReduce::TempStorage smem_reduce;
            __shared__ typename BlockScan::TempStorage smem_scan;

            constexpr int itemsPerIteration = blocksize * itemsPerThread;

            const int numIterations = SDIV(numSequences, itemsPerIteration);
            int myMax = 0;
            BlockPrefixCallbackOp prefix_op(0);

            //full iterations
            for(int iteration = 0; iteration < numIterations; iteration++){
                const int s = iteration * itemsPerIteration + threadIdx.x;
                int sum = 0;

                if(s < numSequences){
                    for(int r = 0; r < numPartitions; r++){
                        sum += numResultsPerSequencePerPartition[r * numSequences + s];
                    }

                    numResultsPerSequence[s] = sum;
                    myMax = max(myMax, sum);
                }

                BlockScan(smem_scan).ExclusiveSum(sum, sum, prefix_op);

                if(s < numSequences){
                    offsets[s] = sum;
                }
            }

            myMax = BlockReduce(smem_reduce).Reduce(myMax, cub::Max{});

            if(threadIdx.x == 0){
                *maxNumResultsPerSequence = myMax;
                offsets[numSequences] = prefix_op.running_total;
            }
        };

        template<class Dummy = void>
        __global__
        void aggregateNumValuesPartitionResultsKernel(
            const int* __restrict__ numResultsPerSequencePerPartition,
            int* __restrict__ numResultsPerSequence,
            int numSequences,
            int numPartitions
        ){
            const int tid = threadIdx.x + blockIdx.x * blockDim.x;
            const int stride = blockDim.x * gridDim.x;

            for(int s = tid; s < numSequences; s += stride){
                int sum = 0;
                
                for(int r = 0; r < numPartitions; r++){
                    sum += numResultsPerSequencePerPartition[r * numSequences + s];
                }

                numResultsPerSequence[s] = sum;
            }

        };
    
        template<class Dummy = void>
        __global__
        void copyToInterleavedKernel(
            const read_number* __restrict__ inputdata,
            const int* __restrict__ inputoffsets,
            const int* __restrict__ inputsegmentsizes,
            const int* __restrict__ inputnumpergpuPS,
            const int* __restrict__ outputbeginoffsets,
            read_number* __restrict__ outputdata,
            int numSequences,
            int partitions
        ){                
            for(int i = blockIdx.x; i < numSequences; i += gridDim.x){
                const int beginoffset = outputbeginoffsets[i];

                int runningOffset = 0;

                for(int r = 0; r < partitions; r++){
                    const int segmentsize = inputsegmentsizes[r * numSequences + i];
                    const int inputOffset = inputoffsets[r * (numSequences+1) + i];
                    const int gpuOffset = inputnumpergpuPS[r];

                    const read_number* myinput = inputdata + gpuOffset + inputOffset;

                    for(int k = threadIdx.x; k < segmentsize; k += blockDim.x){
                        outputdata[beginoffset + runningOffset + k] = myinput[k];
                    }

                    runningOffset += segmentsize;
                }
            }
        }

        template<class Dummy = void>
        __global__ 
        void copyResultsToDestinationKernel(
            const read_number* input,
            int* inputoffsets,
            read_number* output,
            const int* outputoffsets,
            const int* d_numValuesPerSequence,
            int numSequences
        ){
            for(int s = blockIdx.x; s < numSequences; s += gridDim.x){
                const int inputoffset = inputoffsets[s];
                const int outputoffset = outputoffsets[s];
                const int segmentsize = d_numValuesPerSequence[s];

                for(int i = threadIdx.x; i < segmentsize; i += blockDim.x){
                    output[outputoffset + i] = input[inputoffset + i];
                }

                __syncthreads();
                //update d_offsets within this kernel to avoid additional api call
                if(threadIdx.x == 0){
                    inputoffsets[s] = outputoffset;
                }
            }

            //last entry of offsets (total number) is not used for copying. no need for sync
            if(blockIdx.x == 0 && threadIdx.x == 0){
                inputoffsets[numSequences] = outputoffsets[numSequences];
            }

        }

    }


    class MultiGpuMinhasher : public GpuMinhasher{
    public:
        using Key = GpuMinhasher::Key;

        enum class Layout {FirstFit, EvenShare};
    private:
        using DeviceSwitcher = cub::SwitchDevice;

        template<class T>
        using HostBuffer = helpers::SimpleAllocationPinnedHost<T, 5>;
        template<class T>
        using DeviceBuffer = helpers::SimpleAllocationDevice<T, 5>;
        //using DeviceBuffer = helpers::SimpleAllocationPinnedHost<T, 5>;

        struct QueryData{
            enum class Stage{
                None,
                NumValues,
                Retrieve
            };

            Stage previousStage = Stage::None;

            CudaEvent callerEvent{cudaEventDisableTiming};
            HostBuffer<int> pinnedData{};
            std::vector<CudaStream> streams{};
            std::vector<CudaEvent> events{};
            std::vector<std::unique_ptr<MinhasherHandle>> singlegpuMinhasherHandles;

            std::vector<rmm::device_uvector<int>> vec_d_numValuesPerSequence;

            MemoryUsage getMemoryInfo() const{
                MemoryUsage mem{};

                return mem;
            }
        };

    public: 

        MultiGpuMinhasher(Layout layout_, int maxNumKeys_, int maxValuesPerKey, int k, float loadfactor_, std::vector<int> deviceIds_)
            : layout(layout_), maxNumKeys(maxNumKeys_), kmerSize(k), resultsPerMapThreshold(maxValuesPerKey), loadfactor(loadfactor_), deviceIds(deviceIds_)
        {
            for(auto deviceId : deviceIds){
                cub::SwitchDevice sd{deviceId};
                auto mh = std::make_unique<SingleGpuMinhasher>(maxNumKeys, resultsPerMapThreshold, k, loadfactor_);
                sgpuMinhashers.emplace_back(std::move(mh));

                hashFunctionIdsPerGpu.emplace_back();
            }
        }

        int addHashTables(int numAdditionalTables, const int* hashFunctionIds, cudaStream_t stream) override{

            std::vector<int> hashFunctionIdsTmp(hashFunctionIds, hashFunctionIds + numAdditionalTables);

            CudaEvent event;
            event.record(stream);

            const int numDevices = deviceIds.size();
            std::vector<CudaEvent> events;
            for(int g = 0; g < numDevices; g++){
                cub::SwitchDevice sd{deviceIds[g]};
                events.emplace_back(cudaEventDisableTiming);

                CUDACHECK(cudaStreamWaitEvent(cudaStreamPerThread, event, 0));
            }

            int remainingTables = numAdditionalTables;

            if(layout == Layout::EvenShare){
                while(remainingTables > 0){
                    int numZeroAdded = 0;

                    for(int g = 0; g < numDevices; g++){
                        if(remainingTables > 0){
                            cub::SwitchDevice sd{deviceIds[g]};
                            int addedTables = sgpuMinhashers[g]->addHashTables(1, hashFunctionIdsTmp.data(), cudaStreamPerThread);

                            for(int x = 0; x < addedTables; x++){
                                hashTableLocations.push_back(g);
                            }
                            hashFunctionIdsPerGpu[g].insert(hashFunctionIdsPerGpu[g].end(), hashFunctionIdsTmp.begin(), hashFunctionIdsTmp.begin() + addedTables);

                            hashFunctionIdsTmp.erase(hashFunctionIdsTmp.begin(), hashFunctionIdsTmp.begin() + addedTables);
                            remainingTables -= addedTables;

                            if(addedTables == 0){
                                numZeroAdded++;
                            }
                        }
                    }

                    if(numZeroAdded == numDevices){
                        break;
                    }
                }
            }else{
                assert(layout == Layout::FirstFit);

                for(int g = 0; g < numDevices; g++){
                    if(remainingTables > 0){
                        cub::SwitchDevice sd{deviceIds[g]};
                        int addedTables = sgpuMinhashers[g]->addHashTables(remainingTables, hashFunctionIdsTmp.data(), cudaStreamPerThread);

                        for(int x = 0; x < addedTables; x++){
                            hashTableLocations.push_back(g);
                        }
                        hashFunctionIdsPerGpu[g].insert(hashFunctionIdsPerGpu[g].end(), hashFunctionIdsTmp.begin(), hashFunctionIdsTmp.begin() + addedTables);

                        hashFunctionIdsTmp.erase(hashFunctionIdsTmp.begin(), hashFunctionIdsTmp.begin() + addedTables);
                        remainingTables -= addedTables;
                    }
                }
            }

            for(int g = 0; g < numDevices; g++){
                cub::SwitchDevice sd{deviceIds[g]};
                events[g].record(cudaStreamPerThread);
            }

            for(int g = 0; g < numDevices; g++){
                CUDACHECK(cudaStreamWaitEvent(stream, events[g], 0));
            }

            return numAdditionalTables - remainingTables;
        }

        void insert(
            const unsigned int* d_sequenceData2Bit,
            int numSequences,
            const int* d_sequenceLengths,
            std::size_t encodedSequencePitchInInts,
            const read_number* d_readIds,
            const read_number* h_readIds,
            int firstHashfunction,
            int numHashfunctions,
            const int* h_hashFunctionIds,
            cudaStream_t stream,
            rmm::mr::device_memory_resource* /*mr*/
        ) override {
            assert(firstHashfunction + numHashfunctions <= hashTableLocations.size());
            if(numHashfunctions == 0) return;
            if(numSequences == 0) return;

            int oldDeviceId = 0;
            CUDACHECK(cudaGetDevice(&oldDeviceId));

            CudaEvent event;
            event.record(stream);
            std::vector<CudaEvent> events;

            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};
                events.emplace_back(cudaEventDisableTiming);

                CUDACHECK(cudaStreamWaitEvent(cudaStreamPerThread, event, 0));
            }

            std::vector<int> numHashfunctionsPerTargetGpu(deviceIds.size(), 0);
            std::vector<int> firstHashfunctionPerTargetGpu(deviceIds.size(), 0);
            std::vector<std::vector<int>> hashFunctionIdsPerTargetGpu(deviceIds.size());

            for(int g = 0; g < int(deviceIds.size()); g++){
                int countBefore = 0;
                for(int i = 0; i < firstHashfunction + numHashfunctions; i++){
                    if(hashTableLocations[i] == g){
                        if(i < firstHashfunction){
                            countBefore++;
                        }else{
                            numHashfunctionsPerTargetGpu[g]++;
                            hashFunctionIdsPerTargetGpu[g].push_back(h_hashFunctionIds[i - firstHashfunction]);
                        }
                    }
                }
                firstHashfunctionPerTargetGpu[g] = countBefore;
                
                assert(numHashfunctionsPerTargetGpu[g] == int(hashFunctionIdsPerTargetGpu[g].size()));
            }
            
            std::vector<rmm::device_uvector<unsigned int>> vec_d_sequenceData2Bit_target;
            std::vector<rmm::device_uvector<int>> vec_d_sequenceLengths_target;
            std::vector<rmm::device_uvector<read_number>> vec_d_readIds_target;

            //broadcast to all gpus, excluding current gpu
            for(int g = 0; g < int(deviceIds.size()); g++){
                if(numHashfunctionsPerTargetGpu[g] > 0){
                    if(deviceIds[g] != oldDeviceId){
                        cub::SwitchDevice sd{deviceIds[g]};

                        //copy input data to target gpu
                        auto* targetmr = rmm::mr::get_current_device_resource();
                        rmm::device_uvector<unsigned int> d_sequenceData2Bit_target(encodedSequencePitchInInts * numSequences, cudaStreamPerThread, targetmr);
                        rmm::device_uvector<int> d_sequenceLengths_target(numSequences, cudaStreamPerThread, targetmr);
                        rmm::device_uvector<read_number> d_readIds_target(numSequences, cudaStreamPerThread, targetmr);

                        CUDACHECK(cudaMemcpyPeerAsync(
                            d_sequenceData2Bit_target.data(),
                            deviceIds[g],
                            d_sequenceData2Bit,
                            oldDeviceId,
                            sizeof(unsigned int) * encodedSequencePitchInInts * numSequences,
                            cudaStreamPerThread
                        ));

                        CUDACHECK(cudaMemcpyPeerAsync(
                            d_sequenceLengths_target.data(),
                            deviceIds[g],
                            d_sequenceLengths,
                            oldDeviceId,
                            sizeof(int) * numSequences,
                            cudaStreamPerThread
                        ));

                        CUDACHECK(cudaMemcpyPeerAsync(
                            d_readIds_target.data(),
                            deviceIds[g],
                            d_readIds,
                            oldDeviceId,
                            sizeof(read_number) * numSequences,
                            cudaStreamPerThread
                        ));

                        vec_d_sequenceData2Bit_target.push_back(std::move(d_sequenceData2Bit_target));
                        vec_d_sequenceLengths_target.push_back(std::move(d_sequenceLengths_target));
                        vec_d_readIds_target.push_back(std::move(d_readIds_target));
                    }else{
                        vec_d_sequenceData2Bit_target.emplace_back(0, cudaStreamPerThread);
                        vec_d_sequenceLengths_target.emplace_back(0, cudaStreamPerThread);
                        vec_d_readIds_target.emplace_back(0, cudaStreamPerThread);
                    }
                }
            }

            //insert on each gpu
            for(int g = 0; g < int(deviceIds.size()); g++){
                if(numHashfunctionsPerTargetGpu[g] > 0){
                    cub::SwitchDevice sd{deviceIds[g]};

                    const unsigned int* d_seq = d_sequenceData2Bit;
                    const int* d_len = d_sequenceLengths;
                    const read_number* d_ids = d_readIds;

                    if(deviceIds[g] != oldDeviceId){
                        d_seq = vec_d_sequenceData2Bit_target[g].data();
                        d_len = vec_d_sequenceLengths_target[g].data();
                        d_ids = vec_d_readIds_target[g].data();
                    }

                    sgpuMinhashers[g]->insert(
                        d_seq,
                        numSequences,
                        d_len,
                        encodedSequencePitchInInts,
                        d_ids,
                        h_readIds,
                        firstHashfunctionPerTargetGpu[g],
                        numHashfunctionsPerTargetGpu[g],
                        hashFunctionIdsPerTargetGpu[g].data(),
                        cudaStreamPerThread,
                        rmm::mr::get_current_device_resource()
                    );
                    if(deviceIds[g] != oldDeviceId){
                        vec_d_sequenceData2Bit_target[g].release();
                        vec_d_sequenceLengths_target[g].release();
                        vec_d_readIds_target[g].release();
                    }
                    events[g].record(cudaStreamPerThread);
                }
            }

            for(int g = 0; g < int(deviceIds.size()); g++){
                if(numHashfunctionsPerTargetGpu[g] > 0){
                    CUDACHECK(cudaStreamWaitEvent(stream, events[g], 0));
                }
            }
        }

        int checkInsertionErrors(
            int firstHashfunction,
            int numHashfunctions,
            cudaStream_t stream        
        ) override{
            CudaEvent event;
            event.record(stream);

            std::vector<int> numHashfunctionsPerTargetGpu(deviceIds.size(), 0);
            std::vector<int> firstHashfunctionPerTargetGpu(deviceIds.size(), 0);

            for(int g = 0; g < int(deviceIds.size()); g++){
                int countBefore = 0;
                for(int i = 0; i < firstHashfunction + numHashfunctions; i++){
                    if(hashTableLocations[i] == g){
                        if(i < firstHashfunction){
                            countBefore++;
                        }else{
                            numHashfunctionsPerTargetGpu[g]++;
                        }
                    }
                }
                firstHashfunctionPerTargetGpu[g] = countBefore;
            }

            int count = 0;

            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};

                CUDACHECK(cudaStreamWaitEvent(cudaStreamPerThread, event, 0));
                count += sgpuMinhashers[g]->checkInsertionErrors(
                    firstHashfunctionPerTargetGpu[g],
                    numHashfunctionsPerTargetGpu[g],
                    cudaStreamPerThread
                );
            }

            return count;
        }

        // bool tryReplication(){
        //     if(sgpuMinhashers.size() == 1 && usableDeviceIds.size() < deviceIds.size()){
        //         //all hashtables fit into one gpu. try to replace the hash tables on all gpus

        //         std::vector<std::unique_ptr<SingleGpuMinhasher>> replicas{};
        //         bool ok = false;
        //         try{
        //             nvtx::push_range("replicate single gpu minhasher", 0);

        //             for(std::size_t i = 1; i < deviceIds.size(); i++){
        //                 const int targetDeviceId = deviceIds[i];
        //                 helpers::CpuTimer rtimer("make singlegpu minhasher replica");
        //                 replicas.emplace_back(sgpuMinhashers[0]->makeCopy(targetDeviceId));
        //                 rtimer.print();
        //             }
        //             ok = std::all_of(replicas.begin(), replicas.end(), [](const auto& uniqueptr){ return bool(uniqueptr); });

        //             nvtx::pop_range();
        //         }catch(...){
        //             cudaGetLastError();
        //             std::cerr << "error replicating single gpu minhasher. Skipping.\n";
        //         }
        //         if(ok){                    
        //             sgpuMinhashers.insert(sgpuMinhashers.end(), std::make_move_iterator(replicas.begin()), std::make_move_iterator(replicas.end()));

        //             HostBuffer<int> h_currentHashFunctionNumbers(vec_h_currentHashFunctionIds[0].size());
        //             std::copy(vec_h_currentHashFunctionIds[0].begin(), vec_h_currentHashFunctionIds[0].end(), h_currentHashFunctionNumbers.begin());
        //             vec_h_currentHashFunctionIds.push_back(std::move(h_currentHashFunctionNumbers));

        //             usableDeviceIds = deviceIds;

        //             isReplicatedSingleGpu = true;
        //         }

        //         return ok;
        //     }else{
        //         return false;
        //     }
        // }

        MinhasherHandle makeMinhasherHandle() const override{
            auto ptr = std::make_unique<QueryData>();

            const int numMinhashers = sgpuMinhashers.size();

            ptr->streams.resize(numMinhashers);
            ptr->events.resize(numMinhashers);

            for(int i = 0; i < numMinhashers; i++){
                DeviceSwitcher ds(sgpuMinhashers[i]->getDeviceId());

                ptr->streams[i] = std::move(CudaStream{});
                ptr->events[i] = std::move(CudaEvent{cudaEventDisableTiming});

                ptr->singlegpuMinhasherHandles.emplace_back(std::make_unique<MinhasherHandle>(sgpuMinhashers[i]->makeMinhasherHandle()));
            }

            ptr->pinnedData.resize(2 * numMinhashers);

            CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));

            std::unique_lock<SharedMutex> lock(sharedmutex);
            const int handleid = counter++;
            MinhasherHandle h = constructHandle(handleid);

            tempdataVector.emplace_back(std::move(ptr));

            return h;
        }

        void destroyHandle(MinhasherHandle& handle) const override{            

            std::unique_lock<SharedMutex> lock(sharedmutex);

            const int id = handle.getId();
            assert(id < int(tempdataVector.size()));

            const int numMinhashers = sgpuMinhashers.size();
            for(int i = 0; i < numMinhashers; i++){
                sgpuMinhashers[i]->destroyHandle(*tempdataVector[id]->singlegpuMinhasherHandles[i]);
            }
            
            {
                tempdataVector[id] = nullptr;
            }
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
        ) const override{
            QueryData* const queryData = getQueryDataFromHandle(queryHandle);
            queryData->previousStage = QueryData::Stage::NumValues;

            if(numSequences == 0){
                // CUDACHECK(cudaMemsetAsync(d_numValuesPerSequence, 0, sizeof(int) * numSequences, stream));
                totalNumValues = 0;
                return;
            }

            int oldDeviceId = 0;
            CUDACHECK(cudaGetDevice(&oldDeviceId));
            
            rmm::device_uvector<int> d_numValuesPerSequencePerGpu(numSequences * deviceIds.size(), stream, mr);

            CUDACHECK(cudaEventRecord(queryData->callerEvent, stream));

            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};

                CUDACHECK(cudaStreamWaitEvent(queryData->streams[g], queryData->callerEvent, 0));
            }

            std::vector<rmm::device_uvector<unsigned int>> vec_d_sequenceData2Bit_target;
            std::vector<rmm::device_uvector<int>> vec_d_sequenceLengths_target;

            //broadcast to other gpus
            for(int g = 0; g < int(deviceIds.size()); g++){
                if(deviceIds[g] != oldDeviceId){
                    cub::SwitchDevice sd{deviceIds[g]};

                    auto* targetmr = rmm::mr::get_current_device_resource();
                    rmm::device_uvector<unsigned int> d_sequenceData2Bit_target(encodedSequencePitchInInts * numSequences, queryData->streams[g].getStream(), targetmr);
                    rmm::device_uvector<int> d_sequenceLengths_target(numSequences, queryData->streams[g].getStream(), targetmr);

                    CUDACHECK(cudaMemcpyPeerAsync(
                        d_sequenceData2Bit_target.data(),
                        deviceIds[g],
                        d_sequenceData2Bit,
                        oldDeviceId,
                        sizeof(unsigned int) * encodedSequencePitchInInts * numSequences,
                        queryData->streams[g]
                    ));

                    CUDACHECK(cudaMemcpyPeerAsync(
                        d_sequenceLengths_target.data(),
                        deviceIds[g],
                        d_sequenceLengths,
                        oldDeviceId,
                        sizeof(int) * numSequences,
                        queryData->streams[g]
                    ));

                    vec_d_sequenceData2Bit_target.push_back(std::move(d_sequenceData2Bit_target));
                    vec_d_sequenceLengths_target.push_back(std::move(d_sequenceLengths_target));
                }else{
                    vec_d_sequenceData2Bit_target.emplace_back(0, queryData->streams[g].getStream());
                    vec_d_sequenceLengths_target.emplace_back(0, queryData->streams[g].getStream());
                }
            }

            //determine num values on each gpu, and collect results in d_numValuesPerSequencePerGpu
            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};

                const unsigned int* d_seq = d_sequenceData2Bit;
                const int* d_len = d_sequenceLengths;
                if(deviceIds[g] != oldDeviceId){
                    d_seq = vec_d_sequenceData2Bit_target[g].data();
                    d_len = vec_d_sequenceLengths_target[g].data();
                }

                int& myTotalNumValues = queryData->pinnedData[g];
                auto* targetmr = rmm::mr::get_current_device_resource();
                rmm::device_uvector<int> d_numValuesPerSequence_target(numSequences, queryData->streams[g].getStream(), targetmr);

                sgpuMinhashers[g]->determineNumValues(
                    *queryData->singlegpuMinhasherHandles[g].get(),
                    d_seq,
                    encodedSequencePitchInInts,
                    d_len,
                    numSequences,
                    d_numValuesPerSequence_target.data(),
                    myTotalNumValues,
                    queryData->streams[g],
                    targetmr
                );

                if(deviceIds[g] != oldDeviceId){
                    vec_d_sequenceData2Bit_target[g].release();
                    vec_d_sequenceLengths_target[g].release();
                }

                queryData->vec_d_numValuesPerSequence.push_back(std::move(d_numValuesPerSequence_target));
            }

            //gather num values to d_numValuesPerSequencePerGpu
            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};

                CUDACHECK(cudaMemcpyPeerAsync(
                    d_numValuesPerSequencePerGpu.data() + numSequences * g,
                    oldDeviceId,
                    queryData->vec_d_numValuesPerSequence[g].data(),
                    deviceIds[g],
                    sizeof(int) * numSequences,
                    queryData->streams[g]
                ));
            }

            //join streams to wait for pinnedData and d_numValuesPerSequencePerGpu
            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};
                CUDACHECK(cudaStreamSynchronize(queryData->streams[g])); 
            }

            dim3 block = 128;
            dim3 grid = SDIV(numSequences, block.x);
            multigpuminhasherkernels::aggregateNumValuesPartitionResultsKernel
                    <<<grid, block, 0, stream>>>(
                d_numValuesPerSequencePerGpu.data(),
                d_numValuesPerSequence,
                numSequences,
                deviceIds.size()
            );
            CUDACHECKASYNC;

            totalNumValues = 0;
            for(int g = 0; g < int(deviceIds.size()); g++){
                totalNumValues += queryData->pinnedData[g];
            }
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
        ) const override{
            QueryData* const queryData = getQueryDataFromHandle(queryHandle);
            assert(queryData->previousStage == QueryData::Stage::NumValues);
            queryData->previousStage = QueryData::Stage::Retrieve;

            if(numSequences == 0){
                cudaMemsetAsync(d_offsets, 0, sizeof(int) * (numSequences + 1), stream);
                return;
            }

            int oldDeviceId = 0;
            CUDACHECK(cudaGetDevice(&oldDeviceId));

            rmm::device_uvector<read_number> d_allValues(totalNumValues, stream, mr);
            rmm::device_uvector<int> d_numValuesPerSequencePerGpu(numSequences * deviceIds.size(), stream, mr);
            rmm::device_uvector<int> d_offsetsPerSequencePerGpu((numSequences+1) * deviceIds.size(), stream, mr);

            CUDACHECK(cudaEventRecord(queryData->callerEvent, stream));

            CubCallWrapper(mr).cubInclusiveSum(
                d_numValuesPerSequence,
                d_offsets + 1,
                numSequences,
                stream
            );
            CUDACHECK(cudaMemsetAsync(d_offsets, 0, sizeof(int), stream));

            std::vector<rmm::device_uvector<read_number>> vec_d_values_target;
            std::vector<rmm::device_uvector<int>> vec_d_offsets_target;

            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};

                CUDACHECK(cudaStreamWaitEvent(queryData->streams[g], queryData->callerEvent, 0));

                const int totalNumValuesTarget = queryData->pinnedData[g];

                auto* targetmr = rmm::mr::get_current_device_resource();
                rmm::device_uvector<read_number> d_values_target(totalNumValuesTarget, queryData->streams[g].getStream(), targetmr);
                rmm::device_uvector<int> d_offsets_target(numSequences + 1, queryData->streams[g].getStream(), targetmr);

                sgpuMinhashers[g]->retrieveValues(
                    *queryData->singlegpuMinhasherHandles[g].get(),
                    numSequences,
                    totalNumValuesTarget,
                    d_values_target.data(),
                    queryData->vec_d_numValuesPerSequence[g].data(),
                    d_offsets_target.data(), //numSequences + 1
                    queryData->streams[g],
                    targetmr
                );

                vec_d_values_target.push_back(std::move(d_values_target));
                vec_d_offsets_target.push_back(std::move(d_offsets_target));
            }

            int* h_gatherOffsets = queryData->pinnedData.data() + deviceIds.size();
            std::exclusive_scan(
                queryData->pinnedData.data(), 
                queryData->pinnedData.data() + deviceIds.size(), 
                h_gatherOffsets, 
                0
            );

            //gather results from targets
            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};
                const int totalNumValuesTarget = queryData->pinnedData[g];

                CUDACHECK(cudaMemcpyPeerAsync(
                    d_allValues.data() + h_gatherOffsets[g],
                    oldDeviceId,
                    vec_d_values_target[g].data(),
                    deviceIds[g],
                    sizeof(read_number) * totalNumValuesTarget,
                    queryData->streams[g]
                ));
                vec_d_values_target[g].release();
            }

            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};

                CUDACHECK(cudaMemcpyPeerAsync(
                    d_numValuesPerSequencePerGpu.data() + numSequences * g,
                    oldDeviceId,
                    queryData->vec_d_numValuesPerSequence[g].data(),
                    deviceIds[g],
                    sizeof(int) * numSequences,
                    queryData->streams[g]
                ));
                queryData->vec_d_numValuesPerSequence[g].release();
            }

            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};

                CUDACHECK(cudaMemcpyPeerAsync(
                    d_offsetsPerSequencePerGpu.data() + (numSequences+1) * g,
                    oldDeviceId,
                    vec_d_offsets_target[g].data(),
                    deviceIds[g],
                    sizeof(int) * (numSequences+1),
                    queryData->streams[g]
                ));
                vec_d_offsets_target[g].release();
            }

            rmm::device_uvector<int> d_gatherOffsets(deviceIds.size(), stream, mr);
            CUDACHECK(cudaMemcpyAsync(
                d_gatherOffsets.data(), 
                h_gatherOffsets, 
                sizeof(int) * deviceIds.size(), 
                H2D, 
                stream
            ));

            //join per-gpu streams to caller stream to wait for gathered results
            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};
                CUDACHECK(cudaEventRecord(queryData->events[g], queryData->streams[g]));
            }
            for(int g = 0; g < int(deviceIds.size()); g++){
                CUDACHECK(cudaStreamWaitEvent(stream, queryData->events[g]));
            }

            //copy values to output array, interleave results for same sequence
            multigpuminhasherkernels::copyToInterleavedKernel<<<numSequences, 128, 0, stream>>>(
                d_allValues.data(),
                d_offsetsPerSequencePerGpu.data(),
                d_numValuesPerSequencePerGpu.data(),
                d_gatherOffsets.data(),
                d_offsets,
                d_values,
                numSequences,
                deviceIds.size()
            ); CUDACHECKASYNC


            queryData->vec_d_numValuesPerSequence.clear();
        }


        void compact(cudaStream_t stream) override {
            CudaEvent event;
            event.record(stream);

            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};
                CUDACHECK(cudaStreamWaitEvent(cudaStreamPerThread, event, 0));

                sgpuMinhashers[g]->compact(cudaStreamPerThread);

                CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));
            }         
        }

        void constructionIsFinished(cudaStream_t stream) override{
            CudaEvent event;
            event.record(stream);

            std::vector<int> deviceIdsTmp;
            std::vector<std::unique_ptr<SingleGpuMinhasher>> sgpuMinhashersTmp;
            std::vector<std::vector<int>> hashFunctionIdsPerGpuTmp;

            for(int g = 0; g < int(deviceIds.size()); g++){
                cub::SwitchDevice sd{deviceIds[g]};
                CUDACHECK(cudaStreamWaitEvent(cudaStreamPerThread, event, 0));

                sgpuMinhashers[g]->constructionIsFinished(cudaStreamPerThread);

                CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));

                //only keep single-gpu minhashers which are used
                if(sgpuMinhashers[g]->getNumberOfMaps() > 0){
                    deviceIdsTmp.push_back(deviceIds[g]);
                    sgpuMinhashersTmp.push_back(std::move(sgpuMinhashers[g]));
                    hashFunctionIdsPerGpuTmp.push_back(std::move(hashFunctionIdsPerGpu[g]));
                }
            }
            
            std::swap(deviceIds, deviceIdsTmp);
            std::swap(sgpuMinhashers, sgpuMinhashersTmp);
            std::swap(hashFunctionIdsPerGpu, hashFunctionIdsPerGpuTmp);

            // std::cerr << "hashTableLocations:\n";
            // for(int i = 0; i < getNumberOfMaps(); i++){
            //     std::cerr << hashTableLocations[i] << " ";
            // }
            // std::cerr << "\n";

            // for(int g = 0; g < int(deviceIds.size()); g++){
            //     std::cerr << "hashFunctionIdsPerGpu " << g << " (id " << deviceIds[g] << ")" << "\n";
            //     for(auto x : hashFunctionIdsPerGpu[g]){
            //         std::cerr << x << " ";
            //     }
            //     std::cerr << "\n";
            // }

            // for(int g = 0; g < int(deviceIds.size()); g++){
            //     std::cerr << "actual stored hashFunctionIdsPerGpu " << g << " (id " << deviceIds[g] << ")" << "\n";
            //     const int num = sgpuMinhashers[g]->h_currentHashFunctionNumbers.size();
            //     for(int i = 0; i < num; i++){
            //         const int x = sgpuMinhashers[g]->h_currentHashFunctionNumbers[i];
            //         std::cerr << x << " ";
            //     }
            //     std::cerr << "\n";
            // }
        }

        MemoryUsage getMemoryInfo() const noexcept override{
            MemoryUsage mem{};

            for(const auto& minhasher : sgpuMinhashers){
                mem += minhasher->getMemoryInfo();
            }

            return mem;
        }

        MemoryUsage getMemoryInfo(const MinhasherHandle& handle) const noexcept override{
            return tempdataVector[handle.getId()]->getMemoryInfo();
        }

        int getNumResultsPerMapThreshold() const noexcept override{
            return resultsPerMapThreshold;
        }
        
        int getNumberOfMaps() const noexcept override{
            return hashTableLocations.size();
        }

        int getKmerSize() const noexcept override{
            return kmerSize;
        }

        void destroy(){
            for(auto& minhasher : sgpuMinhashers){
                DeviceSwitcher sd(minhasher->getDeviceId());
                minhasher->destroy();
            }
        }

        bool hasGpuTables() const noexcept override {
            return true;
        }

        void setThreadPool(ThreadPool* /*tp*/) override {}

        void setHostMemoryLimitForConstruction(std::size_t /*bytes*/) override{

        }

        void setDeviceMemoryLimitsForConstruction(const std::vector<std::size_t>&) override {

        }

        void writeToStream(std::ostream& /*os*/) const override{
            std::cerr << "MultiGpuMinhasher::writeToStream not supported\n";
        }

        int loadFromStream(std::ifstream& /*is*/, int /*numMapsUpperLimit*/) override{
            std::cerr << "MultiGpuMinhasher::loadFromStream not supported\n";
            return 0;
        } 

        bool canWriteToStream() const noexcept override { return false; };
        bool canLoadFromStream() const noexcept override { return false; };

private:        

        std::uint64_t getKmerMask() const{
            constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;

            return std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - getKmerSize()) * 2);
        }

        constexpr float getLoad() const noexcept{
            return 0.8f;
        }

        QueryData* getQueryDataFromHandle(const MinhasherHandle& queryHandle) const{
            std::shared_lock<SharedMutex> lock(sharedmutex);

            return tempdataVector[queryHandle.getId()].get();
        }

        int getNumberOfOccupiedDevices() const{
            int n = 0;
            for(const auto& vec : hashFunctionIdsPerGpu){
                if(vec.size() > 0){
                    n++;
                }
            }
            return n;
        }
        

        mutable int counter = 0;
        mutable SharedMutex sharedmutex{};

        Layout layout{};
        int maxNumKeys{};
        int kmerSize{};
        int resultsPerMapThreshold{};
        float loadfactor{};
        std::vector<int> deviceIds;
        std::vector<std::unique_ptr<SingleGpuMinhasher>> sgpuMinhashers{};
        mutable std::vector<std::unique_ptr<QueryData>> tempdataVector{};

        std::vector<std::vector<int>> hashFunctionIdsPerGpu{};
        std::vector<int> hashTableLocations{};
    };


}
}




#endif

#endif //#ifdef CARE_HAS_WARPCORE

