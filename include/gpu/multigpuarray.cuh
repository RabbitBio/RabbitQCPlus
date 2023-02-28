#ifndef MULTI_GPU_ARRAY_HPP
#define MULTI_GPU_ARRAY_HPP


#include <hpc_helpers.cuh>
#include <gpu/cudaerrorcheck.cuh>
#include <gpu/singlegpu2darray.cuh>

#include <memorymanagement.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>
#include <limits>
#include <stdexcept>

#include <cub/cub.cuh>
#include <cooperative_groups.h>

#include <thrust/iterator/transform_iterator.h>

#include <rmm/device_uvector.hpp>
#include <rmm/mr/device/per_device_resource.hpp>

namespace cg = cooperative_groups;


namespace MultiGpu2dArrayKernels{

    template<class T, class IndexIterator>
    __global__
    void gatherKernel(
        const T* __restrict__ src, 
        size_t srcRowPitchInBytes, 
        size_t numColumns,
        size_t numRows,
        T* __restrict__ dest, 
        size_t destRowPitchInBytes, 
        IndexIterator indices, 
        size_t numIndices
    ){
        if(numIndices == 0) return;

        auto group = cg::this_grid();

        TwoDimensionalArray<T> array(
            numRows,
            numColumns,
            srcRowPitchInBytes,
            const_cast<T*>(src)
        );

        array.gather(group, dest, destRowPitchInBytes, indices, numIndices);
    }

    template<class T, class IndexIterator>
    __global__
    void gatherKernel(
        const T* __restrict__ src, 
        size_t srcRowPitchInBytes, 
        size_t numRows,
        size_t numColumns,
        T* __restrict__ dest, 
        size_t destRowPitchInBytes, 
        IndexIterator indices, 
        const size_t* __restrict__ numIndicesPtr
    ){
        const size_t numIndices = *numIndicesPtr;

        if(numIndices == 0) return;

        auto group = cg::this_grid();

        TwoDimensionalArray<T> array(
            numRows,
            numColumns,
            srcRowPitchInBytes,
            const_cast<T*>(src)
        );

        array.gather(group, dest, destRowPitchInBytes, indices, numIndices);
    }

    template<class T, class IndexIterator>
    __global__
    void scatterKernel(
        const T* __restrict__ src, 
        size_t srcRowPitchInBytes, 
        size_t numRows,
        size_t numColumns,
        T* __restrict__ dest, 
        size_t destRowPitchInBytes, 
        IndexIterator indices, 
        size_t numIndices
    ){
        if(numIndices == 0) return;

        auto group = cg::this_grid();

        TwoDimensionalArray<T> array(
            numRows,
            numColumns,
            destRowPitchInBytes,
            dest
        );

        array.scatter(group, src, srcRowPitchInBytes, indices, numIndices);
    }

    template<class T, class IndexIterator>
    __global__
    void scatterKernel(
        const T* __restrict__ src, 
        size_t srcRowPitchInBytes, 
        size_t numRows,
        size_t numColumns,
        T* __restrict__ dest, 
        size_t destRowPitchInBytes, 
        IndexIterator indices, 
        const size_t* __restrict__ numIndicesPtr
    ){
        const size_t numIndices = *numIndicesPtr;
        if(numIndices == 0) return;

        auto group = cg::this_grid();

        TwoDimensionalArray<T> array(
            numRows,
            numColumns,
            destRowPitchInBytes,
            dest
        );

        array.scatter(group, src, srcRowPitchInBytes, indices, numIndices);
    }

    template<class T>
    __global__
    void copyKernel(const T* __restrict__ src, const size_t* __restrict__ numElementsPtr, T* __restrict__ dst){
        const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
        const size_t stride = gridDim.x * blockDim.x;

        const size_t numElements = *numElementsPtr;

        for(size_t i = tid; i < numElements; i += stride){
            dst[i] = src[i];
        }
    }

    template<class T>
    __global__
    void copyKernel(const T* __restrict__ src, size_t numElements, T* __restrict__ dst){
        const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
        const size_t stride = gridDim.x * blockDim.x;

        for(size_t i = tid; i < numElements; i += stride){
            dst[i] = src[i];
        }
    }



    template<class IndexType, int maxNumGpus = 8>
    __global__
    void partitionSplitKernel(
        IndexType* __restrict__ splitIndices, // numIndices elements per partition
        size_t* __restrict__ splitDestinationPositions, // numIndices elements per partition
        size_t* __restrict__ numSplitIndicesPerPartition, 
        int numPartitions,
        const size_t* __restrict__ partitionOffsetsPS,
        size_t numIndices,
        const IndexType* __restrict__ indices
    ){

        assert(numPartitions <= maxNumGpus+1);

        using Index_t = size_t;

        auto atomicAggInc = [](auto& group, Index_t* counter){
            Index_t group_res;
            if(group.thread_rank() == 0){
                group_res = atomicAdd((unsigned long long*)counter, (unsigned long long)(group.size()));
            }
            return group.shfl(group_res, 0) + group.thread_rank();
        };

        //save prefixsum in registers
        Index_t reg_partitionOffsetsPS[maxNumGpus+1];

        #pragma unroll
        for(int i = 0; i < maxNumGpus+1; i++){
            if(i < numPartitions+1){
                reg_partitionOffsetsPS[i] = partitionOffsetsPS[i];
            }else{
                reg_partitionOffsetsPS[i] = partitionOffsetsPS[numPartitions];
            }
        }

        const Index_t numIds = numIndices;

        for(Index_t tid = threadIdx.x + blockIdx.x * Index_t(blockDim.x); 
                tid < numIds; 
                tid += Index_t(blockDim.x) * gridDim.x){
            
            const Index_t elementIndex = indices[tid];
            int location = -1;

            #pragma unroll
            for(int i = 0; i < maxNumGpus+1; i++){
                if(i < numPartitions 
                        && reg_partitionOffsetsPS[i] <= elementIndex 
                        && elementIndex < reg_partitionOffsetsPS[i+1]){
                    location = i;
                    break;
                }
            }

#if __CUDACC_VER_MAJOR__ >= 11 && __CUDA_ARCH__ >= 700
            if(location != -1){
                auto g = cg::coalesced_threads();
                //partition into groups of threads with the same location. Cannot use a tiled_partition<N> as input, because for-loop may cause a deadlock:
                //"The implementation may cause the calling thread to wait until all the members of the parent group have invoked the operation before resuming execution."
                auto partitionGroup = cg::labeled_partition(g, location);

                const Index_t j = atomicAggInc(partitionGroup, &numSplitIndicesPerPartition[location]);
                splitIndices[location * numIds + j] = elementIndex;
                splitDestinationPositions[location * numIds + j] = tid;
            }
              
#else
            for(int i = 0; i < numPartitions; ++i){
                if(i == location){
                    auto group = cg::coalesced_threads();
                    const Index_t j = atomicAggInc(group, &numSplitIndicesPerPartition[i]);
                    splitIndices[i * numIds + j] = elementIndex;
                    splitDestinationPositions[i * numIds + j] = tid;
                    break;
                }
            }
#endif 

        }

    }

    template<class T>
    struct LinearAccessFunctor{
        HOSTDEVICEQUALIFIER
        LinearAccessFunctor(T* ptr) : data(ptr){}

        HOSTDEVICEQUALIFIER
        T& operator()(size_t pos){
            return data[pos];
        }

        T* data;
    };

    template<class T>
    struct SubtractFunctor{
        T value;
        __host__ __device__
        SubtractFunctor(T val) : value(val){}

        __host__ __device__
        T operator()(const T& x) const noexcept{
            return x - value;
        }
    };
}




enum class MultiGpu2dArrayLayout {FirstFit, EvenShare};
enum class MultiGpu2dArrayInitMode {MustFitCompletely, CanDiscardRows};

template<class T, class IndexType = int>
class MultiGpu2dArray{
public:

    template<class W>
    using HostBuffer = helpers::SimpleAllocationPinnedHost<W>;
    template<class W>
    using DeviceBuffer = helpers::SimpleAllocationDevice<W>;
    //using DeviceBuffer = helpers::SimpleAllocationPinnedHost<W>;



    MultiGpu2dArray()
    : alignmentInBytes(sizeof(T)){

    }

    MultiGpu2dArray(
        size_t numRows, 
        size_t numColumns, 
        size_t alignmentInBytes,
        std::vector<int> dDeviceIds, // data resides on these gpus
        std::vector<size_t> memoryLimits,
        MultiGpu2dArrayLayout layout,
        MultiGpu2dArrayInitMode initMode = MultiGpu2dArrayInitMode::MustFitCompletely) 
    : alignmentInBytes(alignmentInBytes), 
        dataDeviceIds(dDeviceIds)
    {

        assert(alignmentInBytes > 0);
        //assert(alignmentInBytes % sizeof(T) == 0);

        init(numRows, numColumns, memoryLimits, layout, initMode);
    }

    ~MultiGpu2dArray(){
        destroy();
    }

    MultiGpu2dArray(const MultiGpu2dArray& rhs)
        : MultiGpu2dArray(rhs.numRows, rhs.numColumns, rhs.alignmentInBytes)
    {
        gpuArrays.clear();
        for(const auto& otherarray : rhs.gpuArrays){
            gpuArrays.emplace_back(otherarray);
        }
    }

    MultiGpu2dArray(MultiGpu2dArray&& other) noexcept
        : MultiGpu2dArray()
    {
        swap(*this, other);
    }

    MultiGpu2dArray& operator=(MultiGpu2dArray other){
        swap(*this, other);

        return *this;
    }

    friend void swap(MultiGpu2dArray& l, MultiGpu2dArray& r) noexcept
    {
        std::swap(l.numRows, r.numRows);
        std::swap(l.numColumns, r.numColumns);
        std::swap(l.alignmentInBytes, r.alignmentInBytes);
        std::swap(l.gpuArrays, r.gpuArrays);
        std::swap(l.dataDeviceIds, r.dataDeviceIds);
        std::swap(l.usedDeviceIds, r.usedDeviceIds);
        std::swap(l.h_numRowsPerGpu, r.h_numRowsPerGpu);
        std::swap(l.h_numRowsPerGpuPrefixSum, r.h_numRowsPerGpuPrefixSum);
    }

    void destroy(){
        numRows = 0;
        numColumns = 0;
        alignmentInBytes = 0;
        gpuArrays.clear();
        dataDeviceIds.clear();
        usedDeviceIds.clear();
        h_numRowsPerGpu.destroy();
        h_numRowsPerGpuPrefixSum.destroy();
    }

    struct HandleStruct{
        struct PerDevice{
            PerDevice()
            {        
                cudaGetDevice(&deviceId);
            }

            int deviceId{};

            CudaStream stream{};
            CudaEvent event{cudaEventDisableTiming};
        };

        struct PerCaller{           
            CudaEvent event{cudaEventDisableTiming};
            std::vector<CudaEvent> events{};
            std::vector<CudaStream> streams{};
        };

        HandleStruct(){
            CUDACHECK(cudaGetDevice(&deviceId));            
        }

        HandleStruct(HandleStruct&&) = default;
        HandleStruct(const HandleStruct&) = delete;

        int deviceId = 0;
        std::vector<PerDevice> deviceBuffers{};
        PerCaller callerBuffers{};
        helpers::SimpleAllocationPinnedHost<std::size_t> multisplitTemp{};

        MemoryUsage getMemoryInfo() const{
            MemoryUsage result{};

            return result;
        }
    };

    using Handle = std::shared_ptr<HandleStruct>;

    MemoryUsage getMemoryInfo(const Handle& handle) const{
        MemoryUsage result{};

        return handle->getMemoryInfo();
    }

    Handle makeHandle() const{
        Handle handle = std::make_shared<HandleStruct>();

        const int numDataGpus = usedDeviceIds.size();
        for(int d = 0; d < numDataGpus; d++){
            const int deviceId = usedDeviceIds[d];
            cub::SwitchDevice sd(deviceId);

            typename HandleStruct::PerDevice deviceBuffers; 

            handle->deviceBuffers.emplace_back(std::move(deviceBuffers));
        }

        if(h_numRowsPerGpu.size() > 0){

            auto& callerBuffers = handle->callerBuffers;
            callerBuffers.streams.resize(h_numRowsPerGpu.size());
            
            for(int i = 0; i < int(h_numRowsPerGpu.size()); i++){
                callerBuffers.events.emplace_back(cudaEventDisableTiming);
            }
            CUDACHECKASYNC;
        } 

        handle->multisplitTemp.resize(std::max(usedDeviceIds.size(), dataDeviceIds.size()));

        return handle;
    }

    bool tryReplication(std::vector<std::size_t> memoryLimits){
        if(!isReplicatedSingleGpu && usedDeviceIds.size() == 1 && usedDeviceIds.size() < dataDeviceIds.size()){
            assert(gpuArrays.size() == 1);

            auto memoryUsage = gpuArrays[0]->getMemoryInfo();
            const std::size_t requiredMemory = memoryUsage.device[gpuArrays[0]->getDeviceId()];

            std::vector<std::unique_ptr<Gpu2dArrayManaged<T>>> replicas;

            for(std::size_t i = usedDeviceIds.size(); i < dataDeviceIds.size(); i++){
                if(memoryLimits[i] < requiredMemory){
                    return false;
                }
                replicas.emplace_back(gpuArrays[0]->makeCopy(dataDeviceIds[i]));
            }

            bool ok = std::all_of(replicas.begin(), replicas.end(), [](const auto& uniqueptr){ return bool(uniqueptr); });

            if(ok){
                gpuArrays.insert(gpuArrays.end(), std::make_move_iterator(replicas.begin()), std::make_move_iterator(replicas.end()));
                usedDeviceIds.insert(usedDeviceIds.end(), dataDeviceIds.begin() + usedDeviceIds.size(), dataDeviceIds.end());
            }

            isReplicatedSingleGpu = ok;

            return ok;

        }else{
            return false;
        }
    }

    void gather(
        Handle& handle, 
        T* d_dest, 
        size_t destRowPitchInBytes, 
        const IndexType* d_indices, 
        size_t numIndices, 
        cudaStream_t destStream = 0
    ) const{

        if(getNumRows() == 0) return;

        gather_internal(
            handle,
            d_dest,
            destRowPitchInBytes,
            d_indices,
            numIndices,
            false,
            destStream
        );
    }

    void gather(
        Handle& handle, 
        T* d_dest, 
        size_t destRowPitchInBytes, 
        const IndexType* d_indices, 
        size_t numIndices, 
        bool mayContainInvalidIndices,
        cudaStream_t destStream = 0
    ) const{

        if(getNumRows() == 0) return;

        gather_internal(
            handle,
            d_dest,
            destRowPitchInBytes,
            d_indices,
            numIndices,
            mayContainInvalidIndices,
            destStream
        );
    }

    void scatter(
        Handle& handle, 
        const T* d_src, 
        size_t srcRowPitchInBytes, 
        const IndexType* d_indices, 
        size_t numIndices, 
        cudaStream_t srcStream = 0
    ) const{

        if(getNumRows() == 0) return;

        scatter_internal(
            handle,
            d_src,
            srcRowPitchInBytes,
            d_indices,
            numIndices,
            false,
            srcStream
        );
    }

    void scatter(
        Handle& handle, 
        const T* d_src, 
        size_t srcRowPitchInBytes, 
        const IndexType* d_indices, 
        size_t numIndices, 
        bool mayContainInvalidIndices,
        cudaStream_t srcStream = 0
    ) const{

        if(getNumRows() == 0) return;

        scatter_internal(
            handle,
            d_src,
            srcRowPitchInBytes,
            d_indices,
            numIndices,
            mayContainInvalidIndices,
            srcStream
        );
    }


    void gatherContiguous(Handle& handle, T* d_dest, size_t destRowPitchInBytes, size_t rowBegin, size_t numRowsToGather, cudaStream_t stream = 0) const{
        if(numRowsToGather == 0) return;

        std::size_t rowEnd = rowBegin + numRowsToGather;

        std::vector<std::size_t> numGatherPerDevice(usedDeviceIds.size(), 0);
        std::vector<std::size_t> firstRowPerDevice(usedDeviceIds.size(), 0);

        for(std::size_t i = 0, cur = rowBegin; i < usedDeviceIds.size(); i++){
            if(cur < h_numRowsPerGpuPrefixSum[i+1]){
                auto myEnd = std::min(h_numRowsPerGpuPrefixSum[i+1], rowEnd);
                const auto myNum = myEnd - cur;
                firstRowPerDevice[i] = cur - h_numRowsPerGpuPrefixSum[i];
                numGatherPerDevice[i] = myNum;
                cur += myNum;
            }
        }

        int oldDeviceId = 0;
        CUDACHECK(cudaGetDevice(&oldDeviceId));

        auto& callerBuffers = handle->callerBuffers;
        CUDACHECK(cudaEventRecord(callerBuffers.event, stream));

        #if 1

            for(std::size_t i = 0, offset = 0; i < usedDeviceIds.size(); i++){
                const std::size_t num = numGatherPerDevice[i];
                if(num > 0){
                    const auto& gpuArray = gpuArrays[i];

                    CUDACHECK(cudaStreamWaitEvent(callerBuffers.streams[i], callerBuffers.event));

                    gpuArray->gatherContiguousPeer(
                        oldDeviceId,
                        (T*)(((char*)d_dest) + destRowPitchInBytes * offset),
                        destRowPitchInBytes,
                        firstRowPerDevice[i],
                        num,
                        callerBuffers.streams[i]
                    );

                    CUDACHECK(cudaEventRecord(callerBuffers.events[i], callerBuffers.streams[i]));

                    offset += num;
                }
            }

            for(std::size_t i = 0; i < usedDeviceIds.size(); i++){
                const std::size_t num = numGatherPerDevice[i];
                if(num > 0){                
                    CUDACHECK(cudaStreamWaitEvent(stream, callerBuffers.events[i]));
                }
            }
        #else

        std::vector<rmm::device_uvector<T>> vec_d_gathered;

        for(std::size_t i = 0; i < usedDeviceIds.size(); i++){
            const std::size_t num = numGatherPerDevice[i];
            if(num > 0){
                const auto& gpuArray = gpuArrays[i];
                auto& deviceBuffers = handle->deviceBuffers[i];
                const int deviceId = gpuArray->getDeviceId();

                cub::SwitchDevice sd{deviceId};

                rmm::device_uvector<T> d_gathered(destRowPitchInBytes * num, deviceBuffers.stream.getStream());
                gpuArray->gatherContiguous(
                    d_gathered.data(),
                    destRowPitchInBytes,
                    firstRowPerDevice[i],
                    num,
                    deviceBuffers.stream
                );
                vec_d_gathered.push_back(std::move(d_gathered));
            }else{
                vec_d_gathered.emplace_back(0, stream);
            }
        }

        for(std::size_t i = 0, offset = 0; i < usedDeviceIds.size(); i++){
            const std::size_t num = numGatherPerDevice[i];
            if(num > 0){
                const auto& gpuArray = gpuArrays[i];
                auto& deviceBuffers = handle->deviceBuffers[i];
                const int deviceId = gpuArray->getDeviceId();

                cub::SwitchDevice sd{deviceId};
                CUDACHECK(cudaStreamWaitEvent(deviceBuffers.stream, callerBuffers.event));

                copy(
                    ((char*)d_dest) + destRowPitchInBytes * offset, 
                    oldDeviceId, 
                    vec_d_gathered[i].data(), 
                    deviceId, 
                    destRowPitchInBytes * num, 
                    deviceBuffers.stream
                );

                offset += num;
                vec_d_gathered[i].release();
                CUDACHECK(cudaEventRecord(deviceBuffers.event, deviceBuffers.stream));
                // for(auto id : {0,1}){
                //     cub::SwitchDevice sd{id};
                //     CUDACHECK(cudaDeviceSynchronize());
                // }
            }
        }

        for(std::size_t i = 0; i < usedDeviceIds.size(); i++){
            const std::size_t num = numGatherPerDevice[i];
            if(num > 0){
                
                CUDACHECK(cudaStreamWaitEvent(stream, handle->deviceBuffers[i].event));
            }
        }
        #endif
    }

    void scatterContiguous(Handle& handle, T* d_src, size_t srcRowPitchInBytes, size_t rowBegin, size_t numRowsToScatter, cudaStream_t stream = 0) const{
        if(numRowsToScatter == 0) return;

        std::size_t rowEnd = rowBegin + numRowsToScatter;

        std::vector<std::size_t> numScatterPerDevice(usedDeviceIds.size(), 0);
        std::vector<std::size_t> firstRowPerDevice(usedDeviceIds.size(), 0);

        for(std::size_t i = 0, cur = rowBegin; i < usedDeviceIds.size(); i++){
            if(cur < h_numRowsPerGpuPrefixSum[i+1]){
                auto myEnd = std::min(h_numRowsPerGpuPrefixSum[i+1], rowEnd);
                const auto myNum = myEnd - cur;
                firstRowPerDevice[i] = cur - h_numRowsPerGpuPrefixSum[i];
                numScatterPerDevice[i] = myNum;
                cur += myNum;
            }
        }

        int oldDeviceId = 0;
        CUDACHECK(cudaGetDevice(&oldDeviceId));

        auto& callerBuffers = handle->callerBuffers;
        CUDACHECK(cudaEventRecord(callerBuffers.event, stream));

        #if 1
            for(std::size_t i = 0, offset = 0; i < usedDeviceIds.size(); i++){
                const std::size_t num = numScatterPerDevice[i];
                if(num > 0){
                    const auto& gpuArray = gpuArrays[i];
                    auto& deviceBuffers = handle->deviceBuffers[i];
                    const int deviceId = gpuArray->getDeviceId();

                    CUDACHECK(cudaStreamWaitEvent(callerBuffers.streams[i], callerBuffers.event));

                    gpuArray->scatterContiguousPeer(
                        oldDeviceId,
                        (T*)(((char*)d_src) + srcRowPitchInBytes * offset),
                        srcRowPitchInBytes,
                        firstRowPerDevice[i],
                        num,
                        callerBuffers.streams[i]
                    );

                    CUDACHECK(cudaEventRecord(callerBuffers.events[i], callerBuffers.streams[i]));

                    offset += num;
                }
            }

            for(std::size_t i = 0; i < usedDeviceIds.size(); i++){
                const std::size_t num = numScatterPerDevice[i];
                if(num > 0){                
                    CUDACHECK(cudaStreamWaitEvent(stream, callerBuffers.events[i]));
                }
            }
        #else

        std::vector<rmm::device_uvector<T>> vec_d_scatter;

        //copy chunks to target devices
        for(std::size_t i = 0, offset = 0; i < usedDeviceIds.size(); i++){
            const std::size_t num = numScatterPerDevice[i];
            if(num > 0){
                const auto& gpuArray = gpuArrays[i];
                auto& deviceBuffers = handle->deviceBuffers[i];
                const int deviceId = gpuArray->getDeviceId();

                cub::SwitchDevice sd{deviceId};
                CUDACHECK(cudaStreamWaitEvent(deviceBuffers.stream, callerBuffers.event));
                rmm::device_uvector<T> d_scatter(srcRowPitchInBytes * num, deviceBuffers.stream.getStream());

                copy(
                    d_scatter.data(), 
                    deviceId,
                    ((char*)d_src) + srcRowPitchInBytes * offset,
                    oldDeviceId,
                    srcRowPitchInBytes * num, 
                    deviceBuffers.stream
                );

                offset += num;

                vec_d_scatter.push_back(std::move(d_scatter));
            }else{
                vec_d_scatter.emplace_back(0, stream);
            }
        }

        //scatter on target devices

        for(std::size_t i = 0; i < usedDeviceIds.size(); i++){
            const std::size_t num = numScatterPerDevice[i];
            if(num > 0){
                const auto& gpuArray = gpuArrays[i];
                auto& deviceBuffers = handle->deviceBuffers[i];
                const int deviceId = gpuArray->getDeviceId();

                cub::SwitchDevice sd{deviceId};

                gpuArray->scatterContiguous(
                    vec_d_scatter[i].data(),
                    srcRowPitchInBytes,
                    firstRowPerDevice[i],
                    num,
                    deviceBuffers.stream
                );
                vec_d_scatter[i].release();
                CUDACHECK(cudaEventRecord(deviceBuffers.event, deviceBuffers.stream));
            }
        }

        for(std::size_t i = 0; i < usedDeviceIds.size(); i++){
            const std::size_t num = numScatterPerDevice[i];
            if(num > 0){                
                CUDACHECK(cudaStreamWaitEvent(stream, handle->deviceBuffers[i].event));
            }
        }

        #endif

    }

    size_t getNumRows() const noexcept{
        return numRows;
    }

    size_t getNumColumns() const noexcept{
        return numColumns;
    }

    size_t getAlignmentInBytes() const noexcept{
        return alignmentInBytes;
    }

    bool isSingleGpu() const noexcept{
        return gpuArrays.size() == 1;
    }

    std::vector<size_t> getRowDistribution() const{
        std::vector<size_t> result(h_numRowsPerGpu.size());
        std::copy_n(h_numRowsPerGpu.get(), h_numRowsPerGpu.size(), result.begin());

        return result;
    }

    MemoryUsage getMemoryInfo() const{
        MemoryUsage result{};

        result.host = 0;
        for(const auto& ptr : gpuArrays){
            result += ptr->getMemoryInfo();
        }

        return result;
    }

private:

    void gather_internal(
        Handle& handle, T* d_dest, 
        size_t destRowPitchInBytes, 
        const IndexType* d_indices, 
        size_t numIndices, 
        bool mayContainInvalidIndices,
        cudaStream_t destStream
    ) const{
        if(numIndices == 0) return;

        int destDeviceId = 0;
        CUDACHECK(cudaGetDevice(&destDeviceId));

        if((isSingleGpu() || isReplicatedSingleGpu) && !mayContainInvalidIndices){
            auto it = std::find(usedDeviceIds.begin(), usedDeviceIds.end(), destDeviceId);
            if(it != usedDeviceIds.end()){
                //all data to be gathered resides on the destination device.
                const int position = std::distance(usedDeviceIds.begin(), it);
                
                gpuArrays[position]->gather(d_dest, destRowPitchInBytes, d_indices, numIndices, destStream);
            }else{
                //array was not constructed on current device. use peer access from any gpu

                const auto& gpuArray = gpuArrays[0];
                auto& deviceBuffers = handle->deviceBuffers[0];
                const int deviceId = gpuArray->getDeviceId();

                CUDACHECK(cudaEventRecord(handle->callerBuffers.event, destStream));
                
                CUDACHECK(cudaSetDevice(deviceId));
                CUDACHECK(cudaStreamWaitEvent(deviceBuffers.stream, handle->callerBuffers.event, 0));

                //copy indices to other gpu
                rmm::device_uvector<IndexType> d_indices_target(numIndices, deviceBuffers.stream.getStream());
                copy(
                    d_indices_target.data(), 
                    deviceId, 
                    d_indices, 
                    destDeviceId, 
                    sizeof(IndexType) * numIndices, 
                    deviceBuffers.stream
                );

                rmm::device_buffer d_temp(destRowPitchInBytes * numIndices, deviceBuffers.stream.getStream());
                gpuArrays[0]->gather(reinterpret_cast<T*>(d_temp.data()), destRowPitchInBytes, d_indices_target.data(), numIndices, deviceBuffers.stream);

                //copy gathered data to us
                copy(
                    d_dest, 
                    destDeviceId, 
                    d_temp.data(), 
                    deviceId, 
                    destRowPitchInBytes * numIndices, 
                    deviceBuffers.stream
                );

                CUDACHECK(cudaEventRecord(deviceBuffers.event, deviceBuffers.stream.getStream()));

                //wait on destStream until gathered data is ready
                CUDACHECK(cudaSetDevice(destDeviceId));
                CUDACHECK(cudaStreamWaitEvent(destStream, deviceBuffers.event, 0));
            }
        }else{

            //perform multisplit to distribute indices and filter out invalid indices. then perform local gathers,
            // and scatter into result array

            const int numGpus = gpuArrays.size();
            const int numDistinctGpus = h_numRowsPerGpu.size();

            if(!isReplicatedSingleGpu){
                assert(numGpus == numDistinctGpus);
            }else{
                assert(1 == numDistinctGpus);
            }

            auto& callerBuffers = handle->callerBuffers;

            MultiSplitResult splits = multiSplit(d_indices, numIndices, handle->multisplitTemp.data(), destStream);

            std::vector<std::size_t> numSelectedPrefixSum(numDistinctGpus,0);
            std::partial_sum(
                splits.h_numSelectedPerGpu.data(),
                splits.h_numSelectedPerGpu.data() + numDistinctGpus - 1,
                numSelectedPrefixSum.begin() + 1
            );

            rmm::device_uvector<char> callerBuffers_d_dataCommunicationBuffer(numIndices * destRowPitchInBytes, destStream);

            CUDACHECK(cudaEventRecord(handle->callerBuffers.event, destStream));

            for(int d = 0; d < numDistinctGpus; d++){
                const int num = splits.h_numSelectedPerGpu[d];
                if(num > 0){
                    CUDACHECK(cudaStreamWaitEvent(callerBuffers.streams[d], handle->callerBuffers.event));
                }
            }


            {
                std::vector<IndexType*> all_deviceBuffers_d_selectedIndices(numDistinctGpus, nullptr);
                std::vector<char*> all_deviceBuffers_d_dataCommunicationBuffer(numDistinctGpus, nullptr);
                std::vector<rmm::mr::device_memory_resource*> all_mrs(numDistinctGpus, nullptr);

                //allocate remote buffers
                for(int d = 0; d < numDistinctGpus; d++){
                    const int num = splits.h_numSelectedPerGpu[d];

                    if(num > 0){
                        const auto& gpuArray = gpuArrays[d];
                        auto& deviceBuffers = handle->deviceBuffers[d];
                        const int deviceId = gpuArray->getDeviceId();

                        cub::SwitchDevice sd(deviceId);

                        all_mrs[d] = rmm::mr::get_current_device_resource();
                        all_deviceBuffers_d_selectedIndices[d] = reinterpret_cast<IndexType*>(all_mrs[d]->allocate(sizeof(IndexType) * num, deviceBuffers.stream.getStream()));
                        all_deviceBuffers_d_dataCommunicationBuffer[d] = reinterpret_cast<char*>(all_mrs[d]->allocate(destRowPitchInBytes * num, deviceBuffers.stream.getStream()));
                    }
                }

                //copy selected indices to remote buffer
                for(int d = 0; d < numDistinctGpus; d++){
                    const int num = splits.h_numSelectedPerGpu[d];

                    if(num > 0){
                        const auto& gpuArray = gpuArrays[d];
                        auto& deviceBuffers = handle->deviceBuffers[d];
                        const int deviceId = gpuArray->getDeviceId();

                        cub::SwitchDevice sd(deviceId);

                        //copy selected indices to remote device
                        copy(
                            all_deviceBuffers_d_selectedIndices[d], 
                            deviceId, 
                            splits.d_selectedIndices.data() + d * numIndices, 
                            destDeviceId, 
                            sizeof(IndexType) * num, 
                            deviceBuffers.stream
                        );
                    }
                }

                //gather on remote gpus
                for(int d = 0; d < numDistinctGpus; d++){
                    const int num = splits.h_numSelectedPerGpu[d];

                    if(num > 0){
                        const auto& gpuArray = gpuArrays[d];
                        auto& deviceBuffers = handle->deviceBuffers[d];
                        const int deviceId = gpuArray->getDeviceId();

                        cub::SwitchDevice sd(deviceId);                      

                        auto gatherIndexIterator = thrust::make_transform_iterator(
                            all_deviceBuffers_d_selectedIndices[d],
                            MultiGpu2dArrayKernels::SubtractFunctor<IndexType>(h_numRowsPerGpuPrefixSum[d])
                        );

                        gpuArray->gather(
                            (T*)all_deviceBuffers_d_dataCommunicationBuffer[d], 
                            destRowPitchInBytes, 
                            gatherIndexIterator, 
                            num, 
                            deviceBuffers.stream
                        );
                    }
                }

                //copy remote gathered data to caller
                for(int d = 0; d < numDistinctGpus; d++){
                    const int num = splits.h_numSelectedPerGpu[d];

                    if(num > 0){
                        const auto& gpuArray = gpuArrays[d];
                        auto& deviceBuffers = handle->deviceBuffers[d];
                        const int deviceId = gpuArray->getDeviceId();

                        CUDACHECK(cudaSetDevice(deviceId));

                        copy(
                            callerBuffers_d_dataCommunicationBuffer.data() + destRowPitchInBytes * numSelectedPrefixSum[d], 
                            destDeviceId, 
                            all_deviceBuffers_d_dataCommunicationBuffer[d], 
                            deviceId, 
                            destRowPitchInBytes * num, 
                            deviceBuffers.stream
                        );

                        //free temp allocations of this device after data has been copied
                        all_mrs[d]->deallocate(all_deviceBuffers_d_dataCommunicationBuffer[d], destRowPitchInBytes * num, deviceBuffers.stream.getStream());
                        all_mrs[d]->deallocate(all_deviceBuffers_d_selectedIndices[d], sizeof(IndexType) * num, deviceBuffers.stream.getStream());

                        CUDACHECK(cudaEventRecord(deviceBuffers.event, deviceBuffers.stream));

                        CUDACHECK(cudaSetDevice(destDeviceId));
                        CUDACHECK(cudaStreamWaitEvent(callerBuffers.streams[d], deviceBuffers.event, 0));
                    }
                }

                //scatter into result array
                for(int d = 0; d < numDistinctGpus; d++){
                    const int num = splits.h_numSelectedPerGpu[d];

                    if(num > 0){

                        dim3 block(128, 1, 1);
                        dim3 grid(SDIV(numIndices, block.x), 1, 1);

                        MultiGpu2dArrayKernels::scatterKernel<<<grid, block, 0, callerBuffers.streams[d]>>>(
                            (const T*)(callerBuffers_d_dataCommunicationBuffer.data() + destRowPitchInBytes * numSelectedPrefixSum[d]), 
                            destRowPitchInBytes, 
                            numIndices,
                            getNumColumns(),
                            d_dest, 
                            destRowPitchInBytes, 
                            splits.d_selectedPositions.data() + d * numIndices, 
                            num
                        ); CUDACHECKASYNC; 

                        CUDACHECK(cudaEventRecord(callerBuffers.events[d], callerBuffers.streams[d]));
                    }
                }

                //join streams
                for(int d = 0; d < numDistinctGpus; d++){
                    const int num = splits.h_numSelectedPerGpu[d];

                    if(num > 0){
                        CUDACHECK(cudaStreamWaitEvent(destStream, callerBuffers.events[d], 0));
                    }                
                }             
            }
        }
    }

    void scatter_internal(
        Handle& handle, 
        const T* d_src, 
        size_t srcRowPitchInBytes, 
        const IndexType* d_indices, 
        size_t numIndices, 
        bool mayContainInvalidIndices,
        cudaStream_t srcStream = 0
    ) const{
        assert(!isReplicatedSingleGpu); //not implemented

        if(numIndices == 0) return;

        int srcDeviceId = 0;
        CUDACHECK(cudaGetDevice(&srcDeviceId));

        if(isSingleGpu() && !mayContainInvalidIndices && gpuArrays[0]->getDeviceId() == srcDeviceId){ //if all data should be scattered to same device

            gpuArrays[0]->scatter(d_src, srcRowPitchInBytes, d_indices, numIndices, srcStream);
            return;
        }else{
            const int numGpus = gpuArrays.size();
            const int numDistinctGpus = h_numRowsPerGpu.size();

            if(!isReplicatedSingleGpu){
                assert(numGpus == numDistinctGpus);
            }else{
                assert(1 == numDistinctGpus);
            }

            auto& callerBuffers = handle->callerBuffers;
            auto* srcMR = rmm::mr::get_current_device_resource();

            rmm::device_uvector<char> callerBuffers_d_dataCommunicationBuffer(numIndices * srcRowPitchInBytes, srcStream, srcMR);

            MultiSplitResult splits = multiSplit(d_indices, numIndices, handle->multisplitTemp.data(), srcStream);

            std::vector<std::size_t> numSelectedPrefixSum(numDistinctGpus,0);
            std::partial_sum(
                splits.h_numSelectedPerGpu.data(),
                splits.h_numSelectedPerGpu.data() + numDistinctGpus - 1,
                numSelectedPrefixSum.begin() + 1
            );

            for(int d = 0; d < numDistinctGpus; d++){
                const int numSelected = splits.h_numSelectedPerGpu[d];
                if(numSelected > 0){
                    const auto& gpuArray = gpuArrays[d];
                    const int deviceId = gpuArray->getDeviceId();
                    auto& targetDeviceBuffers = handle->deviceBuffers[d];  
                    auto* mr = rmm::mr::get_per_device_resource(rmm::cuda_device_id{deviceId});

                    //local gather data which should be scattered to current gpuArray
                    dim3 block(128, 1, 1);
                    dim3 grid(SDIV(numSelected, block.x), 1, 1);

                    MultiGpu2dArrayKernels::gatherKernel<<<grid, block, 0, callerBuffers.streams[d]>>>(
                        d_src, 
                        srcRowPitchInBytes, 
                        numIndices,
                        getNumColumns(),
                        (T*)(callerBuffers_d_dataCommunicationBuffer.data() + numSelectedPrefixSum[d] * srcRowPitchInBytes), 
                        srcRowPitchInBytes, 
                        splits.d_selectedPositions.data() + numIndices * d, 
                        splits.d_numSelectedPerGpu.data() + d
                    ); CUDACHECKASYNC;

                    CUDACHECK(cudaEventRecord(callerBuffers.events[d], callerBuffers.streams[d]));

                    cub::SwitchDevice sd{deviceId};
                    CUDACHECK(cudaStreamWaitEvent(targetDeviceBuffers.stream, callerBuffers.events[d], 0));


                    IndexType* d_target_selected_indices = reinterpret_cast<IndexType*>(mr->allocate(sizeof(IndexType) * numSelected, targetDeviceBuffers.stream.getStream()));
                    char* d_target_communicationBuffer = reinterpret_cast<char*>(mr->allocate(sizeof(char) * numSelected * srcRowPitchInBytes, targetDeviceBuffers.stream.getStream()));
                    std::size_t* d_target_numSelected = reinterpret_cast<std::size_t*>(mr->allocate(sizeof(std::size_t), targetDeviceBuffers.stream.getStream()));

                    copy(
                        d_target_numSelected, 
                        deviceId, 
                        splits.d_numSelectedPerGpu.data() + d, 
                        srcDeviceId, 
                        sizeof(size_t), 
                        targetDeviceBuffers.stream
                    );

                    copy(
                        d_target_selected_indices, 
                        deviceId, 
                        splits.d_selectedIndices.data() + numIndices * d, 
                        srcDeviceId, 
                        sizeof(IndexType) * numSelected, 
                        targetDeviceBuffers.stream
                    );

                    copy(
                        d_target_communicationBuffer, 
                        deviceId, 
                        callerBuffers_d_dataCommunicationBuffer.data() + numSelectedPrefixSum[d] * srcRowPitchInBytes, 
                        srcDeviceId, 
                        srcRowPitchInBytes * numSelected, 
                        targetDeviceBuffers.stream
                    );

                    auto scatterIndexIterator = thrust::make_transform_iterator(
                        d_target_selected_indices,
                        MultiGpu2dArrayKernels::SubtractFunctor<IndexType>(h_numRowsPerGpuPrefixSum[d])
                    );

                    gpuArray->scatter(
                        (const T*)d_target_communicationBuffer, 
                        srcRowPitchInBytes, 
                        scatterIndexIterator, 
                        d_target_numSelected,
                        numSelected, 
                        targetDeviceBuffers.stream
                    );

                    mr->deallocate(d_target_selected_indices, sizeof(IndexType) * numSelected, targetDeviceBuffers.stream.getStream());
                    mr->deallocate(d_target_communicationBuffer, sizeof(char) * numSelected * srcRowPitchInBytes, targetDeviceBuffers.stream.getStream());
                    mr->deallocate(d_target_numSelected, sizeof(std::size_t), targetDeviceBuffers.stream.getStream());

                    CUDACHECK(cudaEventRecord(targetDeviceBuffers.event, targetDeviceBuffers.stream));
                }
            }

            for(int d = 0; d < numDistinctGpus; d++){
                const int numSelected = splits.h_numSelectedPerGpu[d];
                if(numSelected > 0){
                    CUDACHECK(cudaStreamWaitEvent(srcStream, handle->deviceBuffers[d].event, 0));
                }
            }
        }
    }

private:
    struct MultiSplitResult{
        rmm::device_uvector<IndexType> d_selectedIndices;
        rmm::device_uvector<std::size_t> d_selectedPositions;
        rmm::device_uvector<std::size_t> d_numSelectedPerGpu;
        std::vector<std::size_t> h_numSelectedPerGpu;
    };

    MultiSplitResult multiSplit(const IndexType* d_indices, std::size_t numIndices, std::size_t* tempmemory, cudaStream_t stream) const{
        auto* srcMR = rmm::mr::get_current_device_resource();

        const int numDistinctGpus = h_numRowsPerGpu.size();

        rmm::device_uvector<IndexType> d_selectedIndices(numIndices * numDistinctGpus, stream, srcMR);
        rmm::device_uvector<std::size_t> d_selectedPositions(numIndices * numDistinctGpus, stream, srcMR);
        rmm::device_uvector<std::size_t> d_numSelectedPerGpu(numDistinctGpus, stream, srcMR);

        CUDACHECK(cudaMemsetAsync(d_numSelectedPerGpu.data(), 0, sizeof(std::size_t) * numDistinctGpus, stream));

        MultiGpu2dArrayKernels::partitionSplitKernel<<<SDIV(numIndices, 128), 128, 0, stream>>>(
            d_selectedIndices.data(),
            d_selectedPositions.data(),
            d_numSelectedPerGpu.data(),
            numDistinctGpus,
            h_numRowsPerGpuPrefixSum.data(),
            numIndices,
            d_indices
        ); CUDACHECKASYNC;

        CUDACHECK(cudaMemcpyAsync(
            tempmemory,
            d_numSelectedPerGpu.data(),
            sizeof(std::size_t) * numDistinctGpus,
            D2H,
            stream
        ));
        CUDACHECK(cudaStreamSynchronize(stream));
        std::vector<std::size_t> h_numSelectedPerGpu(tempmemory, tempmemory + numDistinctGpus);
        
        return MultiSplitResult{
            std::move(d_selectedIndices), 
            std::move(d_selectedPositions), 
            std::move(d_numSelectedPerGpu),
            std::move(h_numSelectedPerGpu)
        };
    }

    void copy(
        void* dst, 
        int dstDevice, 
        const void* src, 
        int srcDevice, 
        size_t count, 
        cudaStream_t stream = 0
    ) const{
        if(dstDevice == srcDevice && dst == src){
            //std::cerr << "copy into same buffer on device " << srcDevice << ". return\n";
            return;
        }else{
            //std::cerr << "copy from device " << srcDevice << " to device " << dstDevice << "\n";
        }

        cudaError_t status = cudaMemcpyPeerAsync(dst, dstDevice, src, srcDevice, count, stream);
        if(status != cudaSuccess){
            cudaGetLastError();
            std::cerr << "dst=" << dst << ", "
            << "dstDevice=" << dstDevice << ", "
            << "src=" << src << ", "
            << "srcDevice=" << srcDevice << ", "
            << "count=" << count << "\n";
        }
        CUDACHECK(status);
    }

    void init(
        size_t numRows, 
        size_t numColumns, 
        std::vector<size_t> memoryLimits, 
        MultiGpu2dArrayLayout layout,
        MultiGpu2dArrayInitMode initMode
    ){
        assert(numRows > 0);
        assert(numColumns > 0);
        assert(memoryLimits.size() == dataDeviceIds.size());

        this->numRows = numRows;
        this->numColumns = numColumns;

        //init gpu arrays

        const int numGpus = dataDeviceIds.size();

        const size_t minbytesPerRow = sizeof(T) * numColumns;
        const size_t rowPitchInBytes = SDIV(minbytesPerRow, getAlignmentInBytes()) * getAlignmentInBytes();

        std::vector<size_t> maxRowsPerGpu(numGpus);

        for(int i = 0; i < numGpus; i++){
            maxRowsPerGpu[i] = memoryLimits[i] / rowPitchInBytes;
        }

        std::vector<size_t> rowsPerGpu(numGpus, 0);
        size_t remaining = numRows;

        if(layout == MultiGpu2dArrayLayout::EvenShare){
            std::size_t divided = numRows / numGpus;

            for(int outer = 0; outer < numGpus; outer++){
                for(int i = 0; i < numGpus; i++){
                    size_t myrows = std::min(maxRowsPerGpu[i] - rowsPerGpu[i], std::min(divided, remaining));
                    rowsPerGpu[i] += myrows;

                    remaining -= myrows;
                }
            }
        }else{
            assert(layout ==MultiGpu2dArrayLayout::FirstFit);

            for(int i = 0; i < numGpus; i++){
                size_t myrows = std::min(maxRowsPerGpu[i], remaining);
                rowsPerGpu[i] = myrows;

                remaining -= myrows;
            }
        }


        if(initMode == MultiGpu2dArrayInitMode::MustFitCompletely){
            if(remaining > 0){
                throw std::invalid_argument("Cannot fit all array elements into provided memory\n");
            }
        }else{
            assert(initMode == MultiGpu2dArrayInitMode::CanDiscardRows);

            this->numRows -= remaining;
        }

        usedDeviceIds.clear();
        std::vector<int> rowsPerUsedGpu;

        for(int i = 0; i < numGpus; i++){
            if(rowsPerGpu[i] > 0){
                cub::SwitchDevice sd(dataDeviceIds[i]);
                auto arrayptr = std::make_unique<Gpu2dArrayManaged<T>>(
                    rowsPerGpu[i], getNumColumns(), getAlignmentInBytes()
                );
                gpuArrays.push_back(std::move(arrayptr));
                usedDeviceIds.push_back(dataDeviceIds[i]);
                rowsPerUsedGpu.push_back(rowsPerGpu[i]);
            }
        }

        const int numUsedDeviceIds = usedDeviceIds.size();

        h_numRowsPerGpu.resize(numUsedDeviceIds);
        h_numRowsPerGpuPrefixSum.resize(numUsedDeviceIds + 1);

        std::copy(rowsPerUsedGpu.begin(), rowsPerUsedGpu.end(), h_numRowsPerGpu.get());
        std::partial_sum(rowsPerUsedGpu.begin(), rowsPerUsedGpu.end(), h_numRowsPerGpuPrefixSum.get() + 1);
        h_numRowsPerGpuPrefixSum[0] = 0;
    }

    bool isReplicatedSingleGpu = false;

    size_t numRows{};
    size_t numColumns{};
    size_t alignmentInBytes{};
    std::vector<std::unique_ptr<Gpu2dArrayManaged<T>>> gpuArrays{};

    std::vector<int> dataDeviceIds;
    std::vector<int> usedDeviceIds;

    HostBuffer<size_t> h_numRowsPerGpu;
    HostBuffer<size_t> h_numRowsPerGpuPrefixSum;
};





#endif
