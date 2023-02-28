#ifndef CARE_MINHASHERADAPTERS_HPP
#define CARE_MINHASHERADAPTERS_HPP

#include <cpuminhasher.hpp>
#include <gpu/gpuminhasher.cuh>

#include <config.hpp>

#include <hpc_helpers.cuh>
#include <sharedmutex.hpp>

#include <cub/cub.cuh>

#include <cstdint>

namespace care{

#if 0
class GpuToCpuMinhasherAdapater : public CpuMinhasher{
private:
    struct GpuData{
        static constexpr int overprovisioningPercent = 5;

        template<class T>
        using DeviceBuffer = helpers::SimpleAllocationDevice<T, overprovisioningPercent>;
        
        template<class T>
        using PinnedBuffer = helpers::SimpleAllocationPinnedHost<T, overprovisioningPercent>;

        DeviceBuffer<char> d_temp{};
        CudaStream stream{};
        std::unique_ptr<gpu::GpuMinhasher::QueryHandle> gpuQueryHandle{};
    };

public:

    GpuToCpuMinhasherAdapater(gpu::GpuMinhasher& wrapped, int deviceId_) 
        : gpuMinhasher{&wrapped}, deviceId{deviceId_}{
    }

    CpuMinhasher::QueryHandle makeQueryHandle() const override{
        cub::SwitchDevice sd(deviceId);

        auto gpuData = std::make_unique<GpuData>();
        gpuData->gpuQueryHandle = std::move(std::make_unique<gpu::GpuMinhasher::QueryHandle>(gpuMinhasher->makeQueryHandle()));

        //std::unique_lock<std::shared_mutex> lock(sharedmutex);
        std::unique_lock<SharedMutex> lock(sharedmutex);
        const int handleid = counter++;
        CpuMinhasher::QueryHandle h = constructHandle(handleid);

        gpuDataVector.emplace_back(std::move(gpuData));
        return h;
    }

    void determineNumValues(
        CpuMinhasher::QueryHandle& queryHandle,
        const unsigned int* h_sequenceData2Bit,
        std::size_t encodedSequencePitchInInts,
        const int* h_sequenceLengths,
        int numSequences,
        int* h_numValuesPerSequence,
        int& totalNumValues
    ) const override{

        cub::SwitchDevice sd(deviceId);

        std::size_t temp_allocation_sizes[3]{};        
        temp_allocation_sizes[0] = sizeof(unsigned int) * encodedSequencePitchInInts * numSequences; // d_sequenceData2Bit
        temp_allocation_sizes[1] = sizeof(int) * numSequences; // d_sequenceLengths
        temp_allocation_sizes[2] = sizeof(int) * numSequences; // d_numValuesPerSequence
        
        void* temp_allocations[3]{};
        std::size_t temp_storage_bytes = 0;
        cudaError_t cubstatus = cub::AliasTemporaries(
            nullptr,
            temp_storage_bytes,
            temp_allocations,
            temp_allocation_sizes
        );
        assert(cubstatus == cudaSuccess);       
        
        GpuData* gpuData = getGpuDataFromHandle(queryHandle);

        gpuData->d_temp.resize(temp_storage_bytes);

        cubstatus = cub::AliasTemporaries(
            gpuData->d_temp.data(),
            temp_storage_bytes,
            temp_allocations,
            temp_allocation_sizes
        );

        unsigned int* const d_sequenceData2Bit = static_cast<unsigned int*>(temp_allocations[0]);
        int* const d_sequenceLengths = static_cast<int*>(temp_allocations[1]);
        int* const d_numValuesPerSequence = static_cast<int*>(temp_allocations[2]);

        CUDACHECK(cudaMemcpyAsync(
            d_sequenceData2Bit, 
            h_sequenceData2Bit, 
            sizeof(unsigned int) * encodedSequencePitchInInts * numSequences,
            H2D,
            gpuData->stream
        ));

        CUDACHECK(cudaMemcpyAsync(
            d_sequenceLengths, 
            h_sequenceLengths, 
            sizeof(int) * numSequences,
            H2D,
            gpuData->stream
        ));

        gpuMinhasher->determineNumValues(
            *gpuData->gpuQueryHandle,
            d_sequenceData2Bit,
            encodedSequencePitchInInts,
            d_sequenceLengths,
            numSequences,
            d_numValuesPerSequence,
            totalNumValues,
            gpuData->stream
        );

        CUDACHECK(cudaMemcpyAsync(
            h_numValuesPerSequence, 
            d_numValuesPerSequence, 
            sizeof(int) * numSequences,
            D2H,
            gpuData->stream
        ));

        CUDACHECK(cudaStreamSynchronize(gpuData->stream));
    }

    void retrieveValues(
        CpuMinhasher::QueryHandle& queryHandle,
        const read_number* h_readIds,
        int numSequences,
        int totalNumValues,
        read_number* h_values,
        int* h_numValuesPerSequence,
        int* h_offsets //numSequences + 1
    ) const override{
        cub::SwitchDevice sd(deviceId);

        std::size_t temp_allocation_sizes[4];        
        temp_allocation_sizes[0] = sizeof(read_number) * numSequences; // d_readIds
        temp_allocation_sizes[1] = sizeof(read_number) * totalNumValues; // d_values
        temp_allocation_sizes[2] = sizeof(int) * numSequences; // d_numValuesPerSequence
        temp_allocation_sizes[3] = sizeof(int) * (numSequences + 1); // d_offsets
        
        void* temp_allocations[4];
        std::size_t temp_storage_bytes = 0;
        cudaError_t cubstatus = cub::AliasTemporaries(
            nullptr,
            temp_storage_bytes,
            temp_allocations,
            temp_allocation_sizes
        );
        assert(cubstatus == cudaSuccess);       
        
        GpuData* gpuData = getGpuDataFromHandle(queryHandle);

        gpuData->d_temp.resize(temp_storage_bytes);

        cubstatus = cub::AliasTemporaries(
            gpuData->d_temp.data(),
            temp_storage_bytes,
            temp_allocations,
            temp_allocation_sizes
        );

        read_number* d_readIds = nullptr;
        read_number* const d_values = static_cast<read_number*>(temp_allocations[1]);
        int* const d_numValuesPerSequence = static_cast<int*>(temp_allocations[2]);
        int* const d_offsets = static_cast<int*>(temp_allocations[3]);

        if(h_readIds != nullptr){
            d_readIds = static_cast<read_number*>(temp_allocations[0]);

            CUDACHECK(cudaMemcpyAsync(
                d_readIds, 
                h_readIds, 
                sizeof(read_number) * numSequences,
                H2D,
                gpuData->stream
            ));
        }

        gpuMinhasher->retrieveValues(
            *gpuData->gpuQueryHandle,
            d_readIds,
            numSequences,
            totalNumValues,
            d_values,
            d_numValuesPerSequence,
            d_offsets,
            gpuData->stream
        );

        CUDACHECK(cudaMemcpyAsync(
            h_numValuesPerSequence, 
            d_numValuesPerSequence, 
            sizeof(int) * numSequences,
            D2H,
            gpuData->stream
        ));

        CUDACHECK(cudaMemcpyAsync(
            h_offsets, 
            d_offsets, 
            sizeof(int) * (numSequences + 1),
            D2H,
            gpuData->stream
        ));

        CUDACHECK(cudaStreamSynchronize(gpuData->stream));

        CUDACHECK(cudaMemcpyAsync(
            h_values, 
            d_values, 
            sizeof(read_number) * (h_offsets[numSequences]),
            D2H,
            gpuData->stream
        ));

        CUDACHECK(cudaStreamSynchronize(gpuData->stream));
    }

    void compact() override{
        gpuMinhasher->compact(cudaStreamPerThread);
        CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));
    }

    MemoryUsage getMemoryInfo() const noexcept override{
        return gpuMinhasher->getMemoryInfo();
    }

    MemoryUsage getMemoryInfo(const CpuMinhasher::QueryHandle& queryHandle) const noexcept{
        GpuData* gpuData = getGpuDataFromHandle(queryHandle);

        return gpuMinhasher->getMemoryInfo(*gpuData->gpuQueryHandle);
    }

    int getNumResultsPerMapThreshold() const noexcept override{
        return gpuMinhasher->getNumResultsPerMapThreshold();
    }
    
    int getNumberOfMaps() const noexcept override{
        return gpuMinhasher->getNumberOfMaps();
    }

    void destroy() override{
        gpuMinhasher->destroy();
    }

private:

    GpuData* getGpuDataFromHandle(const CpuMinhasher::QueryHandle& queryHandle) const{
        std::shared_lock<SharedMutex> lock(sharedmutex);

        return gpuDataVector[queryHandle.getId()].get();
    }

    gpu::GpuMinhasher* gpuMinhasher{};
    int deviceId = 0;
    mutable int counter = 0;
    mutable SharedMutex sharedmutex{};

    mutable std::vector<std::unique_ptr<GpuData>> gpuDataVector{};
};

}
#endif


#endif