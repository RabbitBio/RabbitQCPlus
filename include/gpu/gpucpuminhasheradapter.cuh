#ifndef CARE_GPUCPUMINHASHERADAPTER_CUH
#define CARE_GPUCPUMINHASHERADAPTER_CUH

#include <gpu/gpuminhasher.cuh>
#include <cpuminhasher.hpp>
#include <minhasherhandle.hpp>
#include <gpu/cudaerrorcheck.cuh>

#include <hpc_helpers.cuh>

#include <cub/cub.cuh>

namespace care{
namespace gpu{


#if 0
class GPUCPUMinhasherAdapter : public CpuMinhasher{
public:
    using Key_t = CpuMinhasher::Key;

public:

    GPUCPUMinhasherAdapter() = default;

    GPUCPUMinhasherAdapter(const GpuMinhasher& gpuMinhasher_, MinhasherHandle gpuHandle_, cudaStream_t stream_, cub::CachingDeviceAllocator& cubAllocator_)
        : gpuHandle(gpuHandle_), gpuMinhasher(&gpuMinhasher_), stream(stream_), cubAllocator(&cubAllocator_){

    }

public: //inherited interface

    MinhasherHandle makeMinhasherHandle() const override{
        return gpuMinhasher->makeMinhasherHandle();
    }

    void destroyHandle(MinhasherHandle& handle) const override{
        gpuMinhasher->destroyHandle(handle);
    }

    void determineNumValues(
        MinhasherHandle& queryHandle,
        const unsigned int* h_sequenceData2Bit,
        std::size_t encodedSequencePitchInInts,
        const int* h_sequenceLengths,
        int numSequences,
        int* h_numValuesPerSequence,
        int& totalNumValues
    ) const override{

        unsigned int* d_sequenceData2Bit = nullptr;
        CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_sequenceData2Bit, sizeof(unsigned int) * encodedSequencePitchInInts * numSequences, stream));

        int* d_sequenceLengths;
        CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_sequenceLengths, sizeof(int) * numSequences, stream));

        int* d_numValuesPerSequence;
        CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_numValuesPerSequence, sizeof(int) * numSequences, stream));

        CUDACHECK(cudaMemcpyAsync(
            d_sequenceData2Bit,
            h_sequenceData2Bit,
            sizeof(unsigned int) * encodedSequencePitchInInts * numSequences,
            H2D,
            stream
        ));

        CUDACHECK(cudaMemcpyAsync(
            d_sequenceLengths,
            h_sequenceLengths,
            sizeof(int) * numSequences,
            H2D,
            stream
        ));

        gpuMinhasher->determineNumValues(
            queryHandle,
            d_sequenceData2Bit,
            encodedSequencePitchInInts,
            d_sequenceLengths,
            numSequences,
            d_numValuesPerSequence,
            totalNumValues,
            stream
        );

        CUDACHECK(cudaMemcpyAsync(
            h_numValuesPerSequence,
            d_numValuesPerSequence,
            sizeof(int) * numSequences,
            D2H,
            stream
        ));

        CUDACHECK(cubAllocator->DeviceFree(d_sequenceData2Bit));
        CUDACHECK(cubAllocator->DeviceFree(d_sequenceLengths));
        CUDACHECK(cubAllocator->DeviceFree(d_numValuesPerSequence));

        CUDACHECK(cudaStreamSynchronize(stream));
    }

    void retrieveValues(
        MinhasherHandle& queryHandle,
        const read_number* h_readIds,
        int numSequences,
        int totalNumValues,
        read_number* h_values,
        int* h_numValuesPerSequence,
        int* h_offsets //numSequences + 1
    ) const override{
        read_number* d_readIds = nullptr;
        if(h_readIds != nullptr){
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_readIds, sizeof(read_number) * numSequences, stream));

            CUDACHECK(cudaMemcpyAsync(
                d_readIds,
                h_readIds,
                sizeof(read_number) * numSequences,
                H2D,
                stream
            ));
        }

        read_number* d_values = nullptr;
        CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_values, sizeof(read_number) * totalNumValues, stream));

        int* d_numValuesPerSequence = nullptr;
        CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_numValuesPerSequence, sizeof(int) * numSequences, stream));

        int* d_offsets = nullptr;
        CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_offsets, sizeof(int) * (numSequences + 1), stream));

        //TODO THIS MEMCPY IS ONLY NECCESSARY BECAUSE OF MY STUPID WARPCORE MINHASHER INTERFACE. It needs to be refactored
        CUDACHECK(cudaMemcpyAsync(
            d_numValuesPerSequence,
            h_numValuesPerSequence,
            sizeof(int) * numSequences,
            H2D,
            stream
        ));

        //TODO THIS MEMCPY IS ONLY NECCESSARY BECAUSE OF MY STUPID WARPCORE MINHASHER INTERFACE. It needs to be refactored
        cudaMemcpyAsync(
            d_numValuesPerSequence,
            h_numValuesPerSequence,
            sizeof(int) * numSequences,
            H2D,
            stream
        ); CUERR;

        gpuMinhasher->retrieveValues(
            queryHandle,
            d_readIds,
            numSequences,
            totalNumValues,
            d_values,
            d_numValuesPerSequence,
            d_offsets, //numSequences + 1
            stream
        );

        CUDACHECK(cudaMemcpyAsync(
            h_values,
            d_values,
            sizeof(read_number) * totalNumValues,
            D2H,
            stream
        ));

        CUDACHECK(cudaMemcpyAsync(
            h_numValuesPerSequence,
            d_numValuesPerSequence,
            sizeof(int) * numSequences,
            D2H,
            stream
        ));

        CUDACHECK(cudaMemcpyAsync(
            h_offsets,
            d_offsets,
            sizeof(int) * (numSequences + 1),
            D2H,
            stream
        ));

        CUDACHECK(cubAllocator->DeviceFree(d_offsets));
        CUDACHECK(cubAllocator->DeviceFree(d_numValuesPerSequence));
        CUDACHECK(cubAllocator->DeviceFree(d_values));
        if(d_readIds != nullptr){
            CUDACHECK(cubAllocator->DeviceFree(d_readIds));
        }
        CUDACHECK(cudaStreamSynchronize(stream));
    }

    // void compact() override{
    //     gpuMinhasher->compact(stream);
    //     cudaStreamSynchronize(stream);
    // }

    MemoryUsage getMemoryInfo() const noexcept override{
        return gpuMinhasher->getMemoryInfo();
    }

    MemoryUsage getMemoryInfo(const MinhasherHandle& queryHandle) const noexcept{
        return gpuMinhasher->getMemoryInfo(queryHandle);
    }

    int getNumResultsPerMapThreshold() const noexcept override{
        return gpuMinhasher->getNumResultsPerMapThreshold();
    }
    
    int getNumberOfMaps() const noexcept override{
        return gpuMinhasher->getNumberOfMaps();
    }

    // void destroy() override{
    //     gpuMinhasher->destroy();
    // }

private:
    MinhasherHandle gpuHandle;
    const GpuMinhasher* gpuMinhasher;
    cudaStream_t stream;
    cub::CachingDeviceAllocator* cubAllocator;

};

#endif

} //namespace gpu
} //namespace care




#endif