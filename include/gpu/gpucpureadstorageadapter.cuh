#ifndef CARE_GPUCPUREADSTORAGEADAPTER_CUH
#define CARE_GPUCPUREADSTORAGEADAPTER_CUH

#include <gpu/gpureadstorage.cuh>
#include <cpureadstorage.hpp>
#include <gpu/cudaerrorcheck.cuh>
#include <hpc_helpers.cuh>

#include <cub/cub.cuh>

namespace care{
namespace gpu{

#if 0
    struct GPUCPUReadStorageAdapter : public CpuReadStorage{
    public:

        GPUCPUReadStorageAdapter() = default;

        GPUCPUReadStorageAdapter(const GpuReadStorage& gpuReadStorage_, ReadStorageHandle gpuHandle_, cudaStream_t stream_, cub::CachingDeviceAllocator& cubAllocator_)
            : gpuHandle(gpuHandle_), gpuReadStorage(&gpuReadStorage_), stream(stream_), cubAllocator(&cubAllocator_){

        }

    public: //inherited interface

        void areSequencesAmbiguous(
            bool* result, 
            const read_number* readIds, 
            int numSequences
        ) const override{

            read_number* d_readIds = nullptr;
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_readIds, sizeof(read_number) * numSequences, stream));
            bool* d_result;
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_result, sizeof(bool) * numSequences, stream));

            CUDACHECK(cudaMemcpyAsync(
                d_readIds,
                readIds,
                sizeof(read_number) * numSequences,
                H2D,
                stream
            ));

            gpuReadStorage->areSequencesAmbiguous(
                gpuHandle,
                d_result, 
                d_readIds, 
                numSequences, 
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                result,
                d_result,
                sizeof(bool) * numSequences,
                D2H,
                stream
            ));

            CUDACHECK(cubAllocator->DeviceFree(d_result));
            CUDACHECK(cubAllocator->DeviceFree(d_readIds));

            CUDACHECK(cudaStreamSynchronize(stream));
        }

        void gatherSequences(
            unsigned int* sequence_data,
            size_t outSequencePitchInInts,
            const read_number* readIds,
            int numSequences
        ) const override{
            read_number* d_readIds = nullptr;
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_readIds, sizeof(read_number) * numSequences, stream));
            unsigned int* d_sequence_data;
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_sequence_data, sizeof(unsigned int) * outSequencePitchInInts * numSequences, stream));

            CUDACHECK(cudaMemcpyAsync(
                d_readIds,
                readIds,
                sizeof(read_number) * numSequences,
                H2D,
                stream
            ));

            gpuReadStorage->gatherSequences(
                gpuHandle,
                d_sequence_data,
                outSequencePitchInInts,
                makeAsyncConstBufferWrapper(readIds),
                d_readIds,
                numSequences,
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                sequence_data,
                d_sequence_data,
                sizeof(unsigned int) * outSequencePitchInInts * numSequences,
                D2H,
                stream
            ));

            CUDACHECK(cubAllocator->DeviceFree(d_sequence_data));
            CUDACHECK(cubAllocator->DeviceFree(d_readIds));

            CUDACHECK(cudaStreamSynchronize(stream));
        }

        void gatherQualities(
            char* quality_data,
            size_t out_quality_pitch,
            const read_number* readIds,
            int numSequences
        ) const override{
            read_number* d_readIds = nullptr;
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_readIds, sizeof(read_number) * numSequences, stream));
            char* d_quality_data;
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_quality_data, sizeof(char) * out_quality_pitch * numSequences, stream));

            CUDACHECK(cudaMemcpyAsync(
                d_readIds,
                readIds,
                sizeof(read_number) * numSequences,
                H2D,
                stream
            ));

            gpuReadStorage->gatherQualities(
                gpuHandle,
                d_quality_data,
                out_quality_pitch,
                makeAsyncConstBufferWrapper(readIds),
                d_readIds,
                numSequences,
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                quality_data,
                d_quality_data,
                sizeof(char) * out_quality_pitch * numSequences,
                D2H,
                stream
            ));

            CUDACHECK(cubAllocator->DeviceFree(d_quality_data));
            CUDACHECK(cubAllocator->DeviceFree(d_readIds));

            CUDACHECK(cudaStreamSynchronize(stream));
        }

        void gatherSequenceLengths(
            int* lengths,
            const read_number* readIds,
            int numSequences
        ) const override{
            read_number* d_readIds = nullptr;
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_readIds, sizeof(read_number) * numSequences, stream));
            int* d_lengths;
            CUDACHECK(cubAllocator->DeviceAllocate((void**)&d_lengths, sizeof(int) * numSequences, stream));

            CUDACHECK(cudaMemcpyAsync(
                d_readIds,
                readIds,
                sizeof(read_number) * numSequences,
                H2D,
                stream
            ));

            gpuReadStorage->gatherSequenceLengths(
                gpuHandle,
                d_lengths,
                d_readIds,
                numSequences,    
                stream
            );

            CUDACHECK(cudaMemcpyAsync(
                lengths,
                d_lengths,
                sizeof(int) * numSequences,
                D2H,
                stream
            ));

            CUDACHECK(cubAllocator->DeviceFree(d_lengths));
            CUDACHECK(cubAllocator->DeviceFree(d_readIds));

            CUDACHECK(cudaStreamSynchronize(stream));
        }

        void getIdsOfAmbiguousReads(
            read_number* ids
        ) const override{
            gpuReadStorage->getIdsOfAmbiguousReads(ids);
        }

        std::int64_t getNumberOfReadsWithN() const override{
            return gpuReadStorage->getNumberOfReadsWithN();
        }

        MemoryUsage getMemoryInfo() const override{
            return gpuReadStorage->getMemoryInfo();
        }

        read_number getNumberOfReads() const override{
            return gpuReadStorage->getNumberOfReads();
        }

        bool canUseQualityScores() const override{
            return gpuReadStorage->canUseQualityScores();
        }

        int getSequenceLengthLowerBound() const override{
            return gpuReadStorage->getSequenceLengthLowerBound();
        }

        int getSequenceLengthUpperBound() const override{
            return gpuReadStorage->getSequenceLengthUpperBound();
        }

        bool isPairedEnd() const override{
            return gpuReadStorage->isPairedEnd();
        }

        // void destroy() override{
        //     return gpuReadStorage->destroy();
        // }

    private:
        mutable ReadStorageHandle gpuHandle;
        const GpuReadStorage* gpuReadStorage;
        cudaStream_t stream;
        cub::CachingDeviceAllocator* cubAllocator;
    };

#endif

} //namespace gpu
} //namespace care




#endif