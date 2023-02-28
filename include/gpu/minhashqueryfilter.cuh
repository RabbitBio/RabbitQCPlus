#ifndef CARE_MINHASH_QUERY_FILTER_CUH
#define CARE_MINHASH_QUERY_FILTER_CUH

#include <config.hpp>
#include <hpc_helpers.cuh>

#include <gpu/kernels.hpp>
#include <gpu/cubwrappers.cuh>
#include <gpu/cuda_unique.cuh>

#include <rmm/mr/device/per_device_resource.hpp>
#include <rmm/device_uvector.hpp>

#include <cub/cub.cuh>

namespace care{ 
namespace gpu{


    struct GpuMinhashQueryFilter{
        static void keepDistinctAndNotMatching(
            const read_number* d_dontMatchPerSegment,
            cub::DoubleBuffer<read_number>& d_items,
            cub::DoubleBuffer<int>& d_numItemsPerSegment,
            cub::DoubleBuffer<int>& d_numItemsPerSegmentPrefixSum, //numSegments + 1
            int numSegments,
            int numItems,
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
        ){
            if(numItems <= 0) return;
            if(numSegments <= 0) return;

            GpuSegmentedUnique::unique(
                d_items.Current(),
                numItems,
                d_items.Alternate(),
                d_numItemsPerSegment.Alternate(),
                numSegments,
                d_numItemsPerSegmentPrefixSum.Current(),
                d_numItemsPerSegmentPrefixSum.Current() + 1,
                stream,
                mr
            );

            if(d_dontMatchPerSegment != nullptr){
                //remove self read ids (inplace)
                //--------------------------------------------------------------------
                callFindAndRemoveFromSegmentKernel<read_number,128,4>(
                    d_dontMatchPerSegment,
                    d_items.Alternate(),
                    numSegments,
                    d_numItemsPerSegment.Alternate(),
                    d_numItemsPerSegmentPrefixSum.Current(),
                    stream
                );
            }

            CubCallWrapper(mr).cubInclusiveSum(
                d_numItemsPerSegment.Alternate(),
                d_numItemsPerSegmentPrefixSum.Alternate() + 1,
                numSegments,
                stream
            );
            CUDACHECK(cudaMemsetAsync(d_numItemsPerSegmentPrefixSum.Alternate(), 0, sizeof(int), stream));

            //copy final remaining values into contiguous range
            helpers::lambda_kernel<<<numSegments, 128, 0, stream>>>(
                [
                    d_items_in = d_items.Alternate(),
                    d_items_out = d_items.Current(),
                    numSegments,
                    d_numItemsPerSegment = d_numItemsPerSegment.Alternate(),
                    d_offsets = d_numItemsPerSegmentPrefixSum.Current(),
                    d_newOffsets = d_numItemsPerSegmentPrefixSum.Alternate()
                ] __device__ (){

                    for(int s = blockIdx.x; s < numSegments; s += gridDim.x){
                        const int numValues = d_numItemsPerSegment[s];
                        const int inOffset = d_offsets[s];
                        const int outOffset = d_newOffsets[s];

                        for(int c = threadIdx.x; c < numValues; c += blockDim.x){
                            d_items_out[outOffset + c] = d_items_in[inOffset + c];    
                        }
                    }
                }
            ); CUDACHECKASYNC;

            d_numItemsPerSegment.selector++;
            d_numItemsPerSegmentPrefixSum.selector++;

            // helpers::lambda_kernel<<<1,1,0, stream>>>([
            //     d_offsets  = d_numItemsPerSegmentPrefixSum.Current(),
            //     numSegments
            // ] __device__ (){
            //     printf("final offsets before unique\n");
            //     for(int i = 0; i < numSegments+1; i++){
            //         printf("%d ", d_offsets[i]);
            //     }
            //     printf("\n");
            // }); CUDACHECKASYNC;
            // CUDACHECK(cudaDeviceSynchronize());

            // helpers::lambda_kernel<<<1,1,0, stream>>>([
            //     d_numValuesPerSequence = d_numItemsPerSegment.Current(),
            //     numSegments
            // ] __device__ (){
            //     printf("final numValuesPerSequence\n");
            //     for(int i = 0; i < numSegments; i++){
            //         printf("%d ", d_numValuesPerSequence[i]);
            //     }
            //     printf("\n");
            // }); CUDACHECKASYNC;
            // CUDACHECK(cudaDeviceSynchronize());

            // helpers::lambda_kernel<<<1,1,0, stream>>>([
            //     d_values_out = d_items.Current()
            // ] __device__ (){
            //     printf("final values\n");
            //     for(int r = 0; r < 20; r++){
            //         for(int c = 0; c < 16; c++){
            //             printf("%d ", d_values_out[r * 16 + c]);
            //         }
            //         printf("\n");
            //     }
            //     printf("\n");
            // }); CUDACHECKASYNC;
            // CUDACHECK(cudaDeviceSynchronize());
        }
    
        static void keepDistinct(
            cub::DoubleBuffer<read_number>& d_items,
            cub::DoubleBuffer<int>& d_numItemsPerSegment,
            cub::DoubleBuffer<int>& d_numItemsPerSegmentPrefixSum, //numSegments + 1
            int numSegments,
            int numItems,
            cudaStream_t stream,
            rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
        ){
            keepDistinctAndNotMatching(
                nullptr,
                d_items,
                d_numItemsPerSegment,
                d_numItemsPerSegmentPrefixSum,
                numSegments,
                numItems,
                stream,
                mr
            );
        }
    };



}}




#endif