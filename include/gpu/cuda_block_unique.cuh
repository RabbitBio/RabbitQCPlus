#ifndef CARE_CUDA_BLOCK_UNIQUE_CUH
#define CARE_CUDA_BLOCK_UNIQUE_CUH

#include <cub/cub.cuh>

namespace care{
namespace gpu{


    template<class T, int blocksize, int elemsPerThread>
    struct BlockUnique{
    private:

        using BlockDiscontinuity = cub::BlockDiscontinuity<T, blocksize>;
        using BlockScan = cub::BlockScan<int, blocksize>; 

    public:
        struct TempStorage{
            union{
                union{
                    typename BlockDiscontinuity::TempStorage discontinuity;
                    typename BlockScan::TempStorage scan;
                } cubtemp;

            } data;
        };


        __device__
        BlockUnique(TempStorage& temp) : temp_storage(temp){}

        //return size of output
        __device__
        int execute(T (&input)[elemsPerThread], int validInputSize, T* output){
            constexpr bool isFirstChunk = true;
            return execute_impl<isFirstChunk>(input, validInputSize, T{}, output);
        }

        //return size of output
        __device__
        int execute(T (&input)[elemsPerThread], int validInputSize, T lastElementOfPreviousChunk, T* output){
            constexpr bool isFirstChunk = false;
            return execute_impl<isFirstChunk>(input, validInputSize, lastElementOfPreviousChunk, output);
        }

    private:

        TempStorage& temp_storage;

        template<bool isFirstChunk>
        __device__
        int execute_impl(T (&input)[elemsPerThread], int validInputSize, T lastElementOfPreviousChunk, T* output){
            int prefixsum[elemsPerThread];
            int head_flags[elemsPerThread];

            if(isFirstChunk){
                BlockDiscontinuity(temp_storage.data.cubtemp.discontinuity).FlagHeads(
                    head_flags, 
                    input, 
                    cub::Inequality()
                );
                __syncthreads();
            }else{
                BlockDiscontinuity(temp_storage.data.cubtemp.discontinuity).FlagHeads(
                    head_flags, 
                    input, 
                    cub::Inequality(),
                    lastElementOfPreviousChunk
                );
                __syncthreads();
            }

            #pragma unroll
            for(int i = 0; i < elemsPerThread; i++){
                if(threadIdx.x * elemsPerThread + i >= validInputSize){
                    head_flags[i] = 0;
                }
            }

            int numSelected = 0;
            BlockScan(temp_storage.data.cubtemp.scan).ExclusiveSum(head_flags, prefixsum, numSelected);
            __syncthreads();

            #pragma unroll
            for(int i = 0; i < elemsPerThread; i++){
                if(threadIdx.x * elemsPerThread + i < validInputSize){
                    if(head_flags[i]){
                        output[prefixsum[i]] = input[i];
                    }
                }
            }
            return numSelected;
        }

    };


} //namespace gpu
} //namespace care




#endif