#ifndef MY_CUDA_BLOCK_SELECT_CUH
#define MY_CUDA_BLOCK_SELECT_CUH


#include <cub/cub.cuh>

template<class T, int BLOCK_DIM_X>
struct BlockSelect{
private:
    using BlockScan = cub::BlockScan<int, BLOCK_DIM_X>;

public:
    struct TempStorage{
        union{
            typename BlockScan::TempStorage scan;
            //T gather[BLOCK_DIM_X];
        } storage;
    };

    __device__
    BlockSelect(TempStorage& tmp){
        temp_storage = &tmp;
    }

    //blocked arrangement of inputs
    template<int ITEMS_PER_THREAD, class Flag>
    __device__
    int Flagged(T (&inputs)[ITEMS_PER_THREAD], Flag (&flags)[ITEMS_PER_THREAD], T* output, int num_valid){
        int prefixsum[ITEMS_PER_THREAD];
        int numSelected = 0;

        BlockScan(temp_storage->storage.scan).ExclusiveSum(flags, prefixsum, numSelected);

        #pragma unroll
        for(int i = 0; i < ITEMS_PER_THREAD; i++){
            if(threadIdx.x * ITEMS_PER_THREAD + i < num_valid && flags[i] > 0){
                output[prefixsum[i]] = inputs[i];
            }
        }

        return numSelected;
    }


    //blocked arrangement of inputs
    template<int ITEMS_PER_THREAD, class Flag, class Op>
    __device__
    int ForEachFlagged(T (&inputs)[ITEMS_PER_THREAD], Flag (&flags)[ITEMS_PER_THREAD], int num_valid, Op op){
        int prefixsum[ITEMS_PER_THREAD];
        int numSelected = 0;

        BlockScan(temp_storage->storage.scan).ExclusiveSum(flags, prefixsum, numSelected);

        #pragma unroll
        for(int i = 0; i < ITEMS_PER_THREAD; i++){
            if(threadIdx.x * ITEMS_PER_THREAD + i < num_valid && flags[i] > 0){
                op(inputs[i], prefixsum[i]);
            }
        }

        return numSelected;
    }

    //blocked arrangement of inputs
    template<int ITEMS_PER_THREAD, class Flag, class Op>
    __device__
    int ForEachFlaggedPosition(Flag (&flags)[ITEMS_PER_THREAD], int num_valid, Op op){
        int prefixsum[ITEMS_PER_THREAD];
        int numSelected = 0;

        BlockScan(temp_storage->storage.scan).ExclusiveSum(flags, prefixsum, numSelected);

        #pragma unroll
        for(int i = 0; i < ITEMS_PER_THREAD; i++){
            if(threadIdx.x * ITEMS_PER_THREAD + i < num_valid && flags[i] > 0){
                op(threadIdx.x * ITEMS_PER_THREAD + i, prefixsum[i]);
            }
        }

        return numSelected;
    }

    // //blocked arrangement of inputs
    // template<int ITEMS_PER_THREAD, class Flag>
    // int Flagged(T (&inputs)[ITEMS_PER_THREAD] /*inout*/, Flag (&flags)[ITEMS_PER_THREAD], int num_valid){
    //     int prefixsum[ITEMS_PER_THREAD];
    //     int numSelected = 0;

    //     BlockScan(temp_storage.storage.scan).ExclusiveSum(flags, prefixsum, numSelected);

    //     #pragma unroll
    //     for(int i = 0; i < ITEMS_PER_THREAD; i++){
    //         __syncthreads();
    //         if(threadIdx.x * ITEMS_PER_THREAD + i < num_valid && flags[i] > 0){

    //             temp_storage.storage.gather[prefixsum[i]] = inputs[i];
    //         }
    //         __syncthreads();
    //     }

    //     return numSelected;
    // }

private:
    TempStorage* temp_storage;
};

#endif