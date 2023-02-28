#ifndef TWO_DIMENSIONAL_ARRAY_CUH
#define TWO_DIMENSIONAL_ARRAY_CUH

#include <hpc_helpers.cuh>

#include <cassert>

template<class T>
class TwoDimensionalArray{
public:

    using value_type = T;

    TwoDimensionalArray() = default;

    HOSTDEVICEQUALIFIER
    TwoDimensionalArray(
        size_t numRows,
        size_t numColumns,
        size_t rowPitchInBytes,
        T* arraydata
    ) :
        numRows(numRows),
        numColumns(numColumns),
        rowPitchInBytes(rowPitchInBytes),
        arraydata(arraydata){
    }

    //gather rows from array into dest
    HD_WARNING_DISABLE
    template<class Group, class IndexIterator>
    HOSTDEVICEQUALIFIER
    void gather(Group& group, T* __restrict__ dest, size_t destRowPitchInBytes, IndexIterator indices, size_t numIndices) const{
        #ifdef __CUDA_ARCH__
        //if(group.thread_rank() == 0){printf("destRowPitchInBytes %lu, rowPitchInBytes %lu\n", destRowPitchInBytes, rowPitchInBytes);}
        #endif
        if(destRowPitchInBytes % 4 == 0 && rowPitchInBytes % 4 == 0){
            gather_aligned(group, dest, destRowPitchInBytes, indices, numIndices);
        }else{
            gather_unaligned(group, dest, destRowPitchInBytes, indices, numIndices);
        }
    }

    HD_WARNING_DISABLE
    template<class Group, class IndexIterator>
    HOSTDEVICEQUALIFIER
    void gather_unaligned(Group& group, T* __restrict__ dest, size_t destRowPitchInBytes, IndexIterator indices, size_t numIndices) const{
        #ifdef __CUDA_ARCH__
        //if(group.thread_rank() == 0){printf("unaligned gather %lu\n", numIndices);}
        #endif

        const size_t elementsToCopy = numIndices * numColumns;

        const size_t tid = group.thread_rank();
        const size_t stride = group.size();

        for(size_t i = tid; i < elementsToCopy; i += stride){
            const size_t outputRow = i / numColumns;
            const size_t inputRow = indices[outputRow];
            const size_t column = i % numColumns;

            const T value = ((const T*)(((const char*)arraydata) + inputRow * rowPitchInBytes))[column];
            
            ((T*)(((char*)dest) + outputRow * destRowPitchInBytes))[column] = value;
        }
    }

    //gather rows from array into dest
    HD_WARNING_DISABLE
    template<class Group, class IndexIterator>
    HOSTDEVICEQUALIFIER
    void gather_aligned(Group& group, T* __restrict__ dest, size_t destRowPitchInBytes, IndexIterator indices, size_t numIndices) const{
        #ifdef __CUDA_ARCH__
        //if(group.thread_rank() == 0){printf("aligned gather %lu\n", numIndices);}
        #endif

        using CopyType = int;
        assert(destRowPitchInBytes % sizeof(CopyType) == 0);
        assert(rowPitchInBytes % sizeof(CopyType) == 0);

        const size_t tid = group.thread_rank();
        const size_t stride = group.size();


        const size_t copiesPerRow = SDIV(numColumns * sizeof(T), sizeof(CopyType));
        const size_t numIters = numIndices * copiesPerRow;

        for(size_t i = tid; i < numIters; i += stride){
            const size_t outputRow = i / copiesPerRow;
            const size_t inputRow = indices[outputRow];
            const size_t inputCopyElem = i % copiesPerRow;

            const CopyType value = ((const CopyType*)(((const char*)arraydata) + inputRow * rowPitchInBytes))[inputCopyElem];
            
            ((CopyType*)(((char*)dest) + outputRow * destRowPitchInBytes))[inputCopyElem] = value;
        }
    }

    //scatter rows from src into array
    HD_WARNING_DISABLE
    template<class Group, class IndexIterator>
    HOSTDEVICEQUALIFIER
    void scatter(Group& group, const T* __restrict__ src, size_t srcRowPitchInBytes, IndexIterator indices, size_t numIndices){
        if(srcRowPitchInBytes % 4 == 0 && rowPitchInBytes % 4 == 0){
            scatter_aligned(group, src, srcRowPitchInBytes, indices, numIndices);
        }else{
            scatter_unaligned(group, src, srcRowPitchInBytes, indices, numIndices);
        }
    }

    //scatter rows from src into array
    HD_WARNING_DISABLE
    template<class Group, class IndexIterator>
    HOSTDEVICEQUALIFIER
    void scatter_unaligned(Group& group, const T* __restrict__ src, size_t srcRowPitchInBytes, IndexIterator indices, size_t numIndices){
        #ifdef __CUDA_ARCH__
        //if(group.thread_rank() == 0){printf("unaligned scatter %lu\n", numIndices);}
        #endif

        const size_t tid = group.thread_rank();
        const size_t stride = group.size();

        const size_t elementsToCopy = numIndices * numColumns;

        for(size_t i = tid; i < elementsToCopy; i += stride){
            const size_t inputRow = i / numColumns;
            const size_t outputRow = indices[inputRow];
            const size_t column = i % numColumns;
            
            ((T*)(((char*)arraydata) + outputRow * rowPitchInBytes))[column] 
                = ((const T*)(((const char*)src) + inputRow * srcRowPitchInBytes))[column];
        }
    }

    HD_WARNING_DISABLE
    template<class Group, class IndexIterator>
    HOSTDEVICEQUALIFIER
    void scatter_aligned(Group& group, const T* __restrict__ src, size_t srcRowPitchInBytes, IndexIterator indices, size_t numIndices) const{
        #ifdef __CUDA_ARCH__
        //if(group.thread_rank() == 0){printf("aligned scatter %lu\n", numIndices);}
        #endif

        using CopyType = int;
        assert(srcRowPitchInBytes % sizeof(CopyType) == 0);
        assert(rowPitchInBytes % sizeof(CopyType) == 0);

        const size_t tid = group.thread_rank();
        const size_t stride = group.size();


        const size_t copiesPerRow = SDIV(numColumns * sizeof(T), sizeof(CopyType));
        const size_t numIters = numIndices * copiesPerRow;

        for(size_t i = tid; i < numIters; i += stride){
            const size_t inputRow = i / copiesPerRow;
            const size_t outputRow = indices[inputRow];
            const size_t inputCopyElem = i % copiesPerRow;

            ((CopyType*)(((char*)arraydata) + outputRow * rowPitchInBytes))[inputCopyElem] 
                = ((const CopyType*)(((const char*)src) + inputRow * srcRowPitchInBytes))[inputCopyElem];

        }
    }

private:
    size_t numRows;
    size_t numColumns;
    size_t rowPitchInBytes;
    T* arraydata;
};



#endif
