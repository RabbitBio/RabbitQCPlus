#ifndef CARE_MEMCPY_KERNEL_CUH
#define CARE_MEMCPY_KERNEL_CUH

#include <cuda/std/tuple>
#include <cuda/std/limits>

namespace care{
namespace gpu{

    struct MemcpyParams{
        void* dest{};
        const void* src{};
        std::size_t bytes{};

        MemcpyParams() = default;
        template<class T, class U>
        MemcpyParams(T* dstPointer, const U* srcPointer, std::size_t numBytes)
            : dest(dstPointer), src(srcPointer), bytes(numBytes){}
    };

    template<class VectorCopyType>
    __device__
    void performMemcpy(MemcpyParams params){
        const size_t tid = size_t(threadIdx.x) + size_t(blockIdx.x) * size_t(blockDim.x);
        const size_t stride = size_t(blockDim.x) * size_t(gridDim.x);
        
        //assuming both pointers are aligned to 4 bytes
        
        //check pointer alignment
        const int mod1 = reinterpret_cast<size_t>(params.dest) % sizeof(VectorCopyType);
        const int mod2 = reinterpret_cast<size_t>(params.src) % sizeof(VectorCopyType);
        if(mod1 == 0 && mod2 == 0){
            const size_t numTypedCopies = params.bytes / sizeof(VectorCopyType);
        
            for(size_t i = tid; i < numTypedCopies; i += stride){
                reinterpret_cast<VectorCopyType*>(params.dest)[i] = reinterpret_cast<const VectorCopyType*>(params.src)[i];
            }

            const size_t remaining = params.bytes - numTypedCopies * sizeof(VectorCopyType);
            for(size_t i = tid; i < remaining; i += stride){
                reinterpret_cast<char*>(params.dest)[numTypedCopies * sizeof(VectorCopyType) + i] 
                    = reinterpret_cast<const char*>(params.src)[numTypedCopies * sizeof(VectorCopyType) + i];
            }
        }else{
            for(size_t i = tid; i < params.bytes; i += stride){
                reinterpret_cast<char*>(params.dest)[i] 
                    = reinterpret_cast<const char*>(params.src)[i];
            }
        }
    }

    template<class VectorCopyType, std::size_t I = 0, class... Tp>
    __device__
    typename std::enable_if<I == sizeof...(Tp), void>::type
    performMemcpy(cuda::std::tuple<Tp...>& tuple){
        //end of recursion
    }

    template<class VectorCopyType, std::size_t I = 0, class... Tp>
    __device__
    typename std::enable_if<I < sizeof...(Tp), void>::type
    performMemcpy(cuda::std::tuple<Tp...>& tuple){
        MemcpyParams params = cuda::std::get<I>(tuple);
        performMemcpy<VectorCopyType>(params);
        
        performMemcpy<VectorCopyType, I + 1, Tp...>(tuple);
    }

    template<class CopyType, class MemcpyParamsTuple>
    __global__ 
    void memcpyKernel(MemcpyParamsTuple tuple) {
        performMemcpy<CopyType>(tuple);
    }



}
}


#endif