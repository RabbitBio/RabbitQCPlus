#ifndef CARE_GPU_GROUPMEMCPY_CUH
#define CARE_GPU_GROUPMEMCPY_CUH

#include <cooperative_groups.h>

namespace care{
namespace gpu{

    template<class VectorCopyType, class Group>
    __device__
    void memcpy(Group& group, void* __restrict__ destination, const void* __restrict__ source, size_t bytes){
        const size_t tid = group.thread_rank();
        const size_t stride = group.size();
        char* dest = reinterpret_cast<char*>(destination);
        const char* src = reinterpret_cast<const char*>(source);
                
        //check pointer alignment
        const int mod1 = reinterpret_cast<size_t>(dest) % sizeof(VectorCopyType);
        const int mod2 = reinterpret_cast<size_t>(src) % sizeof(VectorCopyType);

        if(mod1 != mod2){
            //copy byte-wise
            for(size_t i = tid; i < bytes; i += stride){
                reinterpret_cast<char*>(dest)[i] 
                    = reinterpret_cast<const char*>(src)[i];
            }
        }else{
            //copy bytes until dst and src have proper alignment
            if(mod1 > 0){
                const size_t bytesUntilAligned = min(sizeof(VectorCopyType) - mod1, bytes);
                for(size_t i = tid; i < bytesUntilAligned; i += stride){
                    reinterpret_cast<char*>(dest)[i] = reinterpret_cast<const char*>(src)[i];
                }
                dest += bytesUntilAligned;
                src += bytesUntilAligned;
                bytes -= bytesUntilAligned;
            }

            //copy vectorized
            const size_t numTypedCopies = bytes / sizeof(VectorCopyType);        
            for(size_t i = tid; i < numTypedCopies; i += stride){
                reinterpret_cast<VectorCopyType*>(dest)[i] = reinterpret_cast<const VectorCopyType*>(src)[i];
            }

            dest += numTypedCopies * sizeof(VectorCopyType);
            src += numTypedCopies * sizeof(VectorCopyType);
            bytes -= numTypedCopies * sizeof(VectorCopyType);

            //copy leftover bytes
            for(size_t i = tid; i < bytes; i += stride){
                reinterpret_cast<char*>(dest)[i] = reinterpret_cast<const char*>(src)[i];
            }
        }
    }

}
}




#endif

