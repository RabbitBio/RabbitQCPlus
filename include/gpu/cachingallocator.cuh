#ifndef CACHING_ALLOCATOR_CUH
#define CACHING_ALLOCATOR_CUH

#include <thrust/device_malloc_allocator.h>

#include <cub/cub.cuh>

template<class T>
struct ThrustCachingAllocator : thrust::device_malloc_allocator<T> {
    using value_type = T;
    using super_t = thrust::device_malloc_allocator<T>;

    using pointer = typename super_t::pointer;
    using size_type = typename super_t::size_type;
    using reference = typename super_t::reference;
    using const_reference = typename super_t::const_reference;

    int deviceId{};
    cub::CachingDeviceAllocator* cubAllocator{};
    cudaStream_t stream{};

    ThrustCachingAllocator(int deviceId_, cub::CachingDeviceAllocator* cubAllocator_, cudaStream_t stream_)
        : deviceId(deviceId_), cubAllocator(cubAllocator_), stream(stream_){
        
    }

    pointer allocate(size_type n){
        //std::cerr << "alloc" << std::endl;

        T* ptr = nullptr;
        cudaError_t status = cudaSuccess;
        if(cubAllocator){
            status = cubAllocator->DeviceAllocate(deviceId, (void**)&ptr, n * sizeof(T), stream);
        }else{
            status = cudaMalloc(&ptr, n * sizeof(T));
        }
        if(status != cudaSuccess){
            std::cerr << "ThrustCachingAllocator cuda error when allocating " << (n * sizeof(T)) << " bytes: " << cudaGetErrorString(status) << "\n";
            throw std::bad_alloc();
        }
        return thrust::device_pointer_cast(ptr);
    }

    void deallocate(pointer ptr, size_type /*n*/){
    	//std::cerr << "dealloc" << std::endl;
        cudaError_t status = cudaSuccess;
        if(cubAllocator){
            status = cubAllocator->DeviceFree(deviceId, ptr.get());
        }else{
            status = cudaFree(ptr.get());
        }
    	if(status != cudaSuccess){
    		throw std::bad_alloc();
    	}
    }
};

#endif