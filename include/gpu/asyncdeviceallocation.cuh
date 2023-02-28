#ifndef CARE_ASYNCDEVICEALLOCATION_CUH
#define CARE_ASYNCDEVICEALLOCATION_CUH

#include <gpu/cudaerrorcheck.cuh>

namespace care{

    class AsyncDeviceAllocation{
    private:
        void* ptr{};
        std::size_t size{};
        cudaStream_t stream{};
    public:
        AsyncDeviceAllocation() = default;
        AsyncDeviceAllocation(std::size_t size_, cudaStream_t stream_) : size(size_), stream(stream_){
            CUDACHECK(cudaMallocAsync(&ptr, size, stream));
        }

        AsyncDeviceAllocation(const AsyncDeviceAllocation&) = delete;
        AsyncDeviceAllocation(AsyncDeviceAllocation&&) = delete;
        AsyncDeviceAllocation& operator=(const AsyncDeviceAllocation&) = delete;
        AsyncDeviceAllocation& operator=(AsyncDeviceAllocation&&) = delete;

        ~AsyncDeviceAllocation(){
            if(ptr != nullptr){
                destroy(stream);
            }
        }

        void destroy(cudaStream_t stream_){
            CUDACHECK(cudaFreeAsync(ptr, stream_));
            ptr = nullptr;
            size = 0;
        }

        template<class T = void>
        T* get() const noexcept{
            return reinterpret_cast<T*>(ptr);
        }
    };

}

#endif