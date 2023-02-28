#ifndef CARE_ASYNC_RESULT_CUH
#define CARE_ASYNC_RESULT_CUH

#include <hpc_helpers.cuh>
#include <gpu/cudaerrorcheck.cuh>

#include <cassert>

#ifdef __CUDACC__

template<class T>
struct AsyncConstBufferWrapper{
    explicit AsyncConstBufferWrapper(const T* buf) : AsyncConstBufferWrapper(buf, nullptr) {}
    AsyncConstBufferWrapper(const T* buf, cudaEvent_t event) : buffer(buf), readyEvent(event) {}

    void wait() const{
        if(readyEvent != nullptr){
            CUDACHECK(cudaEventSynchronize(readyEvent));
        }
    }

    bool ready() const{
        if(readyEvent != nullptr){
            cudaError_t status = cudaEventQuery(readyEvent);
            assert(status == cudaSuccess || status == cudaErrorNotReady);
            return status == cudaSuccess;
        }else{
            return true;
        }
    }

    //work submitted to stream after linkStream will not be processed until buffer is ready
    void linkStream(cudaStream_t stream) const {
        if(readyEvent != nullptr){
            CUDACHECK(cudaStreamWaitEvent(stream, readyEvent, 0));
        }
    }

    const T* data() const noexcept{
        return buffer;
    }
private:
    const T* buffer = nullptr;
    cudaEvent_t readyEvent = nullptr;
};

template<class T>
AsyncConstBufferWrapper<T> makeAsyncConstBufferWrapper(const T* data, cudaEvent_t event = nullptr){
    return AsyncConstBufferWrapper<T>(data, event);
}





#endif




#endif