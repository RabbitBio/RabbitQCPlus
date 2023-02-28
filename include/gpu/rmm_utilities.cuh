#ifndef CARE_RMM_UTILITIES_CUH
#define CARE_RMM_UTILITIES_CUH


#include <rmm/device_uvector.hpp>
#include <rmm/device_vector.hpp>
#include <rmm/cuda_stream_view.hpp>
#include <rmm/mr/device/device_memory_resource.hpp>
#include <rmm/mr/device/cuda_async_memory_resource.hpp>
#include <rmm/exec_policy.hpp>
#include <thrust/version.h>

#include <gpu/cudaerrorcheck.cuh>

#include <cassert>
#include <cstdint>
#include <mutex>
#include <map>

#include <cub/cub.cuh>

class CubRMMResource : public rmm::mr::device_memory_resource {
public:

    CubRMMResource(      
        cub::CachingDeviceAllocator& pool
    ) : cubAllocator(&pool){

    }
    ~CubRMMResource() override{
    }

    CubRMMResource(CubRMMResource const&) = delete;
    CubRMMResource(CubRMMResource&&)      = delete;
    CubRMMResource& operator=(CubRMMResource const&) = delete;
    CubRMMResource& operator=(CubRMMResource&&) = delete;

    [[nodiscard]] bool supports_streams() const noexcept override { return true; }

    [[nodiscard]] bool supports_get_mem_info() const noexcept override { return false; }

private:
    cub::CachingDeviceAllocator* cubAllocator;

    void* do_allocate(std::size_t bytes, rmm::cuda_stream_view stream) override{
        void* ptr = nullptr;

        if (bytes > 0) {
            CUDACHECK(cubAllocator->DeviceAllocate(&ptr, bytes, stream.value()));
        }

        return ptr;
    }


    void do_deallocate(void* ptr, std::size_t, rmm::cuda_stream_view) override{
        if (ptr != nullptr) {
            CUDACHECK(cubAllocator->DeviceFree(ptr));
        }
    }

    [[nodiscard]] bool do_is_equal(rmm::mr::device_memory_resource const& other) const noexcept override{
        return dynamic_cast<CubRMMResource const*>(&other) != nullptr;
    }

    [[nodiscard]] std::pair<std::size_t, std::size_t> do_get_mem_info(rmm::cuda_stream_view) const override{
        return std::make_pair(0, 0);
    }
};


template<class T>
void resizeUninitialized(rmm::device_uvector<T>& vec, size_t newsize, rmm::cuda_stream_view stream){
    vec = rmm::device_uvector<T>(newsize, stream, vec.memory_resource());
}

template<class T>
void reserve(rmm::device_uvector<T>& vec, size_t newsize, rmm::cuda_stream_view stream){
    vec.resize(newsize, stream);
}

template<class T>
void clear(rmm::device_uvector<T>& vec, rmm::cuda_stream_view stream){
    vec.resize(0, stream);
}

template<class T>
void destroy(rmm::device_uvector<T>& vec, rmm::cuda_stream_view stream){
    vec.resize(0, stream);
    vec.shrink_to_fit(stream);
}

template<class T>
void erase(rmm::device_uvector<T>& vec, T* first, T* last, rmm::cuda_stream_view stream){
    auto currentend = vec.data() + vec.size();
    assert(first >= vec.data());
    assert(last <= currentend);
    assert(first <= last);


    if(last < currentend){
        const std::size_t elementsAfterRangeToErase = currentend - last;

        CUDACHECK(cudaMemcpyAsync(
            first,
            last,
            sizeof(T) * elementsAfterRangeToErase,
            cudaMemcpyDeviceToDevice,
            stream.value()
        ));

    }

    const std::size_t numErased = std::distance(first, last);
    vec.resize(vec.size() - numErased, stream);
}

template<class T>
void append(rmm::device_uvector<T>& vec, const T* rangeBegin, const T* rangeEnd, rmm::cuda_stream_view stream){
    const std::size_t rangesize = std::distance(rangeBegin, rangeEnd);
    if(rangesize > 0){
        const std::size_t oldsize = vec.size();
        vec.resize(oldsize + rangesize, stream);

        CUDACHECK(cudaMemcpyAsync(
            vec.data() + oldsize,
            rangeBegin,
            sizeof(T) * rangesize,
            cudaMemcpyDefault,
            stream.value()
        ));
    }
}


namespace rmm::mr {

class DeviceCheckResourceAdapter final : public device_memory_resource {
public:
    DeviceCheckResourceAdapter(device_memory_resource* upstream) : upstream_(upstream){
        CUDACHECK(cudaGetDevice(&deviceId));
    }

    DeviceCheckResourceAdapter()                            = default;
    ~DeviceCheckResourceAdapter() override                  = default;
    DeviceCheckResourceAdapter(DeviceCheckResourceAdapter const&) = default;
    DeviceCheckResourceAdapter(DeviceCheckResourceAdapter&&)      = default;
    DeviceCheckResourceAdapter& operator=(DeviceCheckResourceAdapter const&) = default;
    DeviceCheckResourceAdapter& operator=(DeviceCheckResourceAdapter&&) = default;

    bool supports_streams() const noexcept override
    {
        return upstream_->supports_streams();
    }

    bool supports_get_mem_info() const noexcept override
    {
        return upstream_->supports_get_mem_info();
    }

    device_memory_resource* get_upstream() const noexcept{
        return upstream_;
    }

private:
    void checkDeviceId() const{
        int currentDevice;
        CUDACHECK(cudaGetDevice(&currentDevice));
        assert(deviceId == currentDevice);
    }

    void* do_allocate(std::size_t bytes, cuda_stream_view stream) override{
        checkDeviceId();

        return upstream_->allocate(bytes, stream);
    }

    void do_deallocate(void* ptr, std::size_t bytes, cuda_stream_view stream) override{
        checkDeviceId();

        upstream_->deallocate(ptr, bytes, stream);
    }

    bool do_is_equal(device_memory_resource const& other) const noexcept override{
        if (this == &other) { return true; }
        auto const* cast = dynamic_cast<DeviceCheckResourceAdapter const*>(&other);
        if (cast != nullptr) { return upstream_->is_equal(*cast->get_upstream()); }
        return upstream_->is_equal(other);
    }

    std::pair<std::size_t, std::size_t> do_get_mem_info(cuda_stream_view stream) const override {
        return upstream_->get_mem_info(stream);
    }

    int deviceId;
    device_memory_resource* upstream_;
};



class StreamCheckResourceAdapter final : public device_memory_resource {
public:
    StreamCheckResourceAdapter(device_memory_resource* upstream) : upstream_(upstream){}

    StreamCheckResourceAdapter()                            = default;
    ~StreamCheckResourceAdapter() override{
        if(!allocations.empty()){
            std::cerr << "~StreamCheckResourceAdapter: " << allocations.size() << "outstanding allocations\n";
        }
    }
    StreamCheckResourceAdapter(StreamCheckResourceAdapter const&) = default;
    StreamCheckResourceAdapter(StreamCheckResourceAdapter&&)      = default;
    StreamCheckResourceAdapter& operator=(StreamCheckResourceAdapter const&) = default;
    StreamCheckResourceAdapter& operator=(StreamCheckResourceAdapter&&) = default;

    bool supports_streams() const noexcept override
    {
        return upstream_->supports_streams();
    }

    bool supports_get_mem_info() const noexcept override
    {
        return upstream_->supports_get_mem_info();
    }

    device_memory_resource* get_upstream() const noexcept{
        return upstream_;
    }

private:
    struct Payload{
        void* ptr;
        std::size_t bytes;

        bool operator<(const Payload& rhs) const{
            if(ptr < rhs.ptr) return true;
            if(ptr > rhs.ptr) return false;
            return bytes < rhs.bytes;
        }
    };

    void track(void* ptr, std::size_t bytes, cuda_stream_view stream){
        Payload p{ptr, bytes};
        std::lock_guard<std::mutex> lg(mutex);

        auto founditer = std::find_if(allocations.begin(), allocations.end(), [ptr](const auto& pair){
            const auto& payload = pair.first;
            return payload.ptr == ptr;
        });
        //ensure that pointer is not double allocated and inserted
        if(founditer != allocations.end()) throw std::runtime_error("Observed invalid device allocation");
        allocations[p] = stream;
    }

    void checkAndRemove(void* ptr, std::size_t bytes, cuda_stream_view stream){
        Payload p{ptr, bytes};
        std::lock_guard<std::mutex> lg(mutex);
        auto founditer = allocations.find(p);
        if(founditer == allocations.end()) throw std::runtime_error("Observed invalid device allocation");
        if(founditer->second != stream) throw std::runtime_error("Observed invalid device allocation");
        allocations.erase(founditer);
    }
    
    void* do_allocate(std::size_t bytes, cuda_stream_view stream) override{
        void* retval = upstream_->allocate(bytes, stream);
        track(retval, bytes, stream);
        return retval;
    }

    void do_deallocate(void* ptr, std::size_t bytes, cuda_stream_view stream) override{
        checkAndRemove(ptr, bytes, stream);

        upstream_->deallocate(ptr, bytes, stream);
    }

    bool do_is_equal(device_memory_resource const& other) const noexcept override{
        if (this == &other) { return true; }
        auto const* cast = dynamic_cast<StreamCheckResourceAdapter const*>(&other);
        if (cast != nullptr) { return upstream_->is_equal(*cast->get_upstream()); }
        return upstream_->is_equal(other);
    }

    std::pair<std::size_t, std::size_t> do_get_mem_info(cuda_stream_view stream) const override {
        return upstream_->get_mem_info(stream);
    }

    std::mutex mutex;
    std::map<Payload, cuda_stream_view> allocations;
    device_memory_resource* upstream_;
};



class CudaAsyncDefaultPoolResource final : public device_memory_resource {
public:
    CudaAsyncDefaultPoolResource()                            = default;
    ~CudaAsyncDefaultPoolResource() override = default;
    CudaAsyncDefaultPoolResource(CudaAsyncDefaultPoolResource const&) = default;
    CudaAsyncDefaultPoolResource(CudaAsyncDefaultPoolResource&&)      = default;
    CudaAsyncDefaultPoolResource& operator=(CudaAsyncDefaultPoolResource const&) = default;
    CudaAsyncDefaultPoolResource& operator=(CudaAsyncDefaultPoolResource&&) = default;

    bool supports_streams() const noexcept override
    {
        return true;
    }

    bool supports_get_mem_info() const noexcept override
    {
        return false;
    }

    cudaMemPool_t pool_handle() const noexcept{
        int deviceId = 0;
        cudaMemPool_t memPool;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetMemPool(&memPool, deviceId));
        return memPool;
    }

private:
    
    void* do_allocate(std::size_t bytes, cuda_stream_view stream) override{
        void* ptr = nullptr;
        CUDACHECK(cudaMallocAsync(&ptr, bytes, stream.value()));
        return ptr;
    }

    void do_deallocate(void* ptr, std::size_t /*bytes*/, cuda_stream_view stream) override{
        CUDACHECK(cudaFreeAsync(ptr, stream));
    }

    bool do_is_equal(device_memory_resource const& other) const noexcept override{
        if (this == &other) { return true; }
        auto const* cast = dynamic_cast<CudaAsyncDefaultPoolResource const*>(&other);
        return cast != nullptr;
    }

    std::pair<std::size_t, std::size_t> do_get_mem_info(cuda_stream_view /*stream*/) const override {
        return {0,0};
    }
};


class FallbackResourceAdapter final : public device_memory_resource {
public:
    FallbackResourceAdapter(device_memory_resource* upstream, device_memory_resource* fallback) 
        : upstream_(upstream), fallback_(fallback){}

    FallbackResourceAdapter()                            = default;
    ~FallbackResourceAdapter() override = default;
    FallbackResourceAdapter(FallbackResourceAdapter const&) = default;
    FallbackResourceAdapter(FallbackResourceAdapter&&)      = default;
    FallbackResourceAdapter& operator=(FallbackResourceAdapter const&) = default;
    FallbackResourceAdapter& operator=(FallbackResourceAdapter&&) = default;

    bool supports_streams() const noexcept override
    {
        return upstream_->supports_streams() && fallback_->supports_streams();
    }

    bool supports_get_mem_info() const noexcept override
    {
        return upstream_->supports_get_mem_info() && fallback_->supports_get_mem_info();
    }

    device_memory_resource* get_upstream() const noexcept{
        return upstream_;
    }

    device_memory_resource* get_fallback() const noexcept{
        return fallback_;
    }

private:
    struct FallbackPayload{
        void* ptr;
        std::size_t bytes;

        bool operator<(const FallbackPayload& rhs) const{
            if(ptr < rhs.ptr) return true;
            if(ptr > rhs.ptr) return false;
            return bytes < rhs.bytes;
        }
    };

    void trackFallback(void* ptr, std::size_t bytes){
        FallbackPayload p{ptr, bytes};
        std::lock_guard<std::mutex> lg(mutex);
        fallbackallocations.insert(p);
    }

    void* do_allocate(std::size_t bytes, cuda_stream_view stream) override{
        void* retval = nullptr;
        try{
            retval = upstream_->allocate(bytes, stream);
        }catch(rmm::bad_alloc& e){
            retval = fallback_->allocate(bytes, stream);
            trackFallback(retval, bytes);
        }
        return retval;
    }

    void do_deallocate(void* ptr, std::size_t bytes, cuda_stream_view stream) override{
        FallbackPayload p{ptr, bytes};
        std::lock_guard<std::mutex> lg(mutex);
        auto founditer = fallbackallocations.find(p);
        if(founditer == fallbackallocations.end()){
            upstream_->deallocate(ptr, bytes, stream);
        }else{
            fallback_->deallocate(ptr, bytes, stream);
        }
    }

    bool do_is_equal(device_memory_resource const& other) const noexcept override{
        if (this == &other) { return true; }
        auto const* cast = dynamic_cast<FallbackResourceAdapter const*>(&other);
        if (cast != nullptr) { return upstream_->is_equal(*cast->get_upstream()); }
        return upstream_->is_equal(other);
    }

    std::pair<std::size_t, std::size_t> do_get_mem_info(cuda_stream_view stream) const override {
        auto [a1,a2] = upstream_->get_mem_info(stream);
        auto [b1,b2] = fallback_->get_mem_info(stream);
        return std::make_pair(a1+b1, a2+b2);
    }

    device_memory_resource* upstream_;
    device_memory_resource* fallback_;
    std::mutex mutex;
    std::set<FallbackPayload> fallbackallocations;
};



}  // namespace rmm::mr


#endif