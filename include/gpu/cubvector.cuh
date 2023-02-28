#ifndef CARE_CUBVECTOR_CUH
#define CARE_CUBVECTOR_CUH

#include <cub/util_allocator.cuh>
#include <hpc_helpers.cuh>

#include <optional>
#include <thread>

namespace care{

    template<class T, int overprovisioningPercent = 0>
    struct CachedDeviceUVector{
    public:
        using value_type = T;
        static_assert(overprovisioningPercent >= 0, "overprovisioningPercent < 0");

        static constexpr size_t getOverprovisionedSize(size_t requiredSize){
            if(overprovisioningPercent <= 0){
                return requiredSize;
            }else{
                const double onePercent = requiredSize / 100.0;
                const size_t extra = onePercent * overprovisioningPercent;
                return requiredSize + std::min(std::size_t(1), extra);
            }            
        }
    private:
        T* data_{};
        size_t size_{};
        size_t capacity_{};
        cub::CachingDeviceAllocator* allocator_{};
    public:
        std::optional<cudaStream_t> previousAllocStream = std::nullopt;

        void setAllocator(cub::CachingDeviceAllocator& alloc){
            allocator_ = &alloc;
        }

        CachedDeviceUVector() {}
        CachedDeviceUVector(cub::CachingDeviceAllocator& alloc) : allocator_(&alloc) {}

        CachedDeviceUVector(size_t size, cudaStream_t stream, cub::CachingDeviceAllocator& alloc)
            : allocator_(&alloc)
        {
            resize(size, stream);
        }

        CachedDeviceUVector(const CachedDeviceUVector& rhs) = delete;
        CachedDeviceUVector& operator=(const CachedDeviceUVector&) = delete;

        CachedDeviceUVector(CachedDeviceUVector&& rhs) noexcept{
            *this = std::move(rhs);
        }

        CachedDeviceUVector& operator=(CachedDeviceUVector&& rhs) noexcept{

            if(data_ != nullptr){
                cudaError_t status = allocator_->DeviceFree(data_);
                throwOnError(status, "operator= DeviceFree");
            }            

            data_ = rhs.data_;
            size_ = rhs.size_;
            capacity_ = rhs.capacity_;
            allocator_ = rhs.allocator_;
            previousAllocStream = rhs.previousAllocStream;

            rhs.data_ = nullptr;
            rhs.size_ = 0;
            rhs.capacity_ = 0;
            rhs.previousAllocStream.reset();

            return *this;
        }

        ~CachedDeviceUVector(){
            destroy();
        }

        friend void swap(CachedDeviceUVector& l, CachedDeviceUVector& r) noexcept{
            using std::swap;

            swap(l.data_, r.data_);
            swap(l.size_, r.size_);
            swap(l.capacity_, r.capacity_);
            swap(l.allocator_, r.allocator_);
        }
        
        void destroy(){
            if(data_ != nullptr){
                cudaError_t status = allocator_->DeviceFree(data_);
                throwOnError(status, "destroy DeviceFree");
            }
            size_ = 0;
            capacity_ = 0;
            data_ = nullptr;
        }

        T& operator[](size_t i){
            return data()[i];
        }

        const T& operator[](size_t i) const{
            return data()[i];
        }

        //return true if reallocation occured
        //memory content is unspecified after operation
        bool reserveUninitialized(size_t newcapacity, cudaStream_t stream){
            cudaStream_t workstream = stream;
            if(capacity_ < newcapacity){
                cudaError_t status = cudaSuccess;

                if(data_ != nullptr){
                    status = allocator_->DeviceFree(data_);
                    throwOnError(status, "reserveUninitialized DeviceFree");
                }
                status = allocator_->DeviceAllocate((void**)&data_, sizeof(T) * newcapacity, workstream);
                throwOnError(status, "reserveUninitialized DeviceAllocate");
                capacity_ = newcapacity;

                if(previousAllocStream.has_value() && stream != *previousAllocStream){
                    //std::cerr << std::this_thread::get_id() << " reserveUninitialized streamchange " << *previousAllocStream << " -> " << stream << "\n";
                } 
                previousAllocStream = stream;

                return true;
            }else{
                return false;
            }            
        }

        //return true if reallocation occured
        //memory content is unspecified after operation
        bool resizeUninitialized(size_t newsize, cudaStream_t stream){
            const size_t newCapacity = getOverprovisionedSize(newsize);
            bool result = reserveUninitialized(newCapacity, stream);            
            size_ = newsize;
            return result;
        }

        //return true if reallocation occured
        bool reserve(size_t newcapacity, cudaStream_t stream){
            cudaStream_t workstream = stream;
            if(capacity_ < newcapacity){
                cudaError_t status = cudaSuccess;

                T* datanew = nullptr;
                status = allocator_->DeviceAllocate((void**)&datanew, sizeof(T) * newcapacity, workstream);
                throwOnError(status, "reserve DeviceAllocate");

                if(data_ != nullptr){
                    cudaError_t status = cudaMemcpyAsync(
                        datanew,
                        data_,
                        sizeof(T) * size_,
                        D2D,
                        workstream
                    );

                    throwOnError(status, "reserve cudaMemcpyAsync");

                    status = allocator_->DeviceFree(data_);

                    throwOnError(status, "reserve DeviceFree");
                }

                data_ = datanew;
                capacity_ = newcapacity;

                if(previousAllocStream.has_value() && stream != *previousAllocStream){
                    //std::cerr << "reserve streamchange\n";
                }
                previousAllocStream = stream;

                return true;
            }else{
                return false;
            }            
        }

        //return true if reallocation occured
        //memory content in range [oldsize, newsize] is unspecified after operation
        bool resize(size_t newsize, cudaStream_t stream){
            const size_t newCapacity = getOverprovisionedSize(newsize);
            bool result = reserve(newCapacity, stream);            
            size_ = newsize;
            return result;
        }

        void clear(){
            size_ = 0;
        }
        

        void erase(T* first, T* last, cudaStream_t stream){
            cudaStream_t workstream = stream;
            assert(first >= begin());
            assert(last <= end());

            if(previousAllocStream.has_value() && stream != *previousAllocStream){
                //std::cerr << "erase not on previousAllocStream\n";
            }

            if(last < end()){
                const std::size_t elementsAfterRangeToErase = end() - last;

                cudaError_t status = cudaMemcpyAsync(
                    first,
                    last,
                    sizeof(T) * elementsAfterRangeToErase,
                    D2D,
                    workstream
                );

                throwOnError(status, "erase cudaMemcpyAsync");
            }

            size_ -= std::distance(first, last);
        }


        //return true if reallocation occured
        bool append(const T* rangeBegin, const T* rangeEnd, cudaStream_t stream){
            cudaStream_t workstream = stream;
            if(previousAllocStream.has_value() && stream != *previousAllocStream){
                //std::cerr << "append not on previousAllocStream\n";
            }

            const std::size_t rangesize = std::distance(rangeBegin, rangeEnd);
            if(rangesize > 0){
                const std::size_t oldsize = size();
                bool realloc = resize(size() + rangesize, stream);
                cudaError_t status = cudaMemcpyAsync(
                    data() + oldsize,
                    rangeBegin,
                    sizeof(T) * rangesize,
                    cudaMemcpyDefault,
                    workstream
                );
                throwOnError(status, "append cudaMemcpyAsync");

                return realloc;
            }
            return false;
        }

        size_t size() const{
            return size_;
        }

        size_t sizeInBytes() const{
            return size() * sizeof(T);
        }

        size_t capacity() const{
            return capacity_;
        }

        size_t capacityInBytes() const{
            return capacity() * sizeof(T);
        }

        T* data() const noexcept{
            return data_;
        }

        T* begin() const noexcept{
            return data();
        }

        T* end() const noexcept{
            return data() + size();
        }

        bool empty() const noexcept{
            return size() == 0;
        }

    private:
        void throwOnError(cudaError_t status, const std::string& info = "") const{
            if(status != cudaSuccess){
                
                std::string msg = cudaGetErrorString(status);
                std::cerr << "CUDA Error: " << msg << " " << info << "\n";
                throw std::runtime_error(msg);
            }
        }
    };

    template<class T, int overprovisioningPercent = 0>
    using CachedDeviceUScalar = CachedDeviceUVector<T, overprovisioningPercent>;

}



#endif