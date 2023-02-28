#ifndef SINGLE_GPU_ARRAY_HPP
#define SINGLE_GPU_ARRAY_HPP

#include <hpc_helpers.cuh>
#include <gpu/cudaerrorcheck.cuh>

#include <2darray.hpp>
#include <memorymanagement.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include <memory>

#include <cooperative_groups.h>
#include <cub/cub.cuh>

#include <thrust/iterator/counting_iterator.h>

namespace cg = cooperative_groups;




namespace Gpu2dArrayManagedKernels{

    template<class T, class IndexIterator>
    __global__
    void scatterkernel(TwoDimensionalArray<T> array, const T* __restrict__ src, size_t srcRowPitchInBytes, IndexIterator indices, const size_t* __restrict__ d_numIndices){
        auto gridGroup = cg::this_grid();

        const size_t numIndices = *d_numIndices;

        if(numIndices == 0) return;

        array.scatter(gridGroup, src, srcRowPitchInBytes, indices, numIndices);
    }

    template<class T, class IndexIterator>
    __global__
    void scatterkernel(TwoDimensionalArray<T> array, const T* __restrict__ src, size_t srcRowPitchInBytes, IndexIterator indices, size_t numIndices){
        auto gridGroup = cg::this_grid();

        if(numIndices == 0) return;

        array.scatter(gridGroup, src, srcRowPitchInBytes, indices, numIndices);
    }

    template<class T, class IndexIterator>
    __global__
    void gatherkernel(TwoDimensionalArray<T> array, T* __restrict__ dest, size_t destRowPitchInBytes, IndexIterator indices, const size_t* __restrict__ d_numIndices){
        auto gridGroup = cg::this_grid();

        const size_t numIndices = *d_numIndices;

        if(numIndices == 0) return;

        array.gather(gridGroup, dest, destRowPitchInBytes, indices, numIndices);
    }

    template<class T, class IndexIterator>
    __global__
    void gatherkernel(TwoDimensionalArray<T> array, T* __restrict__ dest, size_t destRowPitchInBytes, IndexIterator indices, size_t numIndices){
        auto gridGroup = cg::this_grid();

        if(numIndices == 0) return;

        array.gather(gridGroup, dest, destRowPitchInBytes, indices, numIndices);
    }

}

template<class T>
class Gpu2dArrayManaged{
public:

    Gpu2dArrayManaged()
    : alignmentInBytes(sizeof(T)){
        CUDACHECK(cudaGetDevice(&deviceId));
    }

    Gpu2dArrayManaged(size_t numRows, size_t numColumns, size_t alignmentInBytes)
    : alignmentInBytes(alignmentInBytes){

        assert(alignmentInBytes > 0);
        //assert(alignmentInBytes % sizeof(T) == 0);

        init(numRows, numColumns);
    }

    ~Gpu2dArrayManaged(){
        destroy();
    }

    Gpu2dArrayManaged(const Gpu2dArrayManaged& rhs)
        : Gpu2dArrayManaged(rhs.numRows, rhs.numColumns, rhs.alignmentInBytes)
    {
        CUDACHECK(cudaMemcpy(arraydata, rhs.arraydata, numRows * rowPitchInBytes, D2D));
    }

    Gpu2dArrayManaged(Gpu2dArrayManaged&& other) noexcept
        : Gpu2dArrayManaged()
    {
        swap(*this, other);
    }

    Gpu2dArrayManaged& operator=(Gpu2dArrayManaged other){
        swap(*this, other);

        return *this;
    }

    friend void swap(Gpu2dArrayManaged& l, Gpu2dArrayManaged& r) noexcept
    {
        std::swap(l.deviceId, r.deviceId);
        std::swap(l.numRows, r.numRows);
        std::swap(l.numColumns, r.numColumns);
        std::swap(l.rowPitchInBytes, r.rowPitchInBytes);
        std::swap(l.arraydata, r.arraydata);
    }

    void init(size_t numRows, size_t numColumns){
        assert(numRows > 0);
        assert(numColumns > 0);

        this->numRows = numRows;
        this->numColumns = numColumns;

        CUDACHECK(cudaGetDevice(&deviceId));

        rowPitchInBytes = computePitch(numColumns, alignmentInBytes);

        if(numRows > 0 && numColumns > 0){

            CUDACHECK(cudaMalloc(&arraydata, numRows * rowPitchInBytes));
            CUDACHECK(cudaMemset(arraydata, 0, numRows * rowPitchInBytes));

        }
    }

    static size_t computePitch(size_t numberOfColumns, size_t alignment) noexcept{
        const size_t minbytesPerRow = sizeof(T) * numberOfColumns;
        return SDIV(minbytesPerRow, alignment) * alignment;
    }

    void destroy(){
        cub::SwitchDevice sd{deviceId};

        CUDACHECK(cudaFree(arraydata));
        
        numRows = 0;
        numColumns = 0;
        alignmentInBytes = 0;
        rowPitchInBytes = 0;
        arraydata = nullptr;
    }

    std::unique_ptr<Gpu2dArrayManaged> makeCopy(int targetDeviceid) const{
        cub::SwitchDevice sd{targetDeviceid};

        auto result = std::make_unique<Gpu2dArrayManaged>();

        result->numRows = numRows;
        result->numColumns = numColumns;
        result->alignmentInBytes = alignmentInBytes;
        result->rowPitchInBytes = rowPitchInBytes;
        result->arraydata = nullptr;

        cudaError_t status = cudaMalloc(&result->arraydata, numRows * rowPitchInBytes);
        if(status != cudaSuccess){
            return nullptr;
        }
        status = cudaMemcpyAsync(result->arraydata, arraydata, numRows * rowPitchInBytes, D2D, cudaStreamPerThread);
        if(status != cudaSuccess){
            return nullptr;
        }
        status = cudaStreamSynchronize(cudaStreamPerThread);

        if(status != cudaSuccess){
            return nullptr;
        }else{
            return result;
        }
    }

    TwoDimensionalArray<T> wrapper() const noexcept{
        TwoDimensionalArray<T> wrapper(numRows, numColumns, rowPitchInBytes, arraydata);

        return wrapper;
    }

    operator TwoDimensionalArray<T>() const{
        return wrapper();
    }

    void print() const{
        T* tmp;
        CUDACHECK(cudaMallocHost(&tmp, numRows * numColumns));
        gather(tmp, numColumns * sizeof(T), thrust::make_counting_iterator<std::size_t>(0), numRows);
        CUDACHECK(cudaDeviceSynchronize());

        for(size_t i = 0; i < numRows; i++){
            for(size_t k = 0; k < numColumns; k++){
                std::cerr << tmp[i * numColumns + k] << " ";
            }
            std::cerr << "\n";
        }

        CUDACHECK(cudaFreeHost(tmp));
    }

    template<class IndexIterator>
    void gather(T* d_dest, size_t destRowPitchInBytes, IndexIterator d_indices, size_t numIndices, cudaStream_t stream = 0) const{
        if(numIndices == 0) return;
        if(getNumRows() == 0) return;

        dim3 block(128, 1, 1);
        dim3 grid(SDIV(numIndices * numColumns, block.x), 1, 1);

        TwoDimensionalArray<T> array = wrapper();

        Gpu2dArrayManagedKernels::gatherkernel<<<grid, block, 0, stream>>>(
            array, 
            d_dest, 
            destRowPitchInBytes, 
            d_indices, 
            numIndices
        );

        CUDACHECKASYNC;
    }

    template<class IndexIterator>
    void gather(T* d_dest, size_t destRowPitchInBytes, IndexIterator d_indices, const size_t* d_numIndices, size_t maxNumIndices, cudaStream_t stream = 0) const{
        if(maxNumIndices == 0) return;
        if(getNumRows() == 0) return;

        dim3 block(128, 1, 1);
        dim3 grid(SDIV(maxNumIndices * numColumns, block.x), 1, 1);

        TwoDimensionalArray<T> array = wrapper();

        Gpu2dArrayManagedKernels::gatherkernel<<<grid, block, 0, stream>>>(
            array, 
            d_dest, 
            destRowPitchInBytes, 
            d_indices, 
            d_numIndices
        );

        CUDACHECKASYNC;
    }

    void gatherContiguous(T* d_dest, size_t destRowPitchInBytes, size_t rowBegin, size_t numRowsToGather, cudaStream_t stream = 0) const{
        assert(rowBegin + numRowsToGather <= getNumRows());
        if(numRowsToGather == 0) return;
        if(getNumRows() == 0) return;

        if(destRowPitchInBytes == getPitch()){
            CUDACHECK(cudaMemcpyAsync(
                d_dest, 
                ((char*)getGpuData()) + getPitch() * rowBegin, 
                getPitch() * numRowsToGather, 
                D2D, 
                stream
            ));
        }else{
            CUDACHECK(cudaMemcpy2DAsync(
                d_dest, 
                destRowPitchInBytes, 
                ((char*)getGpuData()) + getPitch() * rowBegin, 
                getPitch(), 
                getNumColumns() * sizeof(T), 
                numRowsToGather, 
                D2D, 
                stream
            ));
        }
    }

    //d_dest and stream are on device destDeviceId
    void gatherContiguousPeer(int destDeviceId, T* d_dest, size_t destRowPitchInBytes, size_t rowBegin, size_t numRowsToGather, cudaStream_t stream = 0) const{
        assert(rowBegin + numRowsToGather <= getNumRows());
        if(numRowsToGather == 0) return;
        if(getNumRows() == 0) return;

        if(destRowPitchInBytes == getPitch()){
            peercopy(
                d_dest, 
                destDeviceId, 
                ((char*)getGpuData()) + getPitch() * rowBegin, 
                deviceId,
                getPitch() * numRowsToGather, 
                stream
            );
        }else{
            cudaMemcpy3DPeerParms p;
            std::memset(&p, 0, sizeof(cudaMemcpy3DPeerParms));

            p.dstPtr = make_cudaPitchedPtr((void*)d_dest, destRowPitchInBytes, getNumColumns() * sizeof(T), numRowsToGather);
            p.dstDevice = destDeviceId;
            p.srcPtr = make_cudaPitchedPtr(((char*)getGpuData()) + getPitch() * rowBegin, getPitch(), getNumColumns() * sizeof(T), getNumRows() - rowBegin);
            p.srcDevice = deviceId;
            p.extent = make_cudaExtent(getNumColumns() * sizeof(T), numRowsToGather, 1);

            CUDACHECK(cudaMemcpy3DPeerAsync(&p, stream));
        }
    }


    template<class IndexIterator>
    void scatter(const T* d_src, size_t srcRowPitchInBytes, IndexIterator d_indices, size_t numIndices, cudaStream_t stream = 0) const{
        if(numIndices == 0) return;
        if(getNumRows() == 0) return;

        dim3 block(128, 1, 1);
        dim3 grid(SDIV(numIndices * numColumns, block.x), 1, 1);

        TwoDimensionalArray<T> array = wrapper();

        Gpu2dArrayManagedKernels::scatterkernel<<<grid, block, 0, stream>>>(
            array, 
            d_src, 
            srcRowPitchInBytes, 
            d_indices, 
            numIndices
        );

        CUDACHECKASYNC;
    }

    template<class IndexIterator>
    void scatter(const T* d_src, size_t srcRowPitchInBytes, IndexIterator d_indices, const size_t* d_numIndices, size_t maxNumIndices, cudaStream_t stream = 0) const{
        if(maxNumIndices == 0) return;
        if(getNumRows() == 0) return;

        dim3 block(128, 1, 1);
        dim3 grid(SDIV(maxNumIndices * numColumns, block.x), 1, 1);

        TwoDimensionalArray<T> array = wrapper();

        Gpu2dArrayManagedKernels::scatterkernel<<<grid, block, 0, stream>>>(
            array, 
            d_src, 
            srcRowPitchInBytes, 
            d_indices, 
            d_numIndices
        );

        CUDACHECKASYNC;
    }

    void scatterContiguous(const T* d_src, size_t srcRowPitchInBytes, size_t rowBegin, size_t numRowsToScatter, cudaStream_t stream = 0) const{
        assert(rowBegin + numRowsToScatter <= getNumRows());
        if(numRowsToScatter == 0) return;
        if(getNumRows() == 0) return;

        if(srcRowPitchInBytes == getPitch()){
            CUDACHECK(cudaMemcpyAsync(
                ((char*)getGpuData()) + getPitch() * rowBegin, 
                d_src, 
                getPitch() * numRowsToScatter, 
                D2D, 
                stream
            ));
        }else{
            CUDACHECK(cudaMemcpy2DAsync(
                ((char*)getGpuData()) + getPitch() * rowBegin, 
                getPitch(), 
                d_src,
                srcRowPitchInBytes, 
                getNumColumns() * sizeof(T), 
                numRowsToScatter, 
                D2D, 
                stream
            ));
        }
    }

    //d_src and stream are on device sourceDeviceId
    void scatterContiguousPeer(int sourceDeviceId, const T* d_src, size_t srcRowPitchInBytes, size_t rowBegin, size_t numRowsToScatter, cudaStream_t stream = 0) const{
        assert(rowBegin + numRowsToScatter <= getNumRows());
        if(numRowsToScatter == 0) return;
        if(getNumRows() == 0) return;

        if(srcRowPitchInBytes == getPitch()){
            peercopy(
                ((char*)getGpuData()) + getPitch() * rowBegin, 
                deviceId,
                d_src, 
                sourceDeviceId, 
                getPitch() * numRowsToScatter, 
                stream
            );
        }else{
            cudaMemcpy3DPeerParms p;
            std::memset(&p, 0, sizeof(cudaMemcpy3DPeerParms));

            p.dstDevice = deviceId;
            p.dstPtr = make_cudaPitchedPtr(((char*)getGpuData()) + getPitch() * rowBegin, getPitch(), getNumColumns() * sizeof(T), getNumRows() - rowBegin);
            p.srcDevice = sourceDeviceId;
            p.srcPtr = make_cudaPitchedPtr((void*)d_src, srcRowPitchInBytes, getNumColumns() * sizeof(T), numRowsToScatter);
            p.extent = make_cudaExtent(getNumColumns() * sizeof(T), numRowsToScatter, 1);

            CUDACHECK(cudaMemcpy3DPeerAsync(&p, stream));
        }
    }



    int getDeviceId() const noexcept{
        return deviceId;
    }

    size_t getNumRows() const noexcept{
        return numRows;
    }

    size_t getNumColumns() const noexcept{
        return numColumns;
    }

    size_t getAlignmentInBytes() const noexcept{
        return alignmentInBytes;
    }

    T* getGpuData() const noexcept{
        return arraydata;
    }

    size_t getPitch() const noexcept{
        return rowPitchInBytes;
    }

    MemoryUsage getMemoryInfo() const{
        MemoryUsage result{};

        result.host = 0;
        result.device[getDeviceId()] = getPitch() * getNumRows();

        return result;
    }

private:
    void peercopy(
        void* dst, 
        int dstDevice, 
        const void* src, 
        int srcDevice, 
        size_t count, 
        cudaStream_t stream = 0
    ) const{
        cudaError_t status = cudaMemcpyPeerAsync(dst, dstDevice, src, srcDevice, count, stream);
        if(status != cudaSuccess){
            cudaGetLastError();
            std::cerr << "dst=" << dst << ", "
            << "dstDevice=" << dstDevice << ", "
            << "src=" << src << ", "
            << "srcDevice=" << srcDevice << ", "
            << "count=" << count << "\n";
        }
        CUDACHECK(status);
    }

    int deviceId{};
    size_t numRows{};
    size_t numColumns{};
    size_t alignmentInBytes{};
    size_t rowPitchInBytes{};
    T* arraydata{};   
};




#endif