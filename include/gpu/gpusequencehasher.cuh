#ifndef CARE_GPU_SEQUENCE_HASHER_CUH
#define CARE_GPU_SEQUENCE_HASHER_CUH


#include <config.hpp>
#include <hpc_helpers.cuh>
#include <gpu/cudaerrorcheck.cuh>
#include <gpu/cubwrappers.cuh>
#include <gpu/groupmemcpy.cuh>
#include <sequencehelpers.hpp>

#include <cassert>
#include <cstdint>

#include <cub/cub.cuh>

#include <rmm/mr/device/per_device_resource.hpp>
#include <rmm/mr/device/device_memory_resource.hpp>
#include <rmm/device_uvector.hpp>
#include <rmm/device_buffer.hpp>
#include <rmm/device_scalar.hpp>
#include <gpu/rmm_utilities.cuh>

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <cooperative_groups.h>

namespace cg = cooperative_groups;

namespace care{
namespace gpu{
namespace gpusequencehasher{

    struct GetNumKmers{
        int k;
        __host__ __device__
        constexpr GetNumKmers(int k_) : k(k_){}

        __host__ __device__
        constexpr int operator()(int length) const noexcept{
            return (length >= k) ? (length - k + 1) : 0;
        }
    };

    template<int blocksize, class InputIter, class OutputIter>
    __global__
    void minmaxSingleBlockKernel(InputIter begin, int N, OutputIter minmax){
        using value_type_in = typename std::iterator_traits<InputIter>::value_type;
        using value_type_out = typename std::iterator_traits<OutputIter>::value_type;
        static_assert(std::is_same_v<value_type_in, value_type_out>);

        using value_type = value_type_in;

        using BlockReduce = cub::BlockReduce<value_type, blocksize>;
        __shared__ typename BlockReduce::TempStorage temp1;
        __shared__ typename BlockReduce::TempStorage temp2;

        if(blockIdx.x == 0){

            const int tid = threadIdx.x;
            const int stride = blockDim.x;

            value_type myMin = std::numeric_limits<value_type>::max();
            value_type myMax = 0;

            for(int i = tid; i < N; i += stride){
                const value_type val = *(begin + i);
                myMin = min(myMin, val);
                myMax = max(myMax, val);
            }

            myMin = BlockReduce(temp1).Reduce(myMin, cub::Min{});
            myMax = BlockReduce(temp2).Reduce(myMax, cub::Max{});

            if(tid == 0){
                *(minmax + 0) = myMin;
                *(minmax + 1) = myMax;
            }
        }
    }


    template<class T, class InputLengthIter, class InputOffsetIter, class OutputOffsetIter, int groupsize = 32>
    __global__
    void copyRangesKernel(
        const T* __restrict__ inputData,
        InputOffsetIter inputOffsets,
        InputLengthIter inputLengths,
        int numRanges,
        T* __restrict__ outputData,
        OutputOffsetIter outputOffsets
    ){
        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        const int numGroups = (blockDim.x * gridDim.x) / groupsize;

        for(int s = groupId; s < numRanges; s += numGroups){

            const int inputOffset = inputOffsets[s];
            const int outputOffset = outputOffsets[s];
            const int numElementsToCopy = inputLengths[s];

            care::gpu::memcpy<T>(
                group, 
                outputData + outputOffset, 
                inputData + inputOffset, 
                sizeof(T) * numElementsToCopy
            );
        }
    }


    template<class HashValueType>
    __global__
    void minhashSignatures3264Kernel(
        HashValueType* __restrict__ signatures,
        std::size_t signaturesRowPitchElements,
        bool* __restrict__ valid,
        const unsigned int* __restrict__ sequences2Bit,
        std::size_t sequenceRowPitchElements,
        int numSequences,
        const int* __restrict__ sequenceLengths,
        int k,
        int numHashFuncs,
        const int* __restrict__ hashFunctionNumbers
    ){
                
        constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;
        const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);

        const int tid = threadIdx.x + blockIdx.x * blockDim.x;

        if(tid < numSequences * numHashFuncs){
            const int mySequenceIndex = tid / numHashFuncs;
            const int myNumHashFunc = tid % numHashFuncs;
            const int hashFuncId = hashFunctionNumbers[myNumHashFunc];
            assert(hashFuncId < 64);

            const unsigned int* const mySequence = sequences2Bit + mySequenceIndex * sequenceRowPitchElements;
            const int myLength = sequenceLengths[mySequenceIndex];

            HashValueType* const mySignature = signatures + mySequenceIndex * signaturesRowPitchElements;
            bool* const myValid = valid + mySequenceIndex * numHashFuncs;

            if(myLength >= k){
                std::uint64_t minHashValue = std::numeric_limits<std::uint64_t>::max();

                SequenceHelpers::forEachEncodedCanonicalKmerFromEncodedSequence(
                    mySequence,
                    myLength,
                    k,
                    [&](std::uint64_t kmer, int /*pos*/){
                        using hasher = hashers::MurmurHash<std::uint64_t>;

                        const auto hashvalue = hasher::hash(kmer + hashFuncId);
                        minHashValue = min(minHashValue, hashvalue);
                    }
                );

                mySignature[myNumHashFunc] = HashValueType(minHashValue & kmer_mask);
                myValid[myNumHashFunc] = true;
            }else{
                mySignature[myNumHashFunc] = std::numeric_limits<HashValueType>::max();
                myValid[myNumHashFunc] = false;
            }
        }
    }


    template<class ConstBeginOffsetsIter>
    __global__
    void makeKmersKernel(
        kmer_type* __restrict__ kmersoutput,
        ConstBeginOffsetsIter outputBeginOffsets,
        const unsigned int* __restrict__ sequences2Bit,
        std::size_t sequenceRowPitchElements,
        int numSequences,
        const int* __restrict__ sequenceLengths,
        int k
    ){
        assert(sizeof(kmer_type) * 8 / 2 >= k);

        //constexpr int blocksize = 128;
        constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;
        const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);

        for(int s = blockIdx.x; s < numSequences; s += gridDim.x){
            //compute kmers of sequences[s]

            const auto kmerOffset = outputBeginOffsets[s];

            const unsigned int* const mySequence = sequences2Bit + s * sequenceRowPitchElements;
            const int myLength = sequenceLengths[s];

            const int numKmers = (myLength >= k) ? (myLength - k + 1) : 0;
            // if(threadIdx.x == 0){
            //     printf("s %d, length %d, numkmers %d\n", s, myLength, numKmers);
            // }

            for(int i = threadIdx.x; i < numKmers; i += blockDim.x){
                //compute kmer i

                const int firstIntIndex = i / 16;
                const int secondIntIndex = (i + k - 1) / 16;
                std::uint64_t kmer = 0;
                if(firstIntIndex == secondIntIndex){
                    const std::uint64_t firstInt = mySequence[firstIntIndex];
                    kmer = (firstInt >> 2*(16 - (i+k)%16)) & kmer_mask;
                }else{
                    const std::uint64_t firstInt = mySequence[firstIntIndex];
                    const std::uint64_t secondInt = mySequence[secondIntIndex];
                    const int basesInFirst = 16 - (i % 16);
                    const int basesInSecond = k - basesInFirst;
                    kmer = ((firstInt << 2*basesInSecond) | (secondInt >> 2*(16 - (i+k)%16))) & kmer_mask;
                }

                kmersoutput[kmerOffset + i] = kmer_type(kmer);

                // if(s == 0){
                //     printf("i %d, kmer %lu\n", i, kmer_type(kmer));
                // }
            }
        }
    }


    template<bool debug, class HashValueType, class ConstBeginOffsetsIter>
    __global__
    void getKmerHashes(
        HashValueType* __restrict__ kmerhashesoutput,
        ConstBeginOffsetsIter outputBeginOffsets,
        const unsigned int* __restrict__ sequences2Bit,
        std::size_t sequenceRowPitchElements,
        int numSequences,
        const int* __restrict__ sequenceLengths,
        int k
    ){
        assert(sizeof(kmer_type) * 8 / 2 >= k);

        using hasher = hashers::MurmurHash<std::uint64_t>;

        constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;
        const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);
        const int rcshiftamount = (maximum_kmer_length - k) * 2;

        for(int s = blockIdx.x; s < numSequences; s += gridDim.x){

            const auto outputOffset = outputBeginOffsets[s];
            const unsigned int* const mySequence = sequences2Bit + s * sequenceRowPitchElements;
            const int myLength = sequenceLengths[s];

            const int numKmers = (myLength >= k) ? (myLength - k + 1) : 0;

            if constexpr (debug){
                printf("gpu s = %d\n", s);
            }

            for(int i = threadIdx.x; i < numKmers; i += blockDim.x){
                const std::uint64_t kmer = SequenceHelpers::getEncodedKmerFromEncodedSequence(mySequence, k, i);
                const std::uint64_t rc_kmer = SequenceHelpers::reverseComplementInt2Bit(kmer) >> rcshiftamount;
                const auto smallest = min(kmer, rc_kmer);
                kmerhashesoutput[outputOffset + i] = hasher::hash(smallest) & kmer_mask;

                // if(s < 5 && threadIdx.x == 0){
                //     printf("s = %d, kmer %lu, rckmer %lu, smallestkmer %lu, hash: %lu\n", s, kmer, rc_kmer, smallest, kmerhashesoutput[outputOffset + i]);
                // }

                if constexpr (debug){
                    //printf("(%lu, %lu), %lu , %lu , %lu : %lu\n", kmer, rc_kmer, kmer & kmer_mask, rc_kmer & kmer_mask, smallest, (hasher::hash(smallest) & kmer_mask));
                    printf("%lu , %lu , %lu : %lu\n", kmer & kmer_mask, rc_kmer & kmer_mask, smallest, (hasher::hash(smallest) & kmer_mask));
                }
            }


        }
    }


    template<class HashValueType>
    __global__
    void hashKmersKernel(
        HashValueType* __restrict__ signatures,
        std::size_t signaturesRowPitchElements,
        bool* __restrict__ isValid,
        const kmer_type* __restrict__ kmers,
        const int* __restrict__ kmerBeginOffsets,
        const int* __restrict__ kmerEndOffsets,
        int k,
        int numSequences,
        int numHashFuncs,
        const int* __restrict__ hashFunctionNumbers
    ){

        using hasher = hashers::MurmurHash<std::uint64_t>;
                
        constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;
        const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);
        const int rcshiftamount = (maximum_kmer_length - k) * 2;

        for(int s = blockIdx.x; s < numSequences; s += gridDim.x){
            const int kmersBegin = kmerBeginOffsets[s];
            const int kmersEnd = kmerEndOffsets[s];
            const int numKmers = kmersEnd - kmersBegin;

            HashValueType* const mySignature = signatures + s * signaturesRowPitchElements;
            bool* const myIsValid = isValid + s * numHashFuncs;

            if(numKmers > 0){
                for(int i = threadIdx.x; i < numHashFuncs; i += blockDim.x){

                    const int hashFuncId = hashFunctionNumbers[i];

                    std::uint64_t minHashValue = std::numeric_limits<std::uint64_t>::max();

                    for(int x = 0; x < numKmers; x++){
                        const std::uint64_t kmer = kmers[kmersBegin + x];
                        const std::uint64_t rc_kmer = SequenceHelpers::reverseComplementInt2Bit(kmer) >> rcshiftamount;
                        const auto smallest = min(kmer, rc_kmer);
                        const auto hashvalue = hasher::hash(smallest + hashFuncId);
                        minHashValue = min(minHashValue, hashvalue);
                    }

                    mySignature[i] = minHashValue;
                    myIsValid[i] = true;
                }
            }else{
                for(int i = threadIdx.x; i < numHashFuncs; i += blockDim.x){
                    myIsValid[i] = false;
                }
            }
        }
    }

} //namespace gpusequencehasher

template<class HashValueType>
struct GPUSequenceHasher{

    struct Result{
        Result(
            int numSequences,
            int numHashFuncs,
            cudaStream_t stream, 
            rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
        ) : d_hashvalues(numSequences * numHashFuncs, stream, mr),
            d_isValid(numSequences * numHashFuncs, stream, mr){

            CUDACHECK(cudaMemsetAsync(d_isValid.data(), 0, sizeof(bool) * d_isValid.size(), stream));
        }

        rmm::device_uvector<HashValueType> d_hashvalues;
        rmm::device_uvector<bool> d_isValid;
    };

    struct TopSmallestHashResult{
        TopSmallestHashResult(
            int numSequences,
            int numHashFuncs,
            cudaStream_t stream, 
            rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
        ) : d_hashvalues(numSequences * numHashFuncs, stream, mr),
            d_numPerSequences(numSequences, stream, mr),
            d_numPerSequencesPrefixSum(numSequences, stream, mr){

        }

        rmm::device_uvector<HashValueType> d_hashvalues;
        rmm::device_uvector<int> d_numPerSequences;
        rmm::device_uvector<int> d_numPerSequencesPrefixSum; //size() = numSequences
    };

    struct ComputedKmers{
        ComputedKmers(
            int numSequences,
            cudaStream_t stream, 
            rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
        ) : d_offsets(numSequences + 1, stream, mr),
            d_kmers(0, stream, mr){

        }

        rmm::device_uvector<int> d_offsets;
        rmm::device_uvector<kmer_type> d_kmers;
    };

    TopSmallestHashResult getTopSmallestKmerHashes(
        const unsigned int* __restrict__ d_sequences2Bit,
        std::size_t sequenceRowPitchElements,
        int numSequences,
        const int* __restrict__ d_sequenceLengths,
        int k,
        int numSmallest,
        cudaStream_t stream,
        rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource(),
        bool debug = false
    ){
        assert(sizeof(kmer_type) * 8 / 2 >= k);
        assert(k > 0);

        CubCallWrapper cub(mr);

        rmm::device_uvector<int> d_offsets(numSequences + 1, stream, mr);
        CUDACHECK(cudaMemsetAsync(d_offsets.data(), 0, sizeof(int), stream));

        auto d_numKmersPerSequence = thrust::make_transform_iterator(
            d_sequenceLengths,
            gpusequencehasher::GetNumKmers{k}
        );

        int h_minmaxNumKmersPerSequence[2];
        rmm::device_uvector<int> d_minmaxNumKmersPerSequence(2, stream, mr);

        gpusequencehasher::minmaxSingleBlockKernel<512><<<1, 512, 0, stream>>>(
            d_numKmersPerSequence,
            numSequences,
            d_minmaxNumKmersPerSequence.data()
        ); 
        CUDACHECKASYNC; 

        cub.cubInclusiveSum(d_numKmersPerSequence, d_offsets.data() + 1, numSequences, stream);
        CUDACHECK(cudaMemcpyAsync(&h_minmaxNumKmersPerSequence[0], d_minmaxNumKmersPerSequence.data(), sizeof(int) * 2, D2H, stream));

        const int totalNumKmers = d_offsets.back_element(stream);
        CUDACHECK(cudaStreamSynchronize(stream));

        //compute kmer hashes

        rmm::device_uvector<HashValueType> d_hashes(totalNumKmers, stream);

        if(debug){
            gpusequencehasher::getKmerHashes<true><<<numSequences, 128, 0, stream>>>(
                d_hashes.data(),
                d_offsets.data(),
                d_sequences2Bit,
                sequenceRowPitchElements,
                numSequences,
                d_sequenceLengths,
                k
            );
            CUDACHECKASYNC;
        }else{
            gpusequencehasher::getKmerHashes<false><<<numSequences, 128, 0, stream>>>(
                d_hashes.data(),
                d_offsets.data(),
                d_sequences2Bit,
                sequenceRowPitchElements,
                numSequences,
                d_sequenceLengths,
                k
            );
            CUDACHECKASYNC;
        }
        

        //make kmer hashes unique per sequence

        constexpr int begin_bit = 0;
        const int end_bit = 2 * k;

        rmm::device_uvector<HashValueType> d_uniqueHashes(d_hashes.size(), stream, mr);
        rmm::device_uvector<int> d_numUniquePerSequence(numSequences, stream, mr);

        // helpers::lambda_kernel<<<1,1,0,stream>>>(
        //     [
        //         d_hashes = d_hashes.data(),
        //         numHashes = d_hashes.size(),
        //         d_offsets = d_offsets.data(),
        //         d_numKmersPerSequence,
        //         numSequences
        //     ] __device__(){
        //         for(int s = 0; s < min(10, numSequences); s++){
        //             printf("s = %d, numKmers %d, offset %d\n", s, d_numKmersPerSequence[s], d_offsets[s]);

        //             for(int i = 0; i < d_numKmersPerSequence[s]; i++){
        //                 printf("%lu, ", d_hashes[d_offsets[s] + i]);
        //             }
        //             printf("\n");
        //         }
        //     }
        // );
        // CUDACHECKASYNC;

        GpuSegmentedUnique::unique(
            d_hashes.data(),
            d_hashes.size(),
            d_uniqueHashes.data(),
            d_numUniquePerSequence.data(),
            numSequences,
            h_minmaxNumKmersPerSequence[1],
            d_offsets.data(),
            d_offsets.data() + 1,
            begin_bit,
            end_bit,
            stream,
            mr
        );

        d_hashes.resize(0, stream);
        d_hashes.shrink_to_fit(stream);

        //copy top hashes to output
        TopSmallestHashResult result(numSequences, numSmallest, stream, mr);

        //compute usable num hashes per sequence
        helpers::lambda_kernel<<<SDIV(numSequences, 128), 128, 0, stream>>>(
            [
                d_numUniquePerSequence = d_numUniquePerSequence.data(),
                d_numPerSequences = result.d_numPerSequences.data(),
                numSequences, 
                numSmallest
            ] __device__ (){
                const int tid = threadIdx.x + blockIdx.x * blockDim.x;
                const int stride = blockDim.x * gridDim.x;

                for(int i = tid; i < numSequences; i += stride){
                    d_numPerSequences[i] = min(numSmallest, d_numUniquePerSequence[i]);
                }
            }
        );
        CUDACHECKASYNC;

        cub.cubExclusiveSum(result.d_numPerSequences.data(), result.d_numPerSequencesPrefixSum.data(), numSequences, stream);

        gpusequencehasher::copyRangesKernel<<<SDIV(numSequences * numSmallest, 128), 128, 0, stream>>>(
            d_uniqueHashes.data(),
            d_offsets.data(),
            result.d_numPerSequences.data(),
            numSequences,
            result.d_hashvalues.data(),
            result.d_numPerSequencesPrefixSum.data()
        );

        CUDACHECKASYNC;

        return result;
    }


    ComputedKmers computeKmers(
        const unsigned int* __restrict__ d_sequences2Bit,
        std::size_t sequenceRowPitchElements,
        int numSequences,
        const int* __restrict__ d_sequenceLengths,
        int k,
        cudaStream_t stream,
        rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
    ){
        assert(sizeof(kmer_type) * 8 / 2 >= k);
        assert(k > 0);

        ComputedKmers result(numSequences, stream, mr);

        CUDACHECK(cudaMemsetAsync(result.d_offsets.data(), 0, sizeof(int), stream));

        auto d_numKmersPerSequence = thrust::make_transform_iterator(
            d_sequenceLengths,
            gpusequencehasher::GetNumKmers{k}
        );

        CubCallWrapper cub(mr);
        cub.cubInclusiveSum(d_numKmersPerSequence, result.d_offsets.data() + 1, numSequences, stream);

        const int totalNumKmers = result.d_offsets.back_element(stream);
        CUDACHECK(cudaStreamSynchronize(stream));

        result.d_kmers.resize(totalNumKmers, stream);

        gpusequencehasher::makeKmersKernel<<<numSequences, 128, 0, stream>>>(
            result.d_kmers.data(),
            result.d_offsets.data(),
            d_sequences2Bit,
            sequenceRowPitchElements,
            numSequences,
            d_sequenceLengths,
            k
        );
        CUDACHECKASYNC;

        return result;
    }


    ComputedKmers computeUniqueKmers(
        const unsigned int* __restrict__ d_sequences2Bit,
        std::size_t sequenceRowPitchElements,
        int numSequences,
        const int* __restrict__ d_sequenceLengths,
        int k,
        cudaStream_t stream,
        rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
    ){
        assert(sizeof(kmer_type) * 8 / 2 >= k);
        assert(k > 0);

        ComputedKmers computedKmers = computeKmers(
            d_sequences2Bit,
            sequenceRowPitchElements,
            numSequences,
            d_sequenceLengths,
            k,
            stream,
            mr
        );

        // helpers::lambda_kernel<<<1,1,0,stream>>>(
        //     [offsets = computedKmers.d_offsets.data(), size = int(computedKmers.d_offsets.size())] __device__ (){
        //         printf("before: ");
        //         for(int i = 0; i < min(size, 20); i++){
        //             printf("%d ", offsets[i]);
        //         }
        //         printf("\n");
        //     }
        // );

        auto d_numKmersPerSequence = thrust::make_transform_iterator(
            d_sequenceLengths,
            gpusequencehasher::GetNumKmers{k}
        );

        int h_minmaxNumKmersPerSequence[2];
        rmm::device_uvector<int> d_minmaxNumKmersPerSequence(2, stream, mr);

        gpusequencehasher::minmaxSingleBlockKernel<512><<<1, 512, 0, stream>>>(
            d_numKmersPerSequence,
            numSequences,
            d_minmaxNumKmersPerSequence.data()
        ); CUDACHECKASYNC; 

        CUDACHECK(cudaMemcpyAsync(&h_minmaxNumKmersPerSequence[0], d_minmaxNumKmersPerSequence.data(), sizeof(int) * 2, D2H, stream));
        CUDACHECK(cudaStreamSynchronize(stream));

        constexpr int begin_bit = 0;
        const int end_bit = 2 * k;

        rmm::device_uvector<kmer_type> d_uniqueKmers(computedKmers.d_kmers.size(), stream, mr);
        rmm::device_uvector<int> d_numUniquePerSequence(numSequences, stream, mr);

        GpuSegmentedUnique::unique(
            computedKmers.d_kmers.data(),
            computedKmers.d_kmers.size(),
            d_uniqueKmers.data(),
            d_numUniquePerSequence.data(),
            numSequences,
            h_minmaxNumKmersPerSequence[1],
            computedKmers.d_offsets.data(),
            computedKmers.d_offsets.data() + 1,
            begin_bit,
            end_bit,
            stream,
            mr
        );

        CubCallWrapper cub(mr);
        cub.cubInclusiveSum(d_numUniquePerSequence.data(), computedKmers.d_offsets.data() + 1, numSequences, stream);

        std::swap(computedKmers.d_kmers, d_uniqueKmers);

        // helpers::lambda_kernel<<<1,1,0,stream>>>(
        //     [offsets = computedKmers.d_offsets.data(), size = int(computedKmers.d_offsets.size())] __device__ (){
        //         printf("after: ");
        //         for(int i = 0; i < min(size, 20); i++){
        //             printf("%d ", offsets[i]);
        //         }
        //         printf("\n");
        //     }
        // );

        return computedKmers;
    }

    Result hashUniqueKmers(
        const unsigned int* __restrict__ d_sequences2Bit,
        std::size_t sequenceRowPitchElements,
        int numSequences,
        const int* __restrict__ d_sequenceLengths,
        int k,
        int numHashFuncs,
        const int* __restrict__ d_hashFunctionNumbers,
        cudaStream_t stream,
        rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
    ){
        assert(sizeof(kmer_type) * 8 / 2 >= k);
        assert(k > 0);

        if(numSequences == 0){
            return Result{0,0,stream, mr};
        }

        ComputedKmers uniqueKmers = computeUniqueKmers(
            d_sequences2Bit,
            sequenceRowPitchElements,
            numSequences,
            d_sequenceLengths,
            k,
            stream,
            mr
        );

        Result result(numSequences, numHashFuncs, stream, mr);

        gpusequencehasher::hashKmersKernel<<<numSequences, 128, 0, stream>>>(
            result.d_hashvalues.data(),
            numHashFuncs,
            result.d_isValid.data(),
            uniqueKmers.d_kmers.data(),
            uniqueKmers.d_offsets.data(),
            uniqueKmers.d_offsets.data() + 1,
            k,
            numSequences,
            numHashFuncs,
            d_hashFunctionNumbers
        ); 
        CUDACHECKASYNC;

        // helpers::lambda_kernel<<<1,1,0,stream>>>(
        //     [vec = result.d_hashvalues.data()] __device__ (){
        //         for(int i = 0; i < 48; i++){
        //             printf("%lu ", vec[i]);
        //         }
        //         printf("\n");
        //     }
        // );

        CUDACHECKASYNC;

        return result;
    }

    Result hash(
        const unsigned int* __restrict__ sequences2Bit,
        std::size_t sequenceRowPitchElements,
        int numSequences,
        const int* __restrict__ sequenceLengths,
        int k,
        int numHashFuncs,
        const int* __restrict__ hashFunctionNumbers,
        cudaStream_t stream,
        rmm::mr::device_memory_resource* mr = rmm::mr::get_current_device_resource()
    ){
        if(numSequences == 0){
            return Result{0,0,stream, mr};
        }
        Result result(numSequences, numHashFuncs, stream, mr);

        dim3 block(128,1,1);
        dim3 grid(SDIV(numHashFuncs * numSequences, block.x),1,1);

        gpusequencehasher::minhashSignatures3264Kernel<<<grid, block, 0, stream>>>(
            result.d_hashvalues.data(),
            numHashFuncs,
            result.d_isValid.data(),
            sequences2Bit,
            sequenceRowPitchElements,
            numSequences,
            sequenceLengths,
            k,
            numHashFuncs,
            hashFunctionNumbers
        ); CUDACHECKASYNC;

        return result;

    }


};


} //namespace gpu
} //namespace care



#endif