#ifndef CARE_CPUSEQUENCEHASHER_HPP
#define CARE_CPUSEQUENCEHASHER_HPP


#include <config.hpp>
#include <hpc_helpers.cuh>
#include <sequencehelpers.hpp>

#include <array>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <cassert>

namespace care{

template<class HashValueType>
struct CPUSequenceHasher{

    struct TopSmallestHashResult{
        friend class CPUSequenceHasher;

        TopSmallestHashResult(std::size_t numSmallest) : hashes(numSmallest, std::numeric_limits<HashValueType>::max()) {}

        std::vector<HashValueType> hashes{}; //sorted

        std::size_t size() const noexcept{ return hashes.size(); }
        const HashValueType* cbegin() const noexcept{ return hashes.data(); }
        const HashValueType* cend() const noexcept{ return hashes.data() + size(); }
        const HashValueType* data() const noexcept{ return hashes.data(); }
        HashValueType* begin() noexcept{ return hashes.data(); }
        HashValueType* end() noexcept{ return hashes.data() + size(); }
        HashValueType* data() noexcept{ return hashes.data(); }
        HashValueType& operator[](std::size_t i) noexcept{ assert(i < size()); return hashes[i]; }
        const HashValueType& operator[](std::size_t i) const noexcept{ assert(i < size()); return hashes[i]; }

    private:

        void unique(){
            hashes.erase(
                std::unique(hashes.begin(), hashes.end()),
                hashes.end()
            );
        }

        void insert(HashValueType element){
            //find position of new element
            auto it = std::lower_bound(hashes.begin(), hashes.end(), element);

            if(it != hashes.end()){
                if(*it != element){
                    hashes.pop_back();
                    hashes.insert(it, element);
                }
            }
        }
    };

    TopSmallestHashResult getTopSmallestKmerHashes(
        const unsigned int* sequence,
        const int sequenceLength,
        int kmerLength,
        int numSmallest
    ){
        assert(sizeof(kmer_type) * 8 / 2 >= std::size_t(kmerLength));

        using hasher = hashers::MurmurHash<std::uint64_t>;

        TopSmallestHashResult result(numSmallest);

        SequenceHelpers::forEachEncodedCanonicalKmerFromEncodedSequence(
            sequence,
            sequenceLength,
            kmerLength,
            [&](std::uint64_t kmer, int /*pos*/){
                const auto hashvalue = hasher::hash(kmer);
                result.insert(hashvalue);
            }
        );

        result.unique();

        return result;
    }


    template<class OutputIter>
    OutputIter hashInto(
        OutputIter output,
        const unsigned int* sequence, 
        int sequenceLength, 
        int kmerLength, 
        int numHashFuncs,
        int firstHashFunc
    ){
        constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;
        const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - kmerLength) * 2);

        assert(kmerLength <= maximum_kmer_length);
        
        std::vector<std::uint64_t> hashvalues(numHashFuncs, std::numeric_limits<std::uint64_t>::max());

        if(sequenceLength >= kmerLength){
            SequenceHelpers::forEachEncodedCanonicalKmerFromEncodedSequence(
                sequence,
                sequenceLength,
                kmerLength,
                [&](std::uint64_t kmer, int /*pos*/){
                    using hasher = hashers::MurmurHash<std::uint64_t>;

                    for(int i = 0; i < numHashFuncs; i++){
                        const int hashFuncId = i + firstHashFunc;
                        const auto hashvalue = hasher::hash(kmer + hashFuncId);
                        hashvalues[i] = std::min(hashvalues[i], hashvalue);

                        // if(i == 0){
                        //     std::cerr << pos << ", " << kmer << ", " << (hashvalue & kmer_mask) << "\n";
                        // }
                    }
                }
            );
        }

        return std::transform(hashvalues.begin(), hashvalues.begin() + numHashFuncs, output, [&](auto hash){ return HashValueType(hash & kmer_mask); });
    }

    std::vector<HashValueType> hash(
        const unsigned int* sequence, 
        int sequenceLength, 
        int kmerLength, 
        int numHashFuncs,
        int firstHashFunc
    ){
        std::vector<HashValueType> result(numHashFuncs);

        hashInto(
            result.begin(),
            sequence,
            sequenceLength,
            kmerLength,
            numHashFuncs,
            firstHashFunc
        );
        
        return result;
    }

    template<class OutputIter>
    OutputIter hashWindowedInto(
        OutputIter output,
        const unsigned int* sequence, 
        int sequenceLength, 
        int kmerLength, 
        int windowsize,
        int numHashFuncs,
        int firstHashFunc
    ){
        constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;
        const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - kmerLength) * 2);

        assert(kmerLength <= maximum_kmer_length);

        if(kmerLength > windowsize || sequenceLength < kmerLength) return output;

        const int kmersInWindow = windowsize - kmerLength + 1;

        for(int windowBegin = 0, windowId = 0; windowBegin < sequenceLength - kmerLength + 1; windowBegin += kmersInWindow, windowId++){
            std::vector<std::uint64_t> hashvalues(numHashFuncs, std::numeric_limits<std::uint64_t>::max());

            SequenceHelpers::forEachEncodedCanonicalKmerFromEncodedSequence(
                sequence,
                sequenceLength,
                kmerLength,
                windowBegin, 
                windowBegin + kmersInWindow,
                [&](std::uint64_t kmer, int /*pos*/){
                    using hasher = hashers::MurmurHash<std::uint64_t>;

                    for(int i = 0; i < numHashFuncs; i++){
                        const int hashFuncId = i + firstHashFunc;
                        const auto hashvalue = hasher::hash(kmer + hashFuncId);
                        hashvalues[i] = std::min(hashvalues[i], hashvalue);

                        // if(i == 0){
                        //     std::cerr << pos << ", " << kmer << ", " << (hashvalue & kmer_mask) << "\n";
                        // }
                    }
                }
            );

            output = std::transform(hashvalues.begin(), hashvalues.begin() + numHashFuncs, output, [&](auto hash){ return HashValueType(hash & kmer_mask); });
        }
        
        return output;
    }

    std::vector<HashValueType> hashWindowed(
        const unsigned int* sequence, 
        int sequenceLength,
        int kmerLength, 
        int windowsize,
        int numHashFuncs,
        int firstHashFunc
    ){
        const int kmersInWindow = windowsize - kmerLength + 1;
        const int kmersInSequence = sequenceLength - kmerLength + 1;
        const int numWindows = SDIV(kmersInSequence, kmersInWindow);

        std::vector<HashValueType> result(numWindows * numHashFuncs);

        hashWindowedInto(
            result.begin(),
            sequence,
            sequenceLength,
            kmerLength,
            windowsize,
            numHashFuncs,
            firstHashFunc
        );
        
        return result;
    }
};





}


#endif