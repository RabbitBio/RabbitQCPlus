#ifndef CARE_CONFIG_HPP
#define CARE_CONFIG_HPP

#include <cstdint>
#include <type_traits>


namespace care{

    //unsigned integral type large enough to enumerate all reads. Number of reads must be < numeric_limits<read_number>::max().
    using read_number = std::uint32_t;

    static_assert(std::is_integral<read_number>::value, "read_number must be integral.");
    static_assert(std::is_unsigned<read_number>::value, "read_number must be unsigned.");

    //unsigned integral type of a kmer in the hash map.
    using kmer_type = std::uint64_t;

    static_assert(std::is_integral<kmer_type>::value, "kmer_type must be integral.");
    static_assert(std::is_unsigned<kmer_type>::value, "kmer_type must be unsigned.");

    //unsigned integral type to enumerate hits in a hash table
    using BucketSize = std::uint16_t;
    static_assert(std::is_integral<BucketSize>::value, "BucketSize must be integral.");
    static_assert(std::is_unsigned<BucketSize>::value, "BucketSize must be unsigned.");

    //maximum number of minhash maps
    constexpr int maximum_number_of_maps = 48;
    static_assert(maximum_number_of_maps > 0, "");


    #define MINHASHER_CLEAR_OVEROCCUPIED_BUCKETS
    #define MINHASHER_CLEAR_UNDEROCCUPIED_BUCKETS
    constexpr int MINHASHER_MIN_VALUES_PER_KEY = 2;


    //At least gpuReadStorageHeadroomPerGPU bytes per GPU will not be used by gpuReadStorage
    constexpr std::size_t gpuReadStorageHeadroomPerGPU = std::size_t(1) << 30;


//##################################################

    template<class T> struct max_k;
    template<> struct max_k<std::uint8_t>{static constexpr int value = 4;};
    template<> struct max_k<std::uint16_t>{static constexpr int value = 8;};
    template<> struct max_k<std::uint32_t>{static constexpr int value = 16;};
    template<> struct max_k<std::uint64_t>{static constexpr int value = 32;};
}

#endif
