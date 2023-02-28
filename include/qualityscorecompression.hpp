#ifndef CARE_QUALITY_SCORE_COMPRESSION
#define CARE_QUALITY_SCORE_COMPRESSION

#include <hpc_helpers.cuh>
#include <gpu/cudaerrorcheck.cuh>

#include <cstdint>

#ifdef __CUDACC__

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#endif

struct QualityCompressionHelper{
    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int getNumInts(int lengths, int numBitsPerQual){
        return SDIV(numBitsPerQual * lengths, sizeof(unsigned int) * 8);
    }
};


template<int numBits>
struct QualityCompressorConfig;

template<>
struct QualityCompressorConfig<1>{
    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int bitsPerQual() noexcept{
        return 1;
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int binBoundary(int binNumber) noexcept{
        constexpr int boundaries[2]{13, 256};

        return boundaries[binNumber];
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int binQuality(int binNumber) noexcept{
        constexpr char binsQualities[2]{'(', '5'};

        return binsQualities[binNumber];
    }
};

template<>
struct QualityCompressorConfig<2>{
    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int bitsPerQual() noexcept{
        return 2;
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int binBoundary(int binNumber) noexcept{
        constexpr int boundaries[4]{4,13,20,256};

        return boundaries[binNumber];
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int binQuality(int binNumber) noexcept{
        constexpr char binsQualities[4]{'$', '(', '1', '9'};

        return binsQualities[binNumber];
    }
};



template<class Config>
struct QualityCompressorImpl{
    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int bitsPerQual(){
        return Config::bitsPerQual();
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int numBins() noexcept{
        return 1 << bitsPerQual();
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int getNumInts(int lengths){
        return QualityCompressionHelper::getNumInts(lengths, bitsPerQual());
    }

    static void encodeQualityString(unsigned int* out, const char* quality, int length){
        constexpr int ASCII_BASE = 33;

        const int numInts = getNumInts(length);
        constexpr int numQualsPerInts = sizeof(unsigned int) / bitsPerQual() * 8;

        for(int n = 0; n < numInts; n++){

            unsigned int data = 0;
            const int pend = std::min(length, (n+1) * numQualsPerInts);

            for(int p = n * numQualsPerInts; p < pend; p++){
                const int Q = int(quality[p]) - ASCII_BASE;

                for(unsigned int x = 0; x < (unsigned int)numBins(); x++){
                    if(Q <= Config::binBoundary(x)){
                        data = (data << bitsPerQual()) | x;
                        break;
                    }
                }
            }

            if(n == numInts-1){
                const int remaining = (n+1) * numQualsPerInts - length;
                data <<= remaining * bitsPerQual();
            }

            out[n] = data;
        }
    }

    static void decodeQualityToString(char* quality, const unsigned int* encoded, int length){

        const int numInts = getNumInts(length);
        constexpr int numQualsPerInts = sizeof(unsigned int) / bitsPerQual() * 8;

        for(int n = 0; n < numInts; n++){

            unsigned int data = encoded[n];
            const int pend = std::min(length, (n+1) * numQualsPerInts);

            for(int p = n * numQualsPerInts; p < pend; p++){
                const int posInInt = p - n * numQualsPerInts;

                const unsigned int bits = (data >> ((32 - bitsPerQual()) - (posInInt * bitsPerQual()))) & ((unsigned int)numBins() - 1);
                for(unsigned int x = 0; x < 4; x++){
                    if(bits == x){
                        quality[p] = Config::binQuality(x);
                        break;
                    }
                }
            }
        }
    }

    #ifdef __CUDACC__

    template<class Group>
    DEVICEQUALIFIER INLINEQUALIFIER
    static void encodeQualityString(Group& group, unsigned int* out, const char* quality, int length){
        constexpr int ASCII_BASE = 33;

        const int numInts = getNumInts(length);
        constexpr int numQualsPerInts = sizeof(unsigned int) / bitsPerQual() * 8;

        for(int n = group.thread_rank(); n < numInts; n += group.size()){

            unsigned int data = 0;
            const int pend = min(length, (n+1) * numQualsPerInts);

            for(int p = n * numQualsPerInts; p < pend; p++){
                const int Q = int(quality[p]) - ASCII_BASE;

                for(unsigned int x = 0; x < (unsigned int)numBins(); x++){
                    if(Q <= Config::binBoundary(x)){
                        data = (data << bitsPerQual()) | x;
                        break;
                    }
                }
            }

            if(n == numInts-1){
                const int remaining = (n+1) * numQualsPerInts - length;
                data <<= remaining * bitsPerQual();
            }

            out[n] = data;
        }
    }

    template<class Group>
    DEVICEQUALIFIER INLINEQUALIFIER
    static void decodeQualityToString(Group& group, char* quality, const unsigned int* encoded, int length){
        const int numInts = getNumInts(length);
        constexpr int numQualsPerInts = sizeof(unsigned int) / bitsPerQual() * 8;

        const int iters = SDIV(length, 4);

        for(int i = group.thread_rank(); i < iters; i += group.size()){
            __align__(4) char temp[4];

            #pragma unroll
            for(int k = 0; k < 4; k++){
                const int pos = i * 4 + k;
                if(pos < length){
                    const int whichInt = pos / numQualsPerInts;
                    const int whichBit = pos - whichInt * numQualsPerInts;

                    const unsigned int bits = (encoded[whichInt] >> ((32 - bitsPerQual()) - (whichBit * bitsPerQual()))) & ((unsigned int)numBins() - 1);
                    for(unsigned int x = 0; x < 4; x++){
                        if(bits == x){
                            temp[k] = Config::binQuality(x);
                            break;
                        }
                    }
                }
            }

            if(i < iters - 1){
                ((unsigned int*)quality)[i] = *((const unsigned int*)&temp[0]);
            }else{
                #pragma unroll
                for(int k = 0; k < 4; k++){
                    const int pos = i * 4 + k;
                    if(pos < length){
                        quality[pos] = temp[k];
                    }
                }
            }
        }
    }

    #endif
};





template<int numBits>
struct QualityCompressor{

    using Impl = QualityCompressorImpl<QualityCompressorConfig<numBits>>;

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int bitsPerQual(){
        return Impl::bitsPerQual();
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int getNumInts(int length){
        return Impl::getNumInts(length);
    }

    static void encodeQualityString(unsigned int* out, const char* quality, int length){
        Impl::encodeQualityString(out, quality, length);
    }

    static void decodeQualityToString(char* quality, const unsigned int* encoded, int length){
        Impl::decodeQualityToString(quality, encoded, length);
    }

    #ifdef __CUDACC__

    template<class Group>
    DEVICEQUALIFIER INLINEQUALIFIER
    static void encodeQualityString(Group& group, unsigned int* out, const char* quality, int length){
        Impl::encodeQualityString(group, out, quality, length);
    }

    template<class Group>
    DEVICEQUALIFIER INLINEQUALIFIER
    static void decodeQualityToString(Group& group, char* quality, const unsigned int* encoded, int length){
        Impl::decodeQualityToString(group, quality, encoded, length);
    }

    #endif
};


//8-bit "compression" is identity operation, i.e. memcpy
template<>
struct QualityCompressor<8>{

    //8bits for quality scores means leave them as is

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int bitsPerQual(){
        return 8;
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static constexpr int getNumInts(int lengths){
        return QualityCompressionHelper::getNumInts(lengths, bitsPerQual());
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static void encodeQualityString(unsigned int* out, const char* quality, int length){
        memcpy(out, quality, sizeof(char) * length);
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    static void decodeQualityToString(char* quality, const unsigned int* encoded, int length){
        memcpy(quality, encoded, sizeof(char) * length);
    }

    #ifdef __CUDACC__

    template<class Group>
    DEVICEQUALIFIER INLINEQUALIFIER
    static void encodeQualityString(Group& group, unsigned int* out, const char* quality, int length){
        //copy quality to out
        const int full = length / sizeof(int);
        for(int i = group.thread_rank(); i < full; i += group.size()){
            out[i] = ((const unsigned int*)quality)[i];
        }

        const int remaining = length - full * sizeof(int);
        for(int i = group.thread_rank(); i < remaining; i += group.size()){
            ((char*)out)[full * sizeof(int) + i] = quality[full * sizeof(int) + i];
        }
    }

    template<class Group>
    DEVICEQUALIFIER INLINEQUALIFIER
    static void decodeQualityToString(Group& group, char* quality, const unsigned int* encoded, int length){
        //copy encoded to quality
        const int full = length / sizeof(int);
        for(int i = group.thread_rank(); i < full; i += group.size()){
            ((unsigned int*)quality)[i] = encoded[i];
        }

        const int remaining = length - full * sizeof(int);
        for(int i = group.thread_rank(); i < remaining; i += group.size()){
            quality[full * sizeof(int) + i] = ((const char*)encoded)[full * sizeof(int) + i];
        }
    }

    #endif
};


struct QualityCompressorWrapper{
    int numBits = 8;

    QualityCompressorWrapper(int numQualityBits) : numBits(numQualityBits){}

    //8bits for quality scores means leave them as is

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    int bitsPerQual(){
        return numBits;
    }

    HOSTDEVICEQUALIFIER INLINEQUALIFIER
    int getNumInts(int lengths){
        return QualityCompressionHelper::getNumInts(lengths, bitsPerQual());
    }

    void encodeQualityString(unsigned int* out, const char* quality, int length){
        switch(numBits){
            case 1: QualityCompressor<1>::encodeQualityString(out, quality, length); break;
            case 2: QualityCompressor<2>::encodeQualityString(out, quality, length); break;
            case 8: QualityCompressor<8>::encodeQualityString(out, quality, length); break;
            default: assert(false); break;
        }
    }

    void decodeQualityToString(char* quality, const unsigned int* encoded, int length){
        switch(numBits){
            case 1: QualityCompressor<1>::decodeQualityToString(quality, encoded, length); break;
            case 2: QualityCompressor<2>::decodeQualityToString(quality, encoded, length); break;
            case 8: QualityCompressor<8>::decodeQualityToString(quality, encoded, length); break;
            default: assert(false); break;
        }
    }
};


//kernels

#ifdef __CUDACC__

template<int numBits>
struct QualityCompressionGroupSize;

template<> struct QualityCompressionGroupSize<1>{static constexpr int value = 4;};
template<> struct QualityCompressionGroupSize<2>{static constexpr int value = 8;};
template<> struct QualityCompressionGroupSize<8>{static constexpr int value = 32;};

template<int numBits, class LengthIterator>
__global__
void compressQualityScoresKernel(
    unsigned int* __restrict__ compressed, 
    std::size_t compressedPitchInInts,
    const char* __restrict__ quality, 
    std::size_t qualityPitchInBytes,
    LengthIterator lengths,
    int numSequences
){
    constexpr int groupsize = QualityCompressionGroupSize<numBits>::value;
    
    auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
    const int numGroups = (gridDim.x * blockDim.x) / groupsize;
    const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;

    for(int s = groupId; s < numSequences; s += numGroups){
        const char* const myQuality = quality + qualityPitchInBytes * s;
        unsigned int* const myCompressed = compressed + compressedPitchInInts * s;
        const int myLength = lengths[s];

        QualityCompressor<numBits>::encodeQualityString(group, myCompressed, myQuality, myLength);
    }
}



template<int numBits, class LengthIterator>
void callCompressQualityScoresKernel(
    unsigned int* __restrict__ compressed, 
    std::size_t compressedPitchInInts,
    const char* __restrict__ quality, 
    std::size_t qualityPitchInBytes,
    LengthIterator lengths,
    int numSequences,
    cudaStream_t stream
){
    constexpr int groupsize = QualityCompressionGroupSize<numBits>::value;
    constexpr int blocksize = 128;
    constexpr int groupsPerBlock = blocksize / groupsize;

    const int numBlocks = SDIV(numSequences, groupsPerBlock);

    compressQualityScoresKernel<numBits><<<numBlocks, blocksize, 0, stream>>>(
        compressed,
        compressedPitchInInts,
        quality,
        qualityPitchInBytes,
        lengths,
        numSequences
    );
    CUDACHECKASYNC;
}

template<class LengthIterator>
void callCompressQualityScoresKernel(
    unsigned int* __restrict__ compressed, 
    std::size_t compressedPitchInInts,
    const char* __restrict__ quality, 
    std::size_t qualityPitchInBytes,
    LengthIterator lengths,
    int numSequences,
    int numBits,
    cudaStream_t stream
){
    switch(numBits){
        case 1: callCompressQualityScoresKernel<1>(compressed, compressedPitchInInts, quality, qualityPitchInBytes, lengths, numSequences, stream); break;
        case 2: callCompressQualityScoresKernel<2>(compressed, compressedPitchInInts, quality, qualityPitchInBytes, lengths, numSequences, stream); break;
        case 8: callCompressQualityScoresKernel<8>(compressed, compressedPitchInInts, quality, qualityPitchInBytes, lengths, numSequences, stream); break;
        default:
            std::cerr << "callCompressQualityScoresKernel cannot be called with numBits = " << numBits << "\n";
            assert(false);
            break;
    }
}


template<int numBits>
struct QualityDecompressionGroupSize;

template<> struct QualityDecompressionGroupSize<1>{static constexpr int value = 32;};
template<> struct QualityDecompressionGroupSize<2>{static constexpr int value = 32;};
template<> struct QualityDecompressionGroupSize<8>{static constexpr int value = 32;};

template<int numBits, class LengthIterator>
__global__
void decompressQualityScoresKernel(
    char* __restrict__ quality, 
    std::size_t qualityPitchInBytes,
    const unsigned int* __restrict__ compressed, 
    std::size_t compressedPitchInInts,
    LengthIterator lengths,
    int numSequences
){
    constexpr int groupsize = QualityDecompressionGroupSize<numBits>::value;
    
    auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
    const int numGroups = (gridDim.x * blockDim.x) / groupsize;
    const int groupId = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;

    for(int s = groupId; s < numSequences; s += numGroups){
        const unsigned int* const myCompressed = compressed + compressedPitchInInts * s;
        char* const myQuality = quality + qualityPitchInBytes * s;
        const int myLength = lengths[s];

        QualityCompressor<numBits>::decodeQualityToString(group, myQuality, myCompressed, myLength);
    }
}

template<int numBits, class LengthIterator>
void callDecompressQualityScoresKernel(
    char* __restrict__ quality, 
    std::size_t qualityPitchInBytes,
    const unsigned int* __restrict__ compressed, 
    std::size_t compressedPitchInInts,
    LengthIterator lengths,
    int numSequences,
    cudaStream_t stream
){
    constexpr int groupsize = QualityDecompressionGroupSize<numBits>::value;
    constexpr int blocksize = 128;
    constexpr int groupsPerBlock = blocksize / groupsize;
    const std::size_t smem = 0;

    int deviceId = 0;
    int numSMs = 0;
    int maxBlocksPerSM = 0;
    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        decompressQualityScoresKernel<numBits, LengthIterator>,
        blocksize, 
        smem
    ));

    const int maxBlocks = maxBlocksPerSM * numSMs;

    const int numBlocks = std::min(maxBlocks, SDIV(numSequences, groupsPerBlock));

    decompressQualityScoresKernel<numBits><<<numBlocks, blocksize, 0, stream>>>(
        quality,
        qualityPitchInBytes,
        compressed,
        compressedPitchInInts,
        lengths,
        numSequences
    );
    CUDACHECKASYNC;
}

template<class LengthIterator>
void callDecompressQualityScoresKernel(
    char* __restrict__ quality, 
    std::size_t qualityPitchInBytes,
    const unsigned int* __restrict__ compressed, 
    std::size_t compressedPitchInInts,
    LengthIterator lengths,
    int numSequences,
    int numBits,
    cudaStream_t stream
){
    switch(numBits){
        case 1: callDecompressQualityScoresKernel<1>(quality, qualityPitchInBytes, compressed, compressedPitchInInts, lengths, numSequences, stream); break;
        case 2: callDecompressQualityScoresKernel<2>(quality, qualityPitchInBytes, compressed, compressedPitchInInts, lengths, numSequences, stream); break;
        case 8: callDecompressQualityScoresKernel<8>(quality, qualityPitchInBytes, compressed, compressedPitchInInts, lengths, numSequences, stream); break;
        default:
            std::cerr << "callDecompressQualityScoresKernel cannot be called with numBits = " << numBits << "\n";
            assert(false);
            break;
    }
}




#endif







#endif

