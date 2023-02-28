#ifndef CARE_SEQUENCEHELPERS_HPP
#define CARE_SEQUENCEHELPERS_HPP

#include <config.hpp>

#include "hpc_helpers.cuh"

#include <cstdint>
#include <string>
#include <type_traits>



    template<class T>
    struct EncodedReverseComplement2Bit{
    public:
        template<class I>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr T compute(I uinteger) noexcept{
            static_assert(std::is_same<T, I>::value, "types do no match");

            return compute_(uinteger);
        }
    private:
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr unsigned int compute_(unsigned int n) noexcept{
            n = ((n >> 2)  & 0x33333333u) | ((n & 0x33333333u) << 2);
            n = ((n >> 4)  & 0x0F0F0F0Fu) | ((n & 0x0F0F0F0Fu) << 4);
            n = ((n >> 8)  & 0x00FF00FFu) | ((n & 0x00FF00FFu) << 8);
            n = ((n >> 16) & 0x0000FFFFu) | ((n & 0x0000FFFFu) << 16);
            return ((unsigned int)(-1) - n) >> (8 * sizeof(n) - (16 << 1));
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint64_t compute_(std::uint64_t n) noexcept{
            n = ((n >> 2)  & 0x3333333333333333ull) | ((n & 0x3333333333333333ull) << 2);
            n = ((n >> 4)  & 0x0F0F0F0F0F0F0F0Full) | ((n & 0x0F0F0F0F0F0F0F0Full) << 4);
            n = ((n >> 8)  & 0x00FF00FF00FF00FFull) | ((n & 0x00FF00FF00FF00FFull) << 8);
            n = ((n >> 16) & 0x0000FFFF0000FFFFull) | ((n & 0x0000FFFF0000FFFFull) << 16);
            n = ((n >> 32) & 0x00000000FFFFFFFFull) | ((n & 0x00000000FFFFFFFFull) << 32);
            return ((std::uint64_t)(-1) - n) >> (8 * sizeof(n) - (32 << 1));
        }        
    };

    struct SequenceHelpers{
    public:

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr char complementBaseDecoded(char in) noexcept{
            switch(in){
                case 'A': return 'T';
                case 'C': return 'G';
                case 'G': return 'C';
                case 'T': return 'A';
            }
            return in;
        };

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void reverseComplementSequenceDecoded(char* reverseComplement, const char* sequence, int sequencelength) noexcept{
            for(int i = 0; i < sequencelength; ++i){
                switch(sequence[i]){
                    case 'A': reverseComplement[sequencelength-1-i] = 'T'; break;
                    case 'C': reverseComplement[sequencelength-1-i] = 'G'; break;
                    case 'G': reverseComplement[sequencelength-1-i] = 'C'; break;
                    case 'T': reverseComplement[sequencelength-1-i] = 'A'; break;
                    default : break; // don't change N
                }
            }
        }

        INLINEQUALIFIER
        static std::string reverseComplementSequenceDecoded(const char* sequence, int sequencelength){
            std::string rev;
            rev.resize(sequencelength);

            reverseComplementSequenceDecoded(&rev[0], sequence, sequencelength);

            return rev;
        }

        HD_WARNING_DISABLE
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void reverseComplementSequenceDecodedInplace(char* sequence, int sequencelength) noexcept{

            for(int i = 0; i < sequencelength/2; i++){
                const char front = complementBaseDecoded(sequence[i]);
                const char back = complementBaseDecoded(sequence[sequencelength - 1 - i]);
                sequence[i] = back;
                sequence[sequencelength - 1 - i] = front;
            }

            if(sequencelength % 2 == 1){
                const int middleindex = sequencelength/2;
                sequence[middleindex] = complementBaseDecoded(sequence[middleindex]);
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t encodedbaseA() noexcept{
            return 0;
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t encodedbaseC() noexcept{
            return 1;
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t encodedbaseG() noexcept{
            return 2;
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t encodedbaseT() noexcept{
            return 3;
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr unsigned int basesPerInt2Bit() noexcept{
            return sizeof(unsigned int) * CHAR_BIT / 2;
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr char decodeBase(std::uint8_t enc) noexcept{
            switch(enc){
            case encodedbaseA(): return 'A';
            case encodedbaseC(): return 'C';
            case encodedbaseG(): return 'G';
            case encodedbaseT(): return 'T';
            default: return '_';
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t encodeBase(char c) noexcept{
            switch(c){
            case 'A': return encodedbaseA();
            case 'C': return encodedbaseC();
            case 'G': return encodedbaseG();
            case 'T': return encodedbaseT();
            default: return encodedbaseA();
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr int getEncodedNumInts2Bit(int sequenceLength) noexcept{
            return SDIV(sequenceLength, basesPerInt2Bit());
        }

        HD_WARNING_DISABLE
        template<class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void encodeSequence2Bit(unsigned int* out, const char* sequence, int sequenceLength, IndexTransformation indextrafo) noexcept{

            const int nInts = getEncodedNumInts2Bit(sequenceLength);

            for(int i = 0; i < nInts; i++){
                out[indextrafo(i)] = 0;
            }

            for(int nucIndex = 0; nucIndex < sequenceLength; nucIndex++){
                const int intIndex = nucIndex / basesPerInt2Bit();
                switch(sequence[nucIndex]) {
                case 'A':
                    out[indextrafo(intIndex)] = (out[indextrafo(intIndex)] << 2) | encodedbaseA();
                    break;
                case 'C':
                    out[indextrafo(intIndex)] = (out[indextrafo(intIndex)] << 2) | encodedbaseC();
                    break;
                case 'G':
                    out[indextrafo(intIndex)] = (out[indextrafo(intIndex)] << 2) | encodedbaseG();
                    break;
                case 'T':
                    out[indextrafo(intIndex)] = (out[indextrafo(intIndex)] << 2) | encodedbaseT();
                    break;
                default:
                    out[indextrafo(intIndex)] = (out[indextrafo(intIndex)] << 2) | encodedbaseA();
                    break;
                }
            }
            //pack bits of last integer into higher order bits
            int leftoverbits = 2 * (nInts * basesPerInt2Bit() - sequenceLength);
            if(leftoverbits > 0){
                out[indextrafo(nInts-1)] <<= leftoverbits;
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void encodeSequence2Bit(unsigned int* outencoded, const char* sequence, int sequenceLength) noexcept{
            auto identity = [](auto i){return i;};
            encodeSequence2Bit(outencoded, sequence, sequenceLength, identity);
        }

        HD_WARNING_DISABLE
        template<class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t getEncodedNuc2Bit(const unsigned int* data, int /*sequenceLength*/, int i, IndexTransformation indextrafo) noexcept{
            const int intIndex = i / basesPerInt2Bit();
            const int pos = i % basesPerInt2Bit();
            return ((data[indextrafo(intIndex)] >> (30 - 2*pos)) & 3u);
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static std::uint8_t getEncodedNuc2Bit(const unsigned int* encodedsequence,
                                    int length,
                                    int position) noexcept{
            auto identity = [](auto i){return i;};
            return getEncodedNuc2Bit(encodedsequence, length, position, identity);
        }

        HD_WARNING_DISABLE
        template<class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void decode2BitSequence(char* sequence, const unsigned int* encoded, int sequenceLength, IndexTransformation indextrafo) noexcept{
            for(int i = 0; i < sequenceLength; i++){
                const std::uint8_t base = getEncodedNuc2Bit(encoded, sequenceLength, i, indextrafo);
                sequence[i] = decodeBase(base);
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void decode2BitSequence(char* sequence,
                                    const unsigned int* encodedsequence,
                                    int length) noexcept{
            auto identity = [](auto i){return i;};
            decode2BitSequence(sequence, encodedsequence, length, identity);
        }

        template<class IndexTransformation>
        static std::string get2BitString(const unsigned int* encoded, int sequenceLength, IndexTransformation indextrafo){
            std::string s;
            s.resize(sequenceLength);
            decode2BitSequence(&s[0], encoded, sequenceLength, indextrafo);
            return s;
        }

        static std::string get2BitString(const unsigned int* encodedsequence,
                                    int length){
            auto identity = [](auto i){return i;};
            return get2BitString(encodedsequence, length, identity);
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t getEncodedNucFromInt2Bit(unsigned int data, int pos){
            return ((data >> (30 - 2*pos)) & 3u);
        };

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr unsigned int reverseComplementInt2Bit(unsigned int n) noexcept{
            return EncodedReverseComplement2Bit<unsigned int>::compute(n);
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint64_t reverseComplementInt2Bit(std::uint64_t n) noexcept{
            return EncodedReverseComplement2Bit<std::uint64_t>::compute(n);
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t complementBase2Bit(std::uint8_t n) noexcept{
            return (~n & std::uint8_t{3});
        }

        HD_WARNING_DISABLE
        template<class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void reverseComplementSequenceInplace2Bit(unsigned int* encodedsequence, int sequenceLength, IndexTransformation indextrafo) noexcept{

            const int nInts = getEncodedNumInts2Bit(sequenceLength);
            const int unusedPositions = nInts * basesPerInt2Bit() - sequenceLength;

            for(int i = 0; i < nInts/2; i++){
                const unsigned int front = reverseComplementInt2Bit(encodedsequence[indextrafo(i)]);
                const unsigned int back = reverseComplementInt2Bit(encodedsequence[indextrafo(nInts - 1 - i)]);
                encodedsequence[indextrafo(i)] = back;
                encodedsequence[indextrafo(nInts - 1 - i)] = front;
            }

            if(nInts % 2 == 1){
                const int middleindex = nInts/2;
                encodedsequence[indextrafo(middleindex)] = reverseComplementInt2Bit(encodedsequence[indextrafo(middleindex)]);
            }

            if(unusedPositions > 0){
                for(int i = 0; i < nInts-1; i++){
                    encodedsequence[indextrafo(i)] = (encodedsequence[indextrafo(i)] << (2*unusedPositions))
                                                | (encodedsequence[indextrafo(i+1)] >> (2 * (basesPerInt2Bit()-unusedPositions)));

                }
                encodedsequence[indextrafo(nInts-1)] <<= (2*unusedPositions);
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void reverseComplementSequenceInplace2Bit(unsigned int* encodedsequence,
                                    int length) noexcept{
            auto identity = [](auto i){return i;};
            reverseComplementSequenceInplace2Bit(encodedsequence, length, identity);
        }

        HD_WARNING_DISABLE
        template<class RcIndexTransformation, class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void reverseComplementSequence2Bit(unsigned int* rcencodedsequence,
                                        const unsigned int* encodedsequence,
                                        int sequenceLength,
                                        RcIndexTransformation rcindextrafo,
                                        IndexTransformation indextrafo) noexcept{

            const int nInts = getEncodedNumInts2Bit(sequenceLength);
            for(int i = 0; i < nInts; i++){
                rcencodedsequence[rcindextrafo(i)] = encodedsequence[indextrafo(i)];
            }

            reverseComplementSequenceInplace2Bit(rcencodedsequence, sequenceLength, rcindextrafo);
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void reverseComplementSequence2Bit(unsigned int* rcencodedsequence,
                                    const unsigned int* encodedsequence,
                                    int length) noexcept{
            auto identity = [](auto i){return i;};
            reverseComplementSequence2Bit(rcencodedsequence, encodedsequence, length, identity, identity);
        }

        template<class CopyType = char, class Group>
        DEVICEQUALIFIER
        static constexpr void decodeSequence2Bit(Group& group, const unsigned int* input, int sequencelength, char* output){
            #ifdef __CUDA_ARCH__

            //output must be aligned to sizeof(CopyType) bytes
            constexpr int copyTypeSize = sizeof(CopyType);
            static_assert(copyTypeSize == 1 || copyTypeSize == 2 || copyTypeSize == 4 || copyTypeSize == 8 || copyTypeSize == 16, "Invalid CopyType");

            const int nInts = SequenceHelpers::getEncodedNumInts2Bit(sequencelength);
            constexpr int basesPerInt = SequenceHelpers::basesPerInt2Bit();

            for(int i = group.thread_rank(); i < nInts; i += group.size()){
                unsigned int data = input[i];

                if(i < nInts-1){
                    //not last iteration. int encodes 16 chars
                    __align__(16) char nucs[16];

                    #ifdef __CUDA_ARCH__
                    #pragma unroll
                    #endif
                    for(int p = 0; p < 16; p++){
                        const std::uint8_t encodedBase = SequenceHelpers::getEncodedNucFromInt2Bit(data, p);
                        nucs[p] = SequenceHelpers::decodeBase(encodedBase);
                    }

                    #ifdef __CUDA_ARCH__
                    #pragma unroll
                    #endif
                    for(int p = 0; p < 16 / copyTypeSize; p++){
                        ((CopyType*)output)[(16 / copyTypeSize)*i + p] = *((const CopyType*)&nucs[p]);
                    }
                }else{
                    const int remaining = sequencelength - i * basesPerInt;

                    for(int p = 0; p < remaining; p++){
                        const std::uint8_t encodedBase = SequenceHelpers::getEncodedNucFromInt2Bit(data, p);
                        output[i * basesPerInt + p] = SequenceHelpers::decodeBase(encodedBase);
                    }
                }
            }

            #endif
        }

        template<class Group>
        DEVICEQUALIFIER
        static void encodeSingleSequenceTo2Bit(
            Group& group,
            unsigned int* out,
            const char* in,
            int length
        ){
            #ifdef __CUDA_ARCH__
            if(length > 0){
                const int nInts = SequenceHelpers::getEncodedNumInts2Bit(length);
                constexpr int basesPerInt = SequenceHelpers::basesPerInt2Bit();

                for(int i = group.thread_rank(); i < nInts; i += group.size()){
                    unsigned int data = 0;

                    auto encodeNuc = [&](char nuc){
                        switch(nuc) {
                        case 'A':
                            data = (data << 2) | SequenceHelpers::encodedbaseA();
                            break;
                        case 'C':
                            data = (data << 2) | SequenceHelpers::encodedbaseC();
                            break;
                        case 'G':
                            data = (data << 2) | SequenceHelpers::encodedbaseG();
                            break;
                        case 'T':
                            data = (data << 2) | SequenceHelpers::encodedbaseT();
                            break;
                        default:
                            data = (data << 2) | SequenceHelpers::encodedbaseA();
                            break;
                        }
                    };

                    if(i < nInts - 1){
                        //not last iteration. int encodes 16 chars
                        __align__(16) char nucs[16];
                        ((int4*)nucs)[0] = *((const int4*)&in[i * 16]);

                        #pragma unroll
                        for(int p = 0; p < 16; p++){
                            encodeNuc(nucs[p]);
                        }
                    }else{        
                        for(int nucIndex = i * basesPerInt; nucIndex < length; nucIndex++){
                            encodeNuc(in[nucIndex]);
                        }

                        //pack bits of last integer into higher order bits
                        int leftoverbits = 2 * (nInts * basesPerInt - length);
                        if(leftoverbits > 0){
                            data <<= leftoverbits;
                        }

                    }

                    out[i] = data;
                }
            }
            #endif
        }





        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr int getEncodedNumInts2BitHiLo(int sequenceLength) noexcept{
            return int(2 * SDIV(sequenceLength, sizeof(unsigned int) * CHAR_BIT));
        }

        HD_WARNING_DISABLE
        template<class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void encodeSequence2BitHiLo(unsigned int* out, const char* sequence, int sequenceLength, IndexTransformation indextrafo) noexcept{
            const int nInts = getEncodedNumInts2BitHiLo(sequenceLength);

            for(int i = 0; i < nInts; i++){
                out[indextrafo(i)] = 0;
            }

            unsigned int* const hi = out;
            unsigned int* const lo = out + indextrafo(nInts/2);

            for(int i = 0; i < sequenceLength; i++){
                const int intIndex = i / (CHAR_BIT * sizeof(unsigned int));

                switch(sequence[i]) {
                case 'A':
                    hi[indextrafo(intIndex)] = (hi[indextrafo(intIndex)] << 1) | 0;
                    lo[indextrafo(intIndex)] = (lo[indextrafo(intIndex)] << 1) | 0;
                    break;
                case 'C':
                    hi[indextrafo(intIndex)] = (hi[indextrafo(intIndex)] << 1) | 0;
                    lo[indextrafo(intIndex)] = (lo[indextrafo(intIndex)] << 1) | 1;
                    break;
                case 'G':
                    hi[indextrafo(intIndex)] = (hi[indextrafo(intIndex)] << 1) | 1;
                    lo[indextrafo(intIndex)] = (lo[indextrafo(intIndex)] << 1) | 0;
                    break;
                case 'T':
                    hi[indextrafo(intIndex)] = (hi[indextrafo(intIndex)] << 1) | 1;
                    lo[indextrafo(intIndex)] = (lo[indextrafo(intIndex)] << 1) | 1;
                    break;
                default:
                    hi[indextrafo(intIndex)] = (hi[indextrafo(intIndex)] << 1) | 0;
                    lo[indextrafo(intIndex)] = (lo[indextrafo(intIndex)] << 1) | 0;
                    break;
                }
            }
            //pack bits of last hi integer and lo integer into their higher order bits
            const int leftoverbits = nInts/2 * CHAR_BIT * sizeof(unsigned int) - sequenceLength;
            if(leftoverbits > 0){
                hi[indextrafo(nInts/2-1)] <<= leftoverbits;
                lo[indextrafo(nInts/2-1)] <<= leftoverbits;
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void encodeSequence2BitHiLo(unsigned int* outencoded, const char* sequence, int sequenceLength) noexcept{
            auto identity = [](auto i){return i;};
            encodeSequence2BitHiLo(outencoded, sequence, sequenceLength, identity);
        }

        HD_WARNING_DISABLE
        template<class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr std::uint8_t getEncodedNuc2BitHiLo(const unsigned int* data, int sequenceLength, int i, IndexTransformation indextrafo) noexcept{
            const int nInts = getEncodedNumInts2BitHiLo(sequenceLength);

            const unsigned int* const hi = data;
            const unsigned int* const lo = data + indextrafo(nInts/2);

            const int intIndex = i / (CHAR_BIT * sizeof(unsigned int));
            const int pos = i % (CHAR_BIT * sizeof(unsigned int));
            const unsigned int hibit = (hi[indextrafo(intIndex)] >> (31 - pos)) & 1u;
            const unsigned int lobit = (lo[indextrafo(intIndex)] >> (31 - pos)) & 1u;
            return (hibit << 1) | lobit;
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static std::uint8_t getEncodedNuc2BitHiLo(const unsigned int* encodedsequence,
                                    int length,
                                    int position) noexcept{
            auto identity = [](auto i){return i;};
            return getEncodedNuc2BitHiLo(encodedsequence, length, position, identity);
        }

        HD_WARNING_DISABLE
        template<class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void decode2BitHiLoSequence(char* sequence, const unsigned int* encoded, int sequenceLength, IndexTransformation indextrafo) noexcept{
            for(int i = 0; i < sequenceLength; i++){
                const std::uint8_t base = getEncodedNuc2BitHiLo(encoded, sequenceLength, i, indextrafo);

                switch(base){
                case encodedbaseA(): sequence[i] = 'A'; break;
                case encodedbaseC(): sequence[i] = 'C'; break;
                case encodedbaseG(): sequence[i] = 'G'; break;
                case encodedbaseT(): sequence[i] = 'T'; break;
                //default: sequence[i] = '_'; break; // cannot happen
                }
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void decode2BitHiLoSequence(char* sequence,
                                    const unsigned int* encodedsequence,
                                    int length) noexcept{
            auto identity = [](auto i){return i;};
            decode2BitHiLoSequence(sequence, encodedsequence, length, identity);
        }

        template<class IndexTransformation>
        INLINEQUALIFIER
        static std::string get2BitHiLoString(const unsigned int* encoded, int sequenceLength, IndexTransformation indextrafo){
            std::string s;
            s.resize(sequenceLength);
            decode2BitHiLoSequence(&s[0], encoded, sequenceLength, indextrafo);
            return s;
        }

        INLINEQUALIFIER
        static std::string get2BitHiLoString(const unsigned int* encodedsequence,
                                    int length){
            auto identity = [](auto i){return i;};
            return get2BitHiLoString(encodedsequence, length, identity);
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr unsigned int reverseComplementInt2BitHiLoHalf(unsigned int n) noexcept{
            n = ((n >> 1) & 0x55555555) | ((n << 1) & 0xaaaaaaaa);
            n = ((n >> 2) & 0x33333333) | ((n << 2) & 0xcccccccc);
            n = ((n >> 4) & 0x0f0f0f0f) | ((n << 4) & 0xf0f0f0f0);
            n = ((n >> 8) & 0x00ff00ff) | ((n << 8) & 0xff00ff00);
            n = ((n >> 16) & 0x0000ffff) | ((n << 16) & 0xffff0000);
            return ~n;
        };


        HD_WARNING_DISABLE
        template<class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void reverseComplementSequenceInplace2BitHiLo(unsigned int* encodedsequence, int sequenceLength, IndexTransformation indextrafo) noexcept{
            const int ints = getEncodedNumInts2BitHiLo(sequenceLength);
            const int unusedBitsInt = SDIV(sequenceLength, CHAR_BIT * sizeof(unsigned int)) * CHAR_BIT * sizeof(unsigned int) - sequenceLength;

            unsigned int* const hi = encodedsequence;
            unsigned int* const lo = hi + indextrafo(ints/2);

            const int intsPerHalf = SDIV(sequenceLength, CHAR_BIT * sizeof(unsigned int));
            for(int i = 0; i < intsPerHalf/2; ++i){
                const unsigned int hifront = reverseComplementInt2BitHiLoHalf(hi[indextrafo(i)]);
                const unsigned int hiback = reverseComplementInt2BitHiLoHalf(hi[indextrafo(intsPerHalf - 1 - i)]);
                hi[indextrafo(i)] = hiback;
                hi[indextrafo(intsPerHalf - 1 - i)] = hifront;

                const unsigned int lofront = reverseComplementInt2BitHiLoHalf(lo[indextrafo(i)]);
                const unsigned int loback = reverseComplementInt2BitHiLoHalf(lo[indextrafo(intsPerHalf - 1 - i)]);
                lo[indextrafo(i)] = loback;
                lo[indextrafo(intsPerHalf - 1 - i)] = lofront;
            }
            if(intsPerHalf % 2 == 1){
                const int middleindex = intsPerHalf/2;
                hi[indextrafo(middleindex)] = reverseComplementInt2BitHiLoHalf(hi[indextrafo(middleindex)]);
                lo[indextrafo(middleindex)] = reverseComplementInt2BitHiLoHalf(lo[indextrafo(middleindex)]);
            }

            if(unusedBitsInt != 0){
                for(int i = 0; i < intsPerHalf - 1; ++i){
                    hi[indextrafo(i)] = (hi[indextrafo(i)] << unusedBitsInt) | (hi[indextrafo(i+1)] >> (CHAR_BIT * sizeof(unsigned int) - unusedBitsInt));
                    lo[indextrafo(i)] = (lo[indextrafo(i)] << unusedBitsInt) | (lo[indextrafo(i+1)] >> (CHAR_BIT * sizeof(unsigned int) - unusedBitsInt));
                }

                hi[indextrafo(intsPerHalf - 1)] <<= unusedBitsInt;
                lo[indextrafo(intsPerHalf - 1)] <<= unusedBitsInt;
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void reverseComplementSequenceInplace2BitHiLo(unsigned int* encodedsequence,
                                    int length) noexcept{
            auto identity = [](auto i){return i;};
            reverseComplementSequenceInplace2BitHiLo(encodedsequence, length, identity);
        }

        HD_WARNING_DISABLE
        template<class RcIndexTransformation, class IndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void reverseComplementSequence2BitHiLo(unsigned int* rcencodedsequence,
                                        const unsigned int* encodedsequence,
                                        int sequenceLength,
                                        RcIndexTransformation rcindextrafo,
                                        IndexTransformation indextrafo) noexcept{

            const int nInts = getEncodedNumInts2BitHiLo(sequenceLength);
            for(int i = 0; i < nInts; i++){
                rcencodedsequence[rcindextrafo(i)] = encodedsequence[indextrafo(i)];
            }

            reverseComplementSequenceInplace2BitHiLo(rcencodedsequence, sequenceLength, rcindextrafo);
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void reverseComplementSequence2BitHiLo(unsigned int* rcencodedsequence,
                                    const unsigned int* encodedsequence,
                                    int length) noexcept{
            auto identity = [](auto i){return i;};
            reverseComplementSequence2BitHiLo(rcencodedsequence, encodedsequence, length, identity, identity);
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr unsigned int extractEvenBits(unsigned int x) noexcept{
            x = x & 0x55555555;
            x = (x | (x >> 1)) & 0x33333333;
            x = (x | (x >> 2)) & 0x0F0F0F0F;
            x = (x | (x >> 4)) & 0x00FF00FF;
            x = (x | (x >> 8)) & 0x0000FFFF;
            return x;
        }
        
        HD_WARNING_DISABLE
        template<class InIndexTransformation, class OutIndexTransformation>
        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static constexpr void convert2BitTo2BitHiLo(unsigned int* out,
                                    const unsigned int* in,
                                    int length,
                                    InIndexTransformation inindextrafo,
                                    OutIndexTransformation outindextrafo) noexcept{

            const int inInts = getEncodedNumInts2Bit(length);
            const int outInts = getEncodedNumInts2BitHiLo(length);

            unsigned int* const outHi = out;
            unsigned int* const outLo = out + outindextrafo(outInts/2);

            for(int i = 0; i < outInts; i++){
                out[outindextrafo(i)] = 0;
            }

            for(int i = 0; i < inInts; i++){
                const int inindex = inindextrafo(i);
                const int outindex = outindextrafo(i/2);

                unsigned int even16 = extractEvenBits(in[inindex]);
                unsigned int odd16 = extractEvenBits(in[inindex] >> 1);

                if(i % 2 == 0){
                    outHi[outindex] = odd16 << 16;
                    outLo[outindex] = even16 << 16;
                }else{
                    outHi[outindex] = outHi[outindex] | odd16;
                    outLo[outindex] = outLo[outindex] | even16;
                }
            }
        }

        HOSTDEVICEQUALIFIER INLINEQUALIFIER
        static void convert2BitTo2BitHiLo(unsigned int* out,
                                    const unsigned int* in,
                                    int length) noexcept{
            auto identity = [](auto i){return i;};
            convert2BitTo2BitHiLo(out, in, length, identity, identity);
        }

        //get kmer of length k which starts at position pos of the decoded sequence
        HOSTDEVICEQUALIFIER
        static constexpr std::uint64_t getEncodedKmerFromEncodedSequence(const unsigned int* encodedSequence, int k, int pos){
            assert(k > 0);

            constexpr int maximum_kmer_length = care::max_k<std::uint64_t>::value;
            const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);

            assert(k <= maximum_kmer_length);

            std::uint64_t kmer = 0;

            while(k > 0){
                const int intIndex = pos / 16;
                const std::uint64_t intData = encodedSequence[intIndex];
                const int neededBasesFromInt = std::min(k, 16 - (pos % 16));
                assert(neededBasesFromInt > 0);
                const int intShift = 2*(16 - neededBasesFromInt - (pos % 16));
                kmer = (kmer << (2*neededBasesFromInt)) | (intData >> intShift);
                k -= neededBasesFromInt;
                pos += neededBasesFromInt;
            }           

            return kmer & kmer_mask;
        }

        template<class Func> // Func::operator()(std::uint64_t kmer, int pos)
        HOSTDEVICEQUALIFIER
        static constexpr void forEachEncodedKmerFromEncodedSequence(const unsigned int* encodedSequence, int sequenceLength, int k, int first, int last /*exclusive*/, Func callback){
            if(sequenceLength <= 0 || k > sequenceLength) return;
            last = std::min(last, sequenceLength - k + 1);
            if(first >= last) return;

            assert(k > 0);

            constexpr int maximum_kmer_length = care::max_k<std::uint64_t>::value;
            const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);

            assert(k <= maximum_kmer_length);

            std::uint64_t kmer = getEncodedKmerFromEncodedSequence(encodedSequence, k, first);
            callback(kmer, first);

            for(int pos = first + 1; pos < last; pos++){
                const int nextIntIndex = (pos + k - 1) / basesPerInt2Bit();
                const int nextPositionInInt = (pos + k - 1) % basesPerInt2Bit();

                const std::uint64_t nextBase = encodedSequence[nextIntIndex] >> (30 - 2 * nextPositionInInt);

                kmer = ((kmer << 2) | nextBase) & kmer_mask;
                //assert(kmer == getEncodedKmerFromEncodedSequence(encodedSequence, k, pos));
                callback(kmer, pos);
            }
        }

        template<class Func> // Func::operator()(std::uint64_t kmer, int pos)
        HOSTDEVICEQUALIFIER
        static constexpr void forEachEncodedKmerFromEncodedSequence(const unsigned int* encodedSequence, int sequenceLength, int k, Func callback){
            forEachEncodedKmerFromEncodedSequence(encodedSequence, sequenceLength, k,  0, sequenceLength - k + 1, callback);
        }

        template<class Func> // Func::operator()(std::uint64_t kmer, int pos)
        HOSTDEVICEQUALIFIER
        static constexpr void forEachEncodedKmerFromEncodedSequence_elaborate(const unsigned int* encodedSequence, int sequenceLength, int k, int first, int last /*exclusive*/, Func callback){
            if(sequenceLength <= 0 || k > sequenceLength) return;
            last = std::min(last, sequenceLength - k + 1);
            if(first >= last) return;

            const int newSequenceLength = last + k - 1;

            assert(k > 0);

            constexpr int maximum_kmer_length = care::max_k<std::uint64_t>::value;
            const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);

            assert(k <= maximum_kmer_length);

            std::uint64_t kmer = getEncodedKmerFromEncodedSequence(encodedSequence, k, first);
            kmer >>= 2; //k-1 bases, allows easier loop

            constexpr int basesPerInt = SequenceHelpers::basesPerInt2Bit();

            // first == 1, last == 84, k == 16, length 100
            //itersend1 = min(16, 100) = 16
            // nextSequencePos = 16 -> 0 kmers
            //fullIntIters = 5 -> 80 kmers

            const int itersend1 = std::min(SDIV(first + k-1, basesPerInt) * basesPerInt, newSequenceLength);
            int pos = first;

            //process sequence positions one by one
            // until the next encoded sequence data element is reached
            for(int nextSequencePos = first + k - 1; nextSequencePos < itersend1; nextSequencePos++){
                const int nextIntIndex = nextSequencePos / basesPerInt;
                const int nextPositionInInt = nextSequencePos % basesPerInt;

                const std::uint64_t nextBase = encodedSequence[nextIntIndex] >> (30 - 2 * nextPositionInInt);
                kmer = ((kmer << 2) | nextBase) & kmer_mask;

                callback(kmer, pos);
                pos++;
            }

            const int fullIntIters = (newSequenceLength - itersend1) / basesPerInt;

            //process all fully occupied encoded sequence data elements
            // improves memory access
            for(int iter = 0; iter < fullIntIters; iter++){
                const int intIndex = (itersend1 + iter * basesPerInt) / basesPerInt;
                const unsigned int data = encodedSequence[intIndex];

                // #ifdef __CUDA_ARCH__
                // #pragma unroll
                // #endif
                for(int posInInt = 0; posInInt < basesPerInt; posInInt++){
                    const std::uint64_t nextBase = data >> (30 - 2 * posInInt);
                    kmer = ((kmer << 2) | nextBase) & kmer_mask;

                    callback(kmer, pos);
                    pos++;
                }
            }

            //process remaining positions one by one
            for(int nextSequencePos = fullIntIters * basesPerInt + itersend1; nextSequencePos < newSequenceLength; nextSequencePos++){
                const int nextIntIndex = nextSequencePos / basesPerInt;
                const int nextPositionInInt = nextSequencePos % basesPerInt;

                const std::uint64_t nextBase = encodedSequence[nextIntIndex] >> (30 - 2 * nextPositionInInt);
                kmer = ((kmer << 2) | nextBase) & kmer_mask;

                callback(kmer, pos);
                pos++;
            }
        }

        template<class Func> // Func::operator()(std::uint64_t kmer, int pos)
        HOSTDEVICEQUALIFIER
        static constexpr void forEachEncodedKmerFromEncodedSequence_elaborate(const unsigned int* encodedSequence, int sequenceLength, int k, Func callback){
            forEachEncodedKmerFromEncodedSequence_elaborate(encodedSequence, sequenceLength, k,  0, sequenceLength - k + 1, callback);
        }


        template<class Func> // Func::operator()(std::uint64_t kmer, int pos)
        HOSTDEVICEQUALIFIER
        static constexpr void forEachEncodedCanonicalKmerFromEncodedSequence(const unsigned int* encodedSequence, int sequenceLength, int k, int first, int last /*exclusive*/, Func callback){
            if(sequenceLength <= 0 || k > sequenceLength) return;
            last = std::min(last, sequenceLength - k + 1);
            if(first >= last) return;

            const int newSequenceLength = last + k - 1;

            assert(k > 0);

            constexpr int maximum_kmer_length = care::max_k<std::uint64_t>::value;
            const std::uint64_t kmer_mask = std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - k) * 2);
            const int rcshiftamount = (maximum_kmer_length - k) * 2;

            assert(k <= maximum_kmer_length);

            std::uint64_t kmer = getEncodedKmerFromEncodedSequence(encodedSequence, k, first);
            kmer >>= 2; //k-1 bases, allows easier loop
            std::uint64_t rc_kmer = SequenceHelpers::reverseComplementInt2Bit(kmer);

            auto addBase = [&](std::uint64_t encBase){
                const std::uint64_t revcBase = (~encBase) & 3;
                kmer = ((kmer << 2) | encBase) & kmer_mask;
                rc_kmer = (rc_kmer >> 2) | (revcBase << 62);
            };


            constexpr int basesPerInt = SequenceHelpers::basesPerInt2Bit();

            const int itersend1 = std::min(SDIV(first + k-1, basesPerInt) * basesPerInt, newSequenceLength);
            int pos = first;

            //process sequence positions one by one
            // until the next encoded sequence data element is reached
            for(int nextSequencePos = first + k - 1; nextSequencePos < itersend1; nextSequencePos++){
                const int nextIntIndex = nextSequencePos / basesPerInt;
                const int nextPositionInInt = nextSequencePos % basesPerInt;

                const std::uint64_t nextBase = encodedSequence[nextIntIndex] >> (30 - 2 * nextPositionInInt);
                addBase(nextBase);                
                const std::uint64_t smallest = std::min(kmer, rc_kmer >> rcshiftamount);

                callback(smallest, pos);
                pos++;
            }

            const int fullIntIters = (newSequenceLength - itersend1) / basesPerInt;

            //process all fully occupied encoded sequence data elements
            // improves memory access
            for(int iter = 0; iter < fullIntIters; iter++){
                const int intIndex = (itersend1 + iter * basesPerInt) / basesPerInt;
                const unsigned int data = encodedSequence[intIndex];

                // #ifdef __CUDA_ARCH__
                // #pragma unroll
                // #endif
                for(int posInInt = 0; posInInt < basesPerInt; posInInt++){
                    const std::uint64_t nextBase = data >> (30 - 2 * posInInt);

                    addBase(nextBase);                
                    const std::uint64_t smallest = std::min(kmer, rc_kmer >> rcshiftamount);

                    callback(smallest, pos);
                    pos++;
                }
            }

            //process remaining positions one by one
            for(int nextSequencePos = fullIntIters * basesPerInt + itersend1; nextSequencePos < newSequenceLength; nextSequencePos++){
                const int nextIntIndex = nextSequencePos / basesPerInt;
                const int nextPositionInInt = nextSequencePos % basesPerInt;

                const std::uint64_t nextBase = encodedSequence[nextIntIndex] >> (30 - 2 * nextPositionInInt);
                addBase(nextBase);                
                const std::uint64_t smallest = std::min(kmer, rc_kmer >> rcshiftamount);

                callback(smallest, pos);
                pos++;
            }
        }


        template<class Func> // Func::operator()(std::uint64_t kmer, int pos)
        HOSTDEVICEQUALIFIER
        static constexpr void forEachEncodedCanonicalKmerFromEncodedSequence(const unsigned int* encodedSequence, int sequenceLength, int k, Func callback){
            forEachEncodedCanonicalKmerFromEncodedSequence(encodedSequence, sequenceLength, k,  0, sequenceLength - k + 1, callback);
        }

        template<class Func> // Func::operator()(std::uint64_t kmer, int windowId, int globalPos)
        HOSTDEVICEQUALIFIER
        static constexpr void forEachEncodedCanonicalKmerInWindowsFromEncodedSequence(const unsigned int* encodedSequence, int sequenceLength, int windowsize, int k, Func callback){
            if(k > windowsize) return;

            const int kmersInWindow = windowsize - k + 1;

            auto windowcallback = [&](std::uint64_t kmer, int pos){
                const int windowId = pos / kmersInWindow;

                callback(kmer, windowId, pos);
            };

            for(int windowBegin = 0; windowBegin < sequenceLength; windowBegin += kmersInWindow){
                forEachEncodedCanonicalKmerFromEncodedSequence(
                    encodedSequence, sequenceLength, k, windowBegin, windowBegin + kmersInWindow,
                    windowcallback
                );
            }
        }


        #ifdef __CUDACC__

        template<class Group>
        DEVICEQUALIFIER
        static void reverseAlignedDecodedSequenceWithGroupShfl(Group& group, char* sequence, int sequenceLength){

            auto reverse = [](char4 data){
                char4 s;
                s.x = data.w;
                s.y = data.z;
                s.z = data.y;
                s.w = data.x;
                return s;
            };
        
            auto shiftLeft1 = [](char4 data){
                char4 s;
                s.x = data.y;
                s.y = data.z;
                s.z = data.w;
                s.w = '\0';
                return s;
            };
        
            auto shiftLeft2 = [](char4 data){
                char4 s;
                s.x = data.z;
                s.y = data.w;
                s.z = '\0';
                s.w = '\0';
                return s;
            };
        
            auto shiftLeft3 = [](char4 data){
                char4 s;
                s.x = data.w;
                s.y = '\0';
                s.z = '\0';
                s.w = '\0';
                return s;
            };
        
            //treat [left,right] as "char8", shift to the left by one char. return leftmost 4 chars
            auto handleUnusedPositions1 = [](char4 left, char4 right){
                char4 s;
                s.x = left.y;
                s.y = left.z;
                s.z = left.w;
                s.w = right.x;
                return s;
            };
        
            //treat [left,right] as "char8", shift to the left by two chars. return leftmost 4 chars
            auto handleUnusedPositions2 = [](char4 left, char4 right){
                char4 s;
                s.x = left.z;
                s.y = left.w;
                s.z = right.x;
                s.w = right.y;
                return s;
            };
        
            //treat [left,right] as "char8", shift to the left by three chars. return leftmost 4 chars
            auto handleUnusedPositions3 = [](char4 left, char4 right){
                char4 s;
                s.x = left.w;
                s.y = right.x;
                s.z = right.y;
                s.w = right.z;
                return s;
            };
        
            if(sequenceLength <= 1) return;
        
            const int arrayLength = SDIV(sequenceLength, 4); // 4 bases per int
            const int unusedPositions = arrayLength * 4 - sequenceLength;
            char4* sequenceAsChar4 = (char4*)sequence;
        
            for(int i = group.thread_rank(); i < arrayLength/2; i += group.size()){
                const char4 fdata = ((char4*)sequence)[i];
                const char4 bdata = ((char4*)sequence)[arrayLength - 1 - i];
        
                const char4 front = reverse(fdata);
                const char4 back = reverse(bdata);
                sequenceAsChar4[i] = back;
                sequenceAsChar4[arrayLength - 1 - i] = front;
            }
        
            if(arrayLength % 2 == 1 && group.thread_rank() == 0){
                const int middleindex = arrayLength/2;
                const char4 mdata = ((char4*)sequence)[middleindex];
                sequenceAsChar4[middleindex] = reverse(mdata);
            }
        
            group.sync();
        
            if(unusedPositions > 0){
        
                char4 left;
                char4 right;
                char4 tmp;
        
                const int numIterations = SDIV(arrayLength-1, group.size());
        
                for(int iteration = 0; iteration < numIterations; iteration++){
                    const int index = iteration * group.size() + group.thread_rank();
                    if(index < arrayLength){
                        left = sequenceAsChar4[index];
                    }
                    const int index2 = (iteration+1) * group.size() + group.thread_rank();
                    if(index2 < arrayLength && group.thread_rank() == 0){
                        tmp = sequenceAsChar4[index2];
                    }
                    #if __CUDACC_VER_MAJOR__ < 11
                    //CUDA < 11 does not have shuffle api for char4
                    *((int*)(&right)) = group.shfl_down(*((const int*)(&left)), 1);
                    *((int*)(&tmp)) = group.shfl(*((const int*)(&tmp)), 0);
                    #else
                    right = group.shfl_down(left, 1);
                    tmp = group.shfl(tmp, 0);
                    #endif
                    if(group.thread_rank() == group.size() - 1){
                        right = tmp;
                    }
        
                    if(unusedPositions == 1){
                        char4 result = handleUnusedPositions1(left, right);
                        if(index < arrayLength - 1){
                            sequenceAsChar4[index] = result;
                        }
                    }else if(unusedPositions == 2){
                        char4 result = handleUnusedPositions2(left, right);
                        if(index < arrayLength - 1){
                            sequenceAsChar4[index] = result;
                        }
                    }else{
                        char4 result = handleUnusedPositions3(left, right);
                        if(index < arrayLength - 1){
                            sequenceAsChar4[index] = result;
                        }
                    }
                }
        
                group.sync();
        
                if(group.thread_rank() == 0){
                    if(unusedPositions == 1){
                        sequenceAsChar4[arrayLength-1] = shiftLeft1(sequenceAsChar4[arrayLength-1]);
                    }else if(unusedPositions == 2){
                        sequenceAsChar4[arrayLength-1] = shiftLeft2(sequenceAsChar4[arrayLength-1]);
                    }else{
                        assert(unusedPositions == 3);
                        sequenceAsChar4[arrayLength-1] = shiftLeft3(sequenceAsChar4[arrayLength-1]);
                    }
                }
            }
        }
    
        template<class Group>
        DEVICEQUALIFIER
        static void reverseComplementDecodedSequence(Group& group, char* sequence, int sequenceLength){
            const bool isFourByteAligned = (((unsigned long)sequence) % 4) == 0;
            if(isFourByteAligned){
                for(int i = group.thread_rank(); i < sequenceLength; i += group.size()) {
                    sequence[i] = SequenceHelpers::complementBaseDecoded(sequence[i]);
                }
                group.sync(); // threads may access elements which were written by another thread
                SequenceHelpers::reverseAlignedDecodedSequenceWithGroupShfl(group, sequence, sequenceLength);
                group.sync();
            }else{
                for(int i = group.thread_rank(); i < sequenceLength / 2; i += group.size()) {
                    const char l = SequenceHelpers::complementBaseDecoded(sequence[i]);
                    const char r = SequenceHelpers::complementBaseDecoded(sequence[sequenceLength - i - 1]);
                    sequence[i] = r;
                    sequence[sequenceLength - i - 1] = l;
                }
                if(group.thread_rank() == 0 && sequenceLength % 2 == 1){
                    sequence[sequenceLength / 2] = SequenceHelpers::complementBaseDecoded(sequence[sequenceLength / 2]);
                }
                group.sync();
            }
        }

        #endif
    };


#endif
