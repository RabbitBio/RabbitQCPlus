#ifndef CARE_BITCOMPRESSEDSTRING_HPP
#define CARE_BITCOMPRESSEDSTRING_HPP

#include <algorithm>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <string_view>
#include <cassert>
#include <cmath>

namespace care{


template<class Data_t = std::uint32_t, class Element_t = char>
struct BitCompressedString{
    //using Data_t = std::uint32_t;
    static_assert(std::is_unsigned<Data_t>::value == true, "");

    BitCompressedString() = default;

    BitCompressedString(const std::string_view sv){
        numElements = sv.size();
        auto minmaxiters = std::minmax_element(sv.begin(), sv.end());
        minElement = *minmaxiters.first;
        maxElement = *minmaxiters.second;
        
        assert(minElement <= maxElement);
        assert(numElements <= std::numeric_limits<int>::max());

        const int diff = maxElement - minElement;
        if(diff == 0){
            bitsPerElement = 0;
            return;
        }

        auto nextLargerPowerOfTwo = [](std::uint32_t n){
            n |= n >> 1;
            n |= n >> 2;
            n |= n >> 4;
            n |= n >> 8;
            n |= n >> 16;
            return n + 1;
        };

        auto nextP2 = nextLargerPowerOfTwo(diff);
        assert(nextP2 != std::numeric_limits<Data_t>::max());

        bitsPerElement = int(std::log2(nextP2));

        numData = ((bitsPerElement * size()) + dataTBits() - 1) / dataTBits();

        data = std::make_unique<Data_t[]>(numData);

        for(int i = 0; i < numElements; i++){
            setElement(i, sv[i]);
        }
    }

    BitCompressedString(const BitCompressedString&) = delete;
    BitCompressedString(BitCompressedString&&) = default;
    BitCompressedString& operator=(const BitCompressedString&) = delete;
    BitCompressedString& operator=(BitCompressedString&&) = default;


    void destroy(){
        data = nullptr;
        numElements = 0;
        minElement = 0;
        maxElement = 0;
        bitsPerElement = 0;
        numData = 0;
    }

    bool operator!=(const BitCompressedString& rhs) const noexcept{
        return !(operator==(rhs));
    }

    int getMinElement() const{
        return minElement;
    }

    int getMaxElement() const{
        return maxElement;
    }

    int size() const{
        return numElements;
    }

    int getRawbitsPerElement() const{
        return bitsPerElement;
    }

    Element_t getElement(int index) const{
        assert(index < size());

        if(getMinElement() == getMaxElement()){
            return getMinElement();
        }

        const std::uint64_t firstBit = bitsPerElement * index;
        const std::uint64_t lastBitExcl = bitsPerElement * (index+1);
        const std::uint64_t firstuintindex = firstBit / dataTBits();
        const int begin = firstBit - firstuintindex * dataTBits();
        const int endExcl = lastBitExcl - firstuintindex * dataTBits();
        
        const auto first = data[firstuintindex];
        //prevent oob access
        const auto second = int(firstuintindex) == numData - 1 ? data[firstuintindex] : data[firstuintindex + 1];
        const Data_t elementBits = getBits(first, second, begin, endExcl);

        return Element_t(elementBits) + minElement;
    }

    std::string getString() const{
        std::string result;
        result.reserve(size());

        for(int i = 0; i < size(); i++){
            result.push_back(getElement(i));
        }

        return result;
    }

    std::size_t getSerializedNumBytes() const{
        return sizeof(std::uint8_t) // bitsPerElement
            + sizeof(Element_t) // minElement
            + sizeof(Element_t) // maxElement
            + sizeof(int) // numElements
            + sizeof(int) // numData
            + sizeof(Data_t) * numData; // data
    }

    std::uint8_t* copyToContiguousMemory(std::uint8_t* ptr, std::uint8_t* endPtr) const{
        const std::size_t requiredBytes = getSerializedNumBytes();
        const std::size_t availableBytes = std::distance(ptr, endPtr);

        if(requiredBytes <= availableBytes){                
            std::memcpy(ptr, &bitsPerElement, sizeof(std::uint8_t)); ptr += sizeof(std::uint8_t);
            std::memcpy(ptr, &minElement, sizeof(Element_t)); ptr += sizeof(Element_t);
            std::memcpy(ptr, &maxElement, sizeof(Element_t)); ptr += sizeof(Element_t);
            std::memcpy(ptr, &numElements, sizeof(int)); ptr += sizeof(int);
            std::memcpy(ptr, &numData, sizeof(int)); ptr += sizeof(int);
            std::memcpy(ptr, data.get(), sizeof(Data_t) * numData); ptr += sizeof(Data_t) * numData;

            return ptr;
        }else{
            return nullptr;
        }   
    }

    const std::uint8_t* copyFromContiguousMemory(const std::uint8_t* ptr){
        std::memcpy(&bitsPerElement, ptr, sizeof(std::uint8_t)); ptr += sizeof(std::uint8_t);
        std::memcpy(&minElement, ptr, sizeof(Element_t)); ptr += sizeof(Element_t);
        std::memcpy(&maxElement, ptr, sizeof(Element_t)); ptr += sizeof(Element_t);
        std::memcpy(&numElements, ptr, sizeof(int)); ptr += sizeof(int);
        std::memcpy(&numData, ptr, sizeof(int)); ptr += sizeof(int);

        data = std::make_unique<Data_t[]>(numData);
        std::memcpy(data.get(), ptr, sizeof(Data_t) * numData); ptr += sizeof(Data_t) * numData;

        return ptr;
    }

private:   
    //extracts at most 32 bits out of the 64 bits [l,r]. extracts bits [begin,endExcl) from [l,r]
    //bits are enumerated from left to right
    Data_t getBits(Data_t l, Data_t r, std::size_t begin, std::size_t endExcl) const{

        assert(begin < endExcl && endExcl <= 2 * dataTBits());
        
        const int lbegin = std::min(dataTBits(), begin);
        const int lendExcl = std::min(dataTBits(), endExcl);

        const int rbegin = (begin >= dataTBits()) ? begin - dataTBits() : 0;
        const int rendExcl = (endExcl >= dataTBits()) ? endExcl - dataTBits() : 0;

        Data_t lmask = 0;
        for(int i = lbegin; i < lendExcl; i++){
            lmask = (lmask << 1) | 1;
        }

        Data_t rmask = 0;
        for(int i = rbegin; i < rendExcl; i++){
            rmask = (rmask << 1) | 1;
        }

        const Data_t lpiece = (l >> (dataTBits() - lendExcl)) & lmask;
        const Data_t rpiece = (r >> (dataTBits() - rendExcl)) & rmask;
        Data_t result = lpiece << (bitsPerElement - (lendExcl - lbegin)) | rpiece;
        return result;
    }

    //set bits [begin, endExcl( of the 64 bits [l,r] to the bitsPerElement rightmost bits of toSet
    //bits are enumerated from left to right
    void setBits(Data_t& l, Data_t& r, std::size_t begin, std::size_t endExcl, Data_t toSet) const{

        assert(begin < endExcl && endExcl <= 2 * dataTBits());

        toSet = (toSet & mask());
        
        const int rbegin = (begin >= dataTBits()) ? begin - dataTBits() : 0;
        const int rendExcl = (endExcl >= dataTBits()) ? endExcl - dataTBits() : 0;

        const int rpartition = rendExcl - rbegin; // the rightmost rpartition bits of toSet belong to r

        if(rpartition > 0){
            Data_t rmask = 0;
            for(int i = rbegin; i < rendExcl; i++){
                rmask |= (1u << (dataTBits() - 1 - i));
            }

            r &= ~rmask;

            Data_t rpartitionMask = 0;
            for(int i = 0; i < rpartition; i++){
                rpartitionMask = (rpartitionMask << 1) | 1;
            }

            r |= ((toSet & rpartitionMask) << (dataTBits() - rpartition));
        }

        const int lbegin = std::min(dataTBits(), begin);
        const int lendExcl = std::min(dataTBits(), endExcl);

        Data_t lmask = 0;
        for(int i = lbegin; i < lendExcl; i++){
            lmask |= (1u << (dataTBits() - 1 - i));
        }

        l &= ~lmask;

        toSet = toSet >> rpartition;
        toSet = toSet << (dataTBits() - lendExcl);

        l |= toSet;
    }

    void setElement(int index, Element_t element){
        assert(index < size());
        assert(getMinElement() <= element && element <= getMaxElement());

        if(getMinElement() == getMaxElement()){
            return;
        }

        const int diff = element - minElement;
        const std::uint64_t firstBit = bitsPerElement * index;
        const std::uint64_t lastBitExcl = bitsPerElement * (index+1);
        const std::uint64_t firstuintindex = firstBit / dataTBits();
        const int begin = firstBit - firstuintindex * dataTBits();
        const int endExcl = lastBitExcl - firstuintindex * dataTBits();
        //std::cerr << index << ": Set [" << begin << ", " << endExcl << ") of index " << firstuintindex << "\n";

        auto& first = data[firstuintindex];
        //prevent oob access
        auto& second = int(firstuintindex) == numData - 1 ? data[firstuintindex] : data[firstuintindex + 1];

        setBits(first, second, begin, endExcl, Data_t(diff));
    }

    Data_t mask() const noexcept{
        return (bitsPerElement == sizeof(Data_t) * 8) ? std::numeric_limits<Data_t>::max() : ((1 << bitsPerElement) - 1);
    }

    static constexpr std::size_t dataTBits() {
        return 8 * sizeof(Data_t);
    }

    std::uint8_t bitsPerElement = 0;
    Element_t minElement = 0;
    Element_t maxElement = 0;
    int numElements = 0;
    int numData = 0;
    std::unique_ptr<Data_t[]> data;
};

}












#endif