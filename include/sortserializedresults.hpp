#ifndef CARE_SORTSERIALIZEDRESULTS_HPP
#define CARE_SORTSERIALIZEDRESULTS_HPP

#include <config.hpp>

#include <serializedobjectstorage.hpp>
#include <sortbygeneratedkeys.hpp>

#include <cstdint>
#include <iostream>
#include <limits>
#include <algorithm>

namespace care{

    template<class T> // T type of serialized objects
    void sortSerializedResultsByReadIdAscending(
        SerializedObjectStorage& partialResults,
        std::size_t memoryForSortingInBytes
    ){
        if(partialResults.size() < 2) return;

        //return read id of the object serialized at ptr
        auto extractKey = [](const std::uint8_t* ptr){
            const read_number id = T::parseReadId(ptr);            
            return id;
        };

        //return read id of i-th object serialized in partialResults
        auto keyGenerator = [&](std::size_t i){
            return extractKey(partialResults.getPointer(i));
        };

        bool sortValuesSuccess = false;
      
        try{
            if(partialResults.size() <= std::size_t(std::numeric_limits<std::uint32_t>::max())){
                //std::cerr << "sortValuesByGeneratedKeys<std::uint32_t>\n";
                sortValuesSuccess = sortValuesByGeneratedKeys<std::uint32_t>(
                    memoryForSortingInBytes,
                    partialResults.getOffsetBuffer(),
                    partialResults.size(),
                    keyGenerator
                );
            }else{
                //std::cerr << "sortValuesByGeneratedKeys<std::uint64_t>\n";
                sortValuesSuccess = sortValuesByGeneratedKeys<std::uint64_t>(
                    memoryForSortingInBytes,
                    partialResults.getOffsetBuffer(),
                    partialResults.size(),
                    keyGenerator
                );
            }

        } catch (...){
            std::cerr << "Final fallback\n";
        }

        if(!sortValuesSuccess){        
            auto offsetcomparator = [&](std::size_t elementOffset1, std::size_t elementOffset2){
                const std::uint8_t* ptr1 = partialResults.getDataBuffer() + elementOffset1;
                const std::uint8_t* ptr2 = partialResults.getDataBuffer() + elementOffset2;
                const read_number key1 = extractKey(ptr1);
                const read_number key2 = extractKey(ptr2);

                return key1 < key2;
            };

            std::sort(
                partialResults.getOffsetBuffer(), 
                partialResults.getOffsetBuffer() + partialResults.size(), 
                offsetcomparator
            );
        }
    }

}

#endif