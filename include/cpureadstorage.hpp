#ifndef CARE_CPUREADSTORAGE_HPP
#define CARE_CPUREADSTORAGE_HPP

#include <config.hpp>
#include <memorymanagement.hpp>

#include <cstdint>

namespace care{

class CpuReadStorage{
public:

    virtual ~CpuReadStorage() = default;

    virtual void areSequencesAmbiguous(
        bool* result, 
        const read_number* readIds, 
        int numSequences
    ) const = 0;

    virtual void gatherSequences(
        unsigned int* sequence_data,
        std::size_t outSequencePitchInInts,
        const read_number* readIds,
        int numSequences
    ) const = 0;

    virtual void gatherContiguousSequences(
        unsigned int* sequence_data,
        std::size_t outSequencePitchInInts,
        read_number firstIndex,
        int numSequences
    ) const = 0;

    virtual void gatherQualities(
        char* quality_data,
        std::size_t out_quality_pitch,
        const read_number* readIds,
        int numSequences
    ) const = 0;

    virtual void gatherEncodedQualities(
        unsigned int* encodedQualities,
        std::size_t outputPitchInInts,
        const read_number* readIds,
        int numSequences
    ) const = 0;

    virtual void gatherContiguousEncodedQualities(
        unsigned int* encodedQualities,
        std::size_t outputPitchInInts,
        read_number firstIndex,
        int numSequences
    ) const = 0;

    virtual void gatherSequenceLengths(
        int* lengths,
        const read_number* readIds,
        int numSequences
    ) const = 0;

    virtual void getIdsOfAmbiguousReads(
        read_number* ids
    ) const = 0;

    virtual std::int64_t getNumberOfReadsWithN() const = 0;

    virtual MemoryUsage getMemoryInfo() const = 0;

    virtual read_number getNumberOfReads() const = 0;

    virtual bool canUseQualityScores() const = 0;

    virtual int getSequenceLengthLowerBound() const = 0;

    virtual int getSequenceLengthUpperBound() const = 0;

    virtual bool isPairedEnd() const = 0;

    virtual int getQualityBits() const = 0;

    //virtual void destroy() = 0;
};



}

#endif