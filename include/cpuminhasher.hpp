#ifndef CARE_CPUMINHASHER_HPP
#define CARE_CPUMINHASHER_HPP

#include <config.hpp>
#include <memorymanagement.hpp>
#include <minhasherhandle.hpp>
#include <threadpool.hpp>

#include <cstdint>
#include <vector>

namespace care{

class CpuMinhasher{
public:

    using Key = kmer_type;

    virtual ~CpuMinhasher() = default;

    virtual MinhasherHandle makeMinhasherHandle() const = 0;

    virtual void destroyHandle(MinhasherHandle& handle) const = 0;

    //construction
    virtual void setHostMemoryLimitForConstruction(std::size_t bytes) = 0;
    virtual void setDeviceMemoryLimitsForConstruction(const std::vector<std::size_t>&) = 0;
    virtual int addHashTables(int numAdditionalTables, const int* hashFunctionIds) = 0;
    virtual void setThreadPool(ThreadPool* tp) = 0;
    virtual void compact() = 0;
    virtual void constructionIsFinished() = 0;

    virtual void insert(
        const unsigned int* h_sequenceData2Bit,
        int numSequences,
        const int* h_sequenceLengths,
        std::size_t encodedSequencePitchInInts,
        const read_number* h_readIds,
        int firstHashfunction,
        int numHashfunctions,
        const int* h_hashFunctionNumbers
    ) = 0;

    //return number of hash tables where insert was unsuccessfull
    virtual int checkInsertionErrors(
        int firstHashfunction,
        int numHashfunctions
    ) = 0;

    virtual void determineNumValues(
        MinhasherHandle& queryHandle,
        const unsigned int* h_sequenceData2Bit,
        std::size_t encodedSequencePitchInInts,
        const int* h_sequenceLengths,
        int numSequences,
        int* h_numValuesPerSequence,
        int& totalNumValues
    ) const = 0;

    virtual void retrieveValues(
        MinhasherHandle& queryHandle,
        int numSequences,
        int totalNumValues,
        read_number* h_values,
        const int* h_numValuesPerSequence,
        int* h_offsets //numSequences + 1
    ) const = 0;

    //virtual void compact() = 0;

    virtual MemoryUsage getMemoryInfo() const noexcept = 0;

    virtual MemoryUsage getMemoryInfo(const MinhasherHandle& handle) const noexcept = 0;

    virtual int getNumResultsPerMapThreshold() const noexcept = 0;
    
    virtual int getNumberOfMaps() const noexcept = 0;

    virtual int getKmerSize() const noexcept = 0;

    //virtual void destroy() = 0;

    virtual void writeToStream(std::ostream& os) const = 0;
    virtual int loadFromStream(std::ifstream& is, int numMapsUpperLimit) = 0;
    virtual bool canWriteToStream() const noexcept = 0;
    virtual bool canLoadFromStream() const noexcept = 0;

protected:
    MinhasherHandle constructHandle(int id) const{
        return MinhasherHandle{id};
    }

};



}





#endif