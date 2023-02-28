#ifndef CARE_CONTIGUOUSREADSTORAGE_HPP
#define CARE_CONTIGUOUSREADSTORAGE_HPP

#include <config.hpp>
#include <sequencehelpers.hpp>
#include <lengthstorage.hpp>
#include <memorymanagement.hpp>

#include <threadpool.hpp>
#include <util.hpp>
#include <readstorageconstruction.hpp>
#include <cpureadstorage.hpp>
#include <sharedmutex.hpp>

#include <algorithm>
#include <limits>
#include <cassert>
#include <cstdint>
#include <string>
#include <vector>
#include <omp.h>
#include <map>
#include <fstream>
#include <memory>
#include <cstring>
#include <mutex>
#include <atomic>


namespace care{

namespace cpu{

    struct ContiguousReadStorage : public CpuReadStorage{

        struct Statistics{
            int maximumSequenceLength = 0;
            int minimumSequenceLength = std::numeric_limits<int>::max();

            bool operator==(const Statistics& rhs) const noexcept {
                return maximumSequenceLength == rhs.maximumSequenceLength
                    && minimumSequenceLength == rhs.minimumSequenceLength;
            }
            bool operator!=(const Statistics& rhs) const noexcept{
                return !(operator==(rhs));
            }
        };

        struct GatherHandle{
            std::vector<int> permutation;
        };

        enum class GatherType {Direct, Sorted};

        static constexpr bool has_reverse_complement = false;
        static constexpr int serialization_id = 1;

        using Length_t = int;
        using LengthStore_t = LengthStore<std::uint32_t>;

        std::unique_ptr<unsigned int[]> h_sequence_data = nullptr;
        std::unique_ptr<char[]> h_quality_data = nullptr;
        int sequenceLengthLowerBound = 0;
        int sequenceLengthUpperBound = 0;
        int sequenceDataPitchInInts = 0;
        int sequenceQualitiesPitchInBytes = 0;
        bool useQualityScores = false;
        read_number maximumNumberOfSequences = 0;
        std::size_t sequence_data_bytes = 0;
        std::size_t quality_data_bytes = 0;
        std::vector<read_number> readIdsOfReadsWithUndeterminedBase{}; //sorted in ascending order
        std::mutex mutexUndeterminedBaseReads{};
        Statistics statistics{};
        std::atomic<read_number> numberOfInsertedReads{0};
        LengthStore_t lengthStorage{};

        mutable int counter = 0;
        mutable SharedMutex sharedmutex{};
        //mutable std::vector<std::unique_ptr<TempData>> tempdataVector{};

    public: //inherited interface

    void areSequencesAmbiguous(
        bool* result, 
        const read_number* readIds, 
        int numSequences
    ) const override{
        for(int k = 0; k < numSequences; k++){
            result[k] = readContainsN(readIds[k]);
        }
    }

    void gatherSequences(
        unsigned int* sequence_data,
        size_t outSequencePitchInInts,
        const read_number* readIds,
        int numSequences
    ) const override{
        GatherHandle myhandle{};

        gatherSequenceData(
            myhandle,
            readIds,
            numSequences,
            sequence_data,
            outSequencePitchInInts
        );
    }

    void gatherQualities(
        char* quality_data,
        size_t out_quality_pitch,
        const read_number* readIds,
        int numSequences
    ) const override{
        GatherHandle myhandle{};

        gatherSequenceQualities(
            myhandle,
            readIds,
            numSequences,
            quality_data,
            out_quality_pitch
        );
    }

    void gatherSequenceLengths(
        int* lengths,
        const read_number* readIds,
        int numSequences
    ) const override {
        GatherHandle myhandle{};

        gatherSequenceLengths(
            myhandle,
            readIds,
            numSequences,
            lengths,
            1
        );
    }

    void getIdsOfAmbiguousReads(
        read_number* ids
    ) const override{
        std::copy(readIdsOfReadsWithUndeterminedBase.begin(), readIdsOfReadsWithUndeterminedBase.end(), ids);
    }

    std::int64_t getNumberOfReadsWithN() const override {
        return readIdsOfReadsWithUndeterminedBase.size();
    }

    MemoryUsage getMemoryInfo() const override {
        MemoryUsage info;
        info.host = sequence_data_bytes;
        info.host += lengthStorage.getRawSizeInBytes();
        info.host += sizeof(read_number) * readIdsOfReadsWithUndeterminedBase.capacity();

        if(useQualityScores){
            info.host += quality_data_bytes;
        }         

        return info;
    }

    read_number getNumberOfReads() const override {
        return numberOfInsertedReads;
    }

    bool canUseQualityScores() const override {
        return useQualityScores;
    }

    int getSequenceLengthLowerBound() const override {
        return sequenceLengthLowerBound;
    }

    int getSequenceLengthUpperBound() const override {
        return sequenceLengthUpperBound;
    }

    void destroy() {
        auto deallocVector = [](auto& vec){
            using T = typename std::remove_reference<decltype(vec)>::type;
            T tmp{};
            vec.swap(tmp);
        };
        
        h_sequence_data.reset();
        h_quality_data.reset();
        lengthStorage.destroy();

        deallocVector(readIdsOfReadsWithUndeterminedBase);

    }

    public:

        const unsigned int* getSequenceArray() const noexcept{
            return h_sequence_data.get();
        }

        const char* getQualityArray() const noexcept{
            return h_quality_data.get();
        }

        const LengthStore_t& getLengthStore() const noexcept{
            return lengthStorage;
        }

        std::size_t getSequencePitch() const noexcept{
            return sequenceDataPitchInInts * sizeof(unsigned int);
        }

        std::size_t getQualityPitch() const noexcept{
            return sequenceQualitiesPitchInBytes;
        }

        const read_number* getAmbiguousIds() const noexcept{
            return readIdsOfReadsWithUndeterminedBase.data();
        }



        ContiguousReadStorage() : ContiguousReadStorage(0){}

        ContiguousReadStorage(read_number nSequences) : ContiguousReadStorage(nSequences, false){}

        ContiguousReadStorage(read_number nSequences, bool b) : ContiguousReadStorage(nSequences, b, 0, 0){
        }

        ContiguousReadStorage(read_number nSequences, bool b, int minimum_sequence_length, int maximum_sequence_length){
            init(nSequences, b, minimum_sequence_length, maximum_sequence_length);
        }

        void init(read_number nSequences, bool b, int minimum_sequence_length, int maximum_sequence_length){
            sequenceLengthLowerBound = minimum_sequence_length,
            sequenceLengthUpperBound = maximum_sequence_length,
            sequenceDataPitchInInts = SequenceHelpers::getEncodedNumInts2Bit(maximum_sequence_length),
            sequenceQualitiesPitchInBytes = maximum_sequence_length,
            useQualityScores = b,
            maximumNumberOfSequences = nSequences;

            lengthStorage = std::move(LengthStore_t(sequenceLengthLowerBound, sequenceLengthUpperBound, nSequences));

            h_sequence_data.reset(new unsigned int[std::size_t(maximumNumberOfSequences) * sequenceDataPitchInInts]);
            sequence_data_bytes = sizeof(unsigned int) * std::size_t(maximumNumberOfSequences) * sequenceDataPitchInInts;

            if(useQualityScores){
                h_quality_data.reset(new char[std::size_t(maximumNumberOfSequences) * sequenceQualitiesPitchInBytes]);
                quality_data_bytes = sizeof(char) * std::size_t(maximumNumberOfSequences) * sequenceQualitiesPitchInBytes;
            }

            std::fill_n(&h_sequence_data[0], std::size_t(maximumNumberOfSequences) * sequenceDataPitchInInts, 0);
            std::fill(&h_quality_data[0], &h_quality_data[quality_data_bytes], 0);
        }

        ContiguousReadStorage(const ContiguousReadStorage& other) = delete;
        ContiguousReadStorage& operator=(const ContiguousReadStorage& other) = delete;

        ContiguousReadStorage(ContiguousReadStorage&& other)
            : h_sequence_data(std::move(other.h_sequence_data)),
              h_quality_data(std::move(other.h_quality_data)),
              sequenceLengthLowerBound(other.sequenceLengthLowerBound),
              sequenceLengthUpperBound(other.sequenceLengthUpperBound),
              sequenceDataPitchInInts(other.sequenceDataPitchInInts),
              sequenceQualitiesPitchInBytes(other.sequenceQualitiesPitchInBytes),
              useQualityScores(other.useQualityScores),
              maximumNumberOfSequences(other.maximumNumberOfSequences),
              sequence_data_bytes(other.sequence_data_bytes),
              quality_data_bytes(other.quality_data_bytes),
              readIdsOfReadsWithUndeterminedBase(std::move(other.readIdsOfReadsWithUndeterminedBase)),
              statistics(std::move(other.statistics)),
              numberOfInsertedReads(other.numberOfInsertedReads.load()),
              lengthStorage(std::move(other.lengthStorage)){

            other.numberOfInsertedReads = 0;
            other.statistics = Statistics{};

        }

        ContiguousReadStorage& operator=(ContiguousReadStorage&& other){
            h_sequence_data = std::move(other.h_sequence_data);
            h_quality_data = std::move(other.h_quality_data);
            sequenceLengthLowerBound = other.sequenceLengthLowerBound;
            sequenceLengthUpperBound = other.sequenceLengthUpperBound;
            sequenceDataPitchInInts = other.sequenceDataPitchInInts;
            sequenceQualitiesPitchInBytes = other.sequenceQualitiesPitchInBytes;
            useQualityScores = other.useQualityScores;
            maximumNumberOfSequences = other.maximumNumberOfSequences;
            sequence_data_bytes = other.sequence_data_bytes;
            quality_data_bytes = other.quality_data_bytes;
            readIdsOfReadsWithUndeterminedBase = std::move(other.readIdsOfReadsWithUndeterminedBase);
            statistics = std::move(other.statistics);
            numberOfInsertedReads = other.numberOfInsertedReads.load();
            lengthStorage = std::move(other.lengthStorage);

            other.numberOfInsertedReads = 0;
            other.statistics = Statistics{};

            return *this;
        }

        bool operator==(const ContiguousReadStorage& other) const{
            if(sequenceLengthLowerBound != other.sequenceLengthLowerBound)
                return false;
            if(sequenceLengthUpperBound != other.sequenceLengthUpperBound)
                return false;
            if(sequenceDataPitchInInts != other.sequenceDataPitchInInts)
                return false;
            if(sequenceQualitiesPitchInBytes != other.sequenceQualitiesPitchInBytes)
                return false;
            if(useQualityScores != other.useQualityScores)
                return false;
            if(maximumNumberOfSequences != other.maximumNumberOfSequences)
                return false;
            if(useQualityScores != other.useQualityScores)
                return false;
            if(sequence_data_bytes != other.sequence_data_bytes)
                return false;
            if(quality_data_bytes != other.quality_data_bytes)
                return false;
            if(readIdsOfReadsWithUndeterminedBase != other.readIdsOfReadsWithUndeterminedBase){
                return false;
            }
            if(statistics != other.statistics){
                return false;
            }

            if(lengthStorage != other.lengthStorage){
                return false;
            }

            if(0 != std::memcmp(h_sequence_data.get(), other.h_sequence_data.get(), sequence_data_bytes))
                return false;
            if(0 != std::memcmp(h_quality_data.get(), other.h_quality_data.get(), quality_data_bytes))
                return false;

            return true;
        }

        bool operator!=(const ContiguousReadStorage& other) const{
            return !(*this == other);
        }

        void construct(
            std::vector<std::string> inputfiles,
            bool useQualityScores,
            read_number expectedNumberOfReads,
            int expectedMinimumReadLength,
            int expectedMaximumReadLength,
            int threads,
            bool showProgress
        ){

            auto makeInserterFunc = [this](){
                return [&, this](ThreadPool* tp, read_number* indices, Read* reads, int count){
                    this->setReads(tp, indices, reads, count);
                };
            };
            auto makeReadContainsNFunc = [this](){
                return [&, this](read_number readId, bool contains){
                    this->setReadContainsN(readId, contains);
                };
            };

            constructReadStorageFromFiles(
                inputfiles,
                useQualityScores,
                expectedNumberOfReads,
                expectedMinimumReadLength,
                expectedMaximumReadLength,
                threads,
                showProgress,
                makeInserterFunc,
                makeReadContainsNFunc
            );
     
        }

        void constructPaired(
            std::vector<std::string> inputfiles,
            bool useQualityScores,
            read_number expectedNumberOfReads,
            int expectedMinimumReadLength,
            int expectedMaximumReadLength,
            int threads,
            bool showProgress
        ){

            auto makeInserterFunc = [this](){
                return [&, this](ThreadPool* tp, read_number* indices, Read* reads, int count){
                    this->setReads(tp, indices, reads, count);
                };
            };
            auto makeReadContainsNFunc = [this](){
                return [&, this](read_number readId, bool contains){
                    this->setReadContainsN(readId, contains);
                };
            };

            constructReadStorageFromPairedEndFiles(
                inputfiles,
                useQualityScores,
                expectedNumberOfReads,
                expectedMinimumReadLength,
                expectedMaximumReadLength,
                threads,
                showProgress,
                makeInserterFunc,
                makeReadContainsNFunc
            );
     
        }

        std::size_t size() const{
            //assert(std::size_t(maximumNumberOfSequences) * maximum_allowed_sequence_bytes == sequence_data_bytes);

            std::size_t result = 0;
            result += sequence_data_bytes;
            result += lengthStorage.getRawSizeInBytes();

            if(useQualityScores){
                //assert(std::size_t(maximumNumberOfSequences) * sequenceLengthUpperBound * sizeof(char) == quality_data_bytes);
                result += quality_data_bytes;
            }

            return result;
        }

    	void resize(read_number nReads){
    		assert(getMaximumNumberOfSequences() >= nReads);

            maximumNumberOfSequences = nReads;
    	}

        Statistics getStatistics() const{
            return statistics;
        }

        

        void setReadContainsN(read_number readId, bool contains){

            std::lock_guard<std::mutex> l(mutexUndeterminedBaseReads);

            auto pos = std::lower_bound(readIdsOfReadsWithUndeterminedBase.begin(),
                                                readIdsOfReadsWithUndeterminedBase.end(),
                                                readId);

            if(contains){
                //if readId is not already in the vector, insert it
                if((pos == readIdsOfReadsWithUndeterminedBase.end()) || (pos != readIdsOfReadsWithUndeterminedBase.end() && *pos != readId)){                    
                    readIdsOfReadsWithUndeterminedBase.insert(pos, readId);
                }
            }else{
                if(pos != readIdsOfReadsWithUndeterminedBase.end() && *pos == readId){
                    //remove mark
                    readIdsOfReadsWithUndeterminedBase.erase(pos);
                }
            }
        }

        bool readContainsN(read_number readId) const{

            auto pos = std::lower_bound(readIdsOfReadsWithUndeterminedBase.begin(),
                                                readIdsOfReadsWithUndeterminedBase.end(),
                                                readId);
            bool b2 = readIdsOfReadsWithUndeterminedBase.end() != pos && *pos == readId;

            return b2;
        }

        void printAmbig(){
            for(auto x : readIdsOfReadsWithUndeterminedBase){
                std::cerr << x << " ";
            }
            std::cerr << "\n";
        }

        

private:
        void insertSequence(read_number readNumber, const std::string& sequence){
            auto identity = [](auto i){return i;};

            const int sequencelength = sequence.length();

            unsigned int* dest = &h_sequence_data[std::size_t(readNumber) * sequenceDataPitchInInts];
            SequenceHelpers::encodeSequence2Bit(
                dest,
                sequence.c_str(),
                sequence.length(),
                identity
            );

            statistics.minimumSequenceLength = std::min(statistics.minimumSequenceLength, sequencelength);
            statistics.maximumSequenceLength = std::max(statistics.maximumSequenceLength, sequencelength);

            lengthStorage.setLength(readNumber, Length_t(sequence.length()));

            read_number prev_value = numberOfInsertedReads;
            while(prev_value < readNumber+1 && !numberOfInsertedReads.compare_exchange_weak(prev_value, readNumber+1)){
                ;
            }
        }
public:

        void setReads(
            ThreadPool* threadPool, 
            const read_number* indices, 
            const Read* reads, 
            int numReads
        ){
            if(numReads == 0) return;
            
            //TIMERSTARTCPU(internalinit);
            #ifndef NDEBUG
            auto lengthInRange = [&](Length_t length){
                return getSequenceLengthLowerBound() <= length && length <= getSequenceLengthUpperBound();
            };
            assert(numReads > 0);
            assert(std::all_of(indices, indices + numReads, [&](auto i){ return i < getMaximumNumberOfSequences();}));
            assert(std::all_of(reads, reads + numReads, [&](const auto& r){ return lengthInRange(Length_t(r.sequence.length()));}));
            #endif
            
            if(canUseQualityScores()){
                assert(std::all_of(reads, reads + numReads, [&](const auto& r){ return r.sequence.length() == r.quality.length();}));
            }

            auto minmax = std::minmax_element(reads, reads + numReads, [](const auto& r1, const auto& r2){
                return r1.sequence.length() < r2.sequence.length();
            });

            statistics.minimumSequenceLength = std::min(statistics.minimumSequenceLength, int(minmax.first->sequence.length()));
            statistics.maximumSequenceLength = std::max(statistics.maximumSequenceLength, int(minmax.second->sequence.length()));

            read_number maxIndex = *std::max_element(indices, indices + numReads);

            read_number prev_value = numberOfInsertedReads;
            while(prev_value < maxIndex+1 && !numberOfInsertedReads.compare_exchange_weak(prev_value, maxIndex+1)){
                ;
            }


            ThreadPool::ParallelForHandle pforHandle;

            auto body = [&](auto begin, auto end, int /*threadid*/){
                for(auto i = begin; i < end; i++){
                    const auto& myRead = reads[i];
                    const read_number myReadNumber = indices[i];

                    insertRead(myReadNumber, myRead.sequence, myRead.quality);
                }
            };

            threadPool->parallelFor(pforHandle, 0, numReads, body);

        }


        void insertRead(read_number readNumber, const std::string& sequence){
            assert(readNumber < getMaximumNumberOfSequences());
            assert(int(sequence.length()) <= sequenceLengthUpperBound);

    		if(useQualityScores){
    			insertRead(readNumber, sequence, std::string(sequence.length(), 'A'));
    		}else{
    			insertSequence(readNumber, sequence);
    		}
    	}

        void insertRead(read_number readNumber, const std::string& sequence, const std::string& quality){
            assert(readNumber < getMaximumNumberOfSequences());
            assert(int(sequence.length()) <= sequenceLengthUpperBound);
            
            

    		insertSequence(readNumber, sequence);

    		if(useQualityScores){
                assert(int(quality.length()) <= sequenceLengthUpperBound);
                assert(sequence.length() == quality.length());
                
                std::memcpy(&h_quality_data[std::size_t(readNumber) * std::size_t(sequenceLengthUpperBound)],
                            quality.c_str(),
                            sizeof(char) * quality.length());
    		}
    	}

        const char* fetchQuality_ptr(read_number readNumber) const{
            if(useQualityScores){
                return &h_quality_data[std::size_t(readNumber) * std::size_t(sequenceLengthUpperBound)];
            }else{
                return nullptr;
            }
        }

        const unsigned int* fetchSequenceData_ptr(read_number readNumber) const{
        	return &h_sequence_data[std::size_t(readNumber) * std::size_t(sequenceDataPitchInInts)];
        }

        int fetchSequenceLength(read_number readNumber) const{
            return lengthStorage.getLength(readNumber);
        }

        template<class T, class GatherType>
        void gatherImpl(
                GatherHandle& handle,
                GatherType gatherType,
                const T* source,
                size_t sourcePitchElements,
                const read_number* readIds,
                int numReadIds,
                T* destination,
                size_t destinationPitchElements) const noexcept{
            
            if(numReadIds == 0){
                return;
            }

            if(gatherType == GatherType::Sorted){
                handle.permutation.resize(numReadIds);
                //handle.data.resize(sourcePitchElement * sizeof(T) * numReadIds);

                std::iota(
                    handle.permutation.begin(), 
                    handle.permutation.end(),
                    0
                );

                std::sort(
                    handle.permutation.begin(), 
                    handle.permutation.end(),
                    [&](const auto& l, const auto& r){
                        return readIds[l] < readIds[r];
                    }
                );
            }

            constexpr int prefetch_distance = 4;

            for(int i = 0; i < numReadIds && i < prefetch_distance; ++i) {
                const int index = gatherType == GatherType::Sorted ? handle.permutation[i] : i;
                const read_number nextReadId = readIds[index];
                const T* const nextData = source + sourcePitchElements * nextReadId;
                __builtin_prefetch(nextData, 0, 0);
            }

            for(int i = 0; i < numReadIds; i++){
                if(i + prefetch_distance < numReadIds) {
                    const int index = gatherType == GatherType::Sorted ? handle.permutation[i + prefetch_distance] : i + prefetch_distance;
                    const read_number nextReadId = readIds[index];
                    const T* const nextData = source + sourcePitchElements * nextReadId;
                    __builtin_prefetch(nextData, 0, 0);
                }

                const int index = gatherType == GatherType::Sorted ? handle.permutation[i] : i;
                const read_number readId = readIds[index];
                const T* const data = source + sourcePitchElements * readId;

                T* const destData = destination + destinationPitchElements * index;
                std::copy_n(data, sourcePitchElements, destData);
            }
        }


        void gatherSequenceData(
                GatherHandle& handle,
                const read_number* readIds,
                int numReadIds,
                unsigned int* destination,
                int destinationPitchElements) const noexcept{

            gatherImpl(
                handle,
                GatherType::Direct,
                &h_sequence_data[0],
                sequenceDataPitchInInts,
                readIds,
                numReadIds,
                destination,
                destinationPitchElements
            );
        }

        void gatherSequenceDataSpecial(
                GatherHandle& handle,
                const read_number* readIds,
                int numReadIds,
                unsigned int* destination,
                int destinationPitchElements) const noexcept{

            gatherImpl(
                handle,
                GatherType::Sorted,
                &h_sequence_data[0],
                sequenceDataPitchInInts,
                readIds,
                numReadIds,
                destination,
                destinationPitchElements
            );
        }

        void gatherSequenceQualities(
                GatherHandle& handle,
                const read_number* readIds,
                int numReadIds,
                char* destination,
                int destinationPitchElements) const noexcept{

            gatherImpl(
                handle,
                GatherType::Direct,
                &h_quality_data[0],
                sequenceQualitiesPitchInBytes,
                readIds,
                numReadIds,
                destination,
                destinationPitchElements
            );
        }

        void gatherSequenceQualitiesSpecial(
                GatherHandle& handle,
                const read_number* readIds,
                int numReadIds,
                char* destination,
                int destinationPitchElements) const noexcept{

            gatherImpl(
                handle,
                GatherType::Sorted,
                &h_quality_data[0],
                sequenceQualitiesPitchInBytes,
                readIds,
                numReadIds,
                destination,
                destinationPitchElements
            );
        }

        void gatherSequenceLengths(
                GatherHandle& /*handle*/,
                const read_number* readIds,
                int numReadIds,
                int* destination,
                int destinationPitchElements = 1) const noexcept{

            for(int i = 0; i < numReadIds; i++){
                int* const destLength = destination + i * destinationPitchElements;
                *destLength = fetchSequenceLength(readIds[i]);
            }            
        }

        void gatherSequenceLengthsSpecial(
                GatherHandle& handle,
                const read_number* readIds,
                int numReadIds,
                int* destination,
                int destinationPitchElements = 1) const noexcept{

            gatherSequenceLengths(
                handle,
                readIds,
                numReadIds,
                destination,
                destinationPitchElements
            );
        }

       	std::uint64_t getMaximumNumberOfSequences() const{
    		return maximumNumberOfSequences;
    	}

        int getMaximumAllowedSequenceBytes() const{
            return sequenceDataPitchInInts * sizeof(unsigned int);
        }



        void saveToFile(const std::string& filename) const{
            std::ofstream stream(filename, std::ios::binary);

            read_number inserted = getNumberOfReads();

            stream.write(reinterpret_cast<const char*>(&inserted), sizeof(read_number));
            stream.write(reinterpret_cast<const char*>(&sequenceLengthUpperBound), sizeof(int));
            stream.write(reinterpret_cast<const char*>(&sequenceLengthLowerBound), sizeof(int));
            stream.write(reinterpret_cast<const char*>(&useQualityScores), sizeof(bool));
            stream.write(reinterpret_cast<const char*>(&statistics), sizeof(Statistics));

            auto pos = stream.tellp();

            stream.seekp(sizeof(std::size_t) * 4, std::ios_base::cur);

            std::size_t lengthsBytes = lengthStorage.writeToStream(stream);
            std::size_t sequencesBytes = sequence_data_bytes + sizeof(std::size_t);      
            std::size_t qualitiesBytes = quality_data_bytes + sizeof(std::size_t);

            //sequence data
            stream.write(reinterpret_cast<const char*>(&sequence_data_bytes), sizeof(std::size_t));
            stream.write(reinterpret_cast<const char*>(&h_sequence_data[0]), sequence_data_bytes); 
            //quality data
            stream.write(reinterpret_cast<const char*>(&quality_data_bytes), sizeof(std::size_t));
            stream.write(reinterpret_cast<const char*>(&h_quality_data[0]), quality_data_bytes);

            std::size_t numUndeterminedReads = readIdsOfReadsWithUndeterminedBase.size();
            stream.write(reinterpret_cast<const char*>(&numUndeterminedReads), sizeof(size_t));
            stream.write(reinterpret_cast<const char*>(readIdsOfReadsWithUndeterminedBase.data()), numUndeterminedReads * sizeof(read_number));

            std::size_t ambigBytes = sizeof(std::size_t) + numUndeterminedReads * sizeof(read_number);

            stream.seekp(pos);
            stream.write(reinterpret_cast<const char*>(&lengthsBytes), sizeof(std::size_t));
            stream.write(reinterpret_cast<const char*>(&sequencesBytes), sizeof(std::size_t));
            stream.write(reinterpret_cast<const char*>(&qualitiesBytes), sizeof(std::size_t));
            stream.write(reinterpret_cast<const char*>(&ambigBytes), sizeof(std::size_t));

        }

        void loadFromFile(const std::string& filename){
            std::ifstream stream(filename, std::ios::binary);
            if(!stream)
                throw std::runtime_error("Cannot open file " + filename);

            destroy();

            read_number loaded_inserted = 0;
            int loaded_sequenceLengthLowerBound = 0;
            int loaded_sequenceLengthUpperBound = 0;
            bool loaded_hasQualityScores = false;

            std::size_t loaded_sequence_data_bytes = 0;
            std::size_t loaded_quality_data_bytes = 0;            

            stream.read(reinterpret_cast<char*>(&loaded_inserted), sizeof(read_number));
            stream.read(reinterpret_cast<char*>(&loaded_sequenceLengthLowerBound), sizeof(int));
            stream.read(reinterpret_cast<char*>(&loaded_sequenceLengthUpperBound), sizeof(int));            
            stream.read(reinterpret_cast<char*>(&loaded_hasQualityScores), sizeof(bool));

            init(loaded_inserted, canUseQualityScores(), loaded_sequenceLengthLowerBound, loaded_sequenceLengthUpperBound);

            numberOfInsertedReads = loaded_inserted;

            stream.read(reinterpret_cast<char*>(&statistics), sizeof(Statistics));

            std::size_t lengthsBytes = 0;
            std::size_t sequencesBytes = 0;      
            std::size_t qualitiesBytes = 0;       
            std::size_t ambigBytes = 0;

            stream.read(reinterpret_cast<char*>(&lengthsBytes), sizeof(std::size_t));
            stream.read(reinterpret_cast<char*>(&sequencesBytes), sizeof(std::size_t));
            stream.read(reinterpret_cast<char*>(&qualitiesBytes), sizeof(std::size_t));
            stream.read(reinterpret_cast<char*>(&ambigBytes), sizeof(std::size_t));

            lengthStorage.readFromStream(stream);

            stream.read(reinterpret_cast<char*>(&loaded_sequence_data_bytes), sizeof(std::size_t));
            stream.read(reinterpret_cast<char*>(&h_sequence_data[0]), loaded_sequence_data_bytes);

            if(canUseQualityScores() && loaded_hasQualityScores){
                //std::cerr << "load qualities\n";

                stream.read(reinterpret_cast<char*>(&loaded_quality_data_bytes), sizeof(std::size_t));            
                stream.read(reinterpret_cast<char*>(&h_quality_data[0]), loaded_quality_data_bytes);
            }else if(canUseQualityScores() && !loaded_hasQualityScores){
                    //std::cerr << "no q in bin file\n";
                    throw std::runtime_error("Quality scores expected in preprocessed reads file to load, but none are present. Abort.");
            }else if(!canUseQualityScores() && loaded_hasQualityScores){
                    //std::cerr << "skip qualities\n";
                    stream.ignore(qualitiesBytes);
            }else{
                //!canUseQualityScores() && !loaded_hasQualityScores
                //std::cerr << "no q in file, and no q required. Ok\n";
                stream.ignore(qualitiesBytes);
            }
            

            std::size_t numUndeterminedReads = 0;
            stream.read(reinterpret_cast<char*>(&numUndeterminedReads), sizeof(std::size_t));
            readIdsOfReadsWithUndeterminedBase.resize(numUndeterminedReads);
            stream.read(reinterpret_cast<char*>(readIdsOfReadsWithUndeterminedBase.data()), numUndeterminedReads * sizeof(read_number));
        }
    };



}



}

#endif
