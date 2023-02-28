#if 0

#ifndef CARE_SINGLEHASHCPUMINHASHER_HPP
#define CARE_SINGLEHASHCPUMINHASHER_HPP


#include <cpuminhasher.hpp>
#include <cpureadstorage.hpp>
#include <groupbykey.hpp>
#include <cpusequencehasher.hpp>

#include <config.hpp>

#include <cpuhashtable.hpp>

#include <options.hpp>
#include <util.hpp>
#include <hpc_helpers.cuh>
#include <filehelpers.hpp>
#include <sequencehelpers.hpp>
#include <memorymanagement.hpp>
#include <threadpool.hpp>
#include <sharedmutex.hpp>


#include <cassert>
#include <array>
#include <vector>
#include <memory>
#include <limits>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>

namespace care{


    class SingleHashCpuMinhasher : public CpuMinhasher{
    public:
        using Key_t = CpuMinhasher::Key;
        using Value_t = read_number;
    private:
        using HashTable = CpuReadOnlyMultiValueHashTable<kmer_type, read_number>;
        using Range_t = std::pair<const Value_t*, const Value_t*>;

        struct QueryData{

            enum class Stage{
                None,
                NumValues,
                Retrieve
            };

            Stage previousStage = Stage::None;
            std::vector<Range_t> ranges{};
            SetUnionHandle suHandle{};

            MemoryUsage getMemoryInfo() const{
                MemoryUsage info{};
                info.host += sizeof(Range_t) * ranges.capacity();
    
                return info;
            }

            void destroy(){
            }
        };

        struct KVMetadata{
            read_number offset;
            BucketSize numValues;

            bool operator==(const KVMetadata& rhs) const noexcept{
                if(offset != rhs.offset) return false;
                if(numValues != rhs.numValues) return false;
                return true;
            }

            bool operator!=(const KVMetadata& rhs) const noexcept{
                return !operator==(rhs);
            }
        };

        
    public:

        SingleHashCpuMinhasher() : SingleHashCpuMinhasher(0, 255, 16, 0.8f){

        }

        SingleHashCpuMinhasher(int maxNumKeys_, int maxValuesPerKey, int k, float loadfactor_)
            : loadfactor(loadfactor_), maxNumKeys(maxNumKeys_), kmerSize(k), resultsPerMapThreshold(maxValuesPerKey){

        }

        SingleHashCpuMinhasher(const SingleHashCpuMinhasher&) = delete;
        SingleHashCpuMinhasher(SingleHashCpuMinhasher&&) = default;
        SingleHashCpuMinhasher& operator=(const SingleHashCpuMinhasher&) = delete;
        SingleHashCpuMinhasher& operator=(SingleHashCpuMinhasher&&) = default;        


        void constructFromReadStorage(
            const ProgramOptions& programOptions,
            std::uint64_t /*nReads*/,
            const CpuReadStorage& cpuReadStorage
        ){
            auto& readStorage = cpuReadStorage;

            const int requestedNumberOfMaps = programOptions.numHashFunctions;
            numSmallest = requestedNumberOfMaps;

            const std::uint64_t numReads = readStorage.getNumberOfReads();
            const int maximumSequenceLength = readStorage.getSequenceLengthUpperBound();
            const std::size_t encodedSequencePitchInInts = SequenceHelpers::getEncodedNumInts2Bit(maximumSequenceLength);

            const MemoryUsage memoryUsageOfReadStorage = readStorage.getMemoryInfo();
            std::size_t totalLimit = programOptions.memoryTotalLimit;
            if(totalLimit > memoryUsageOfReadStorage.host){
                totalLimit -= memoryUsageOfReadStorage.host;
            }else{
                totalLimit = 0;
            }
            if(totalLimit == 0){
                throw std::runtime_error("Not enough memory available for hash tables. Abort!");
            }
            std::size_t maxMemoryForTables = getAvailableMemoryInKB() * 1024;
            // std::cerr << "available: " << maxMemoryForTables 
            //         << ",memoryForHashtables: " << programOptions.memoryForHashtables
            //         << ", memoryTotalLimit: " << programOptions.memoryTotalLimit
            //         << ", rsHostUsage: " << memoryUsageOfReadStorage.host << "\n";

            maxMemoryForTables = std::min(maxMemoryForTables, 
                                    std::min(programOptions.memoryForHashtables, totalLimit));

            std::cerr << "maxMemoryForTables = " << maxMemoryForTables << " bytes\n";

            setMemoryLimitForConstruction(maxMemoryForTables);

            kvtable = std::make_unique<DoublePassMultiValueHashTable<kmer_type, read_number>>((numReads * numSmallest) / 8, loadfactor);

            ThreadPool tpForHashing(programOptions.threads);
            setThreadPool(&tpForHashing);


            std::size_t numKeys = 0;

            constexpr int batchsize = 1000000;
            //constexpr int batchsize = 20;
            const int numIterations = SDIV(numReads, batchsize);

            std::vector<read_number> currentReadIds(batchsize);
            std::vector<unsigned int> sequencedata(batchsize * encodedSequencePitchInInts);
            std::vector<int> sequencelengths(batchsize);

            helpers::CpuTimer firstpasstimer("firstpass");

            for(int iteration = 0; iteration < numIterations; iteration++){
                const read_number beginid = iteration * batchsize;
                const read_number endid = std::min((iteration + 1) * batchsize, int(numReads));
                const read_number currentbatchsize = endid - beginid;

                std::iota(currentReadIds.begin(), currentReadIds.end(), beginid);

                readStorage.gatherSequences(
                    sequencedata.data(),
                    encodedSequencePitchInInts,
                    currentReadIds.data(),
                    currentbatchsize
                );

                readStorage.gatherSequenceLengths(
                    sequencelengths.data(),
                    currentReadIds.data(),
                    currentbatchsize
                );

                std::vector<kmer_type> tmpkeys(currentbatchsize * numSmallest);
                std::vector<read_number> tmpids(currentbatchsize * numSmallest);

                auto numNewKeys = computeHashesAndReadIds(
                    tmpkeys.data(),
                    tmpids.data(),
                    sequencedata.data(),
                    currentbatchsize,
                    sequencelengths.data(),
                    encodedSequencePitchInInts,
                    currentReadIds.data()
                );

                numKeys += numNewKeys;
                tmpkeys.erase(tmpkeys.begin() + numNewKeys, tmpkeys.end());
                tmpids.erase(tmpids.begin() + numNewKeys, tmpids.end());

                kvtable->firstPassInsert(tmpkeys.data(), tmpids.data(), tmpkeys.size());
            }

            firstpasstimer.print();


            //kvtable->firstPassDone(2, 75);
            kvtable->firstPassDone(2, resultsPerMapThreshold);

            helpers::CpuTimer secondPassTimer("secondpass");

            for(int iteration = 0; iteration < numIterations; iteration++){
                const read_number beginid = iteration * batchsize;
                const read_number endid = std::min((iteration + 1) * batchsize, int(numReads));
                const read_number currentbatchsize = endid - beginid;

                std::iota(currentReadIds.begin(), currentReadIds.end(), beginid);

                readStorage.gatherSequences(
                    sequencedata.data(),
                    encodedSequencePitchInInts,
                    currentReadIds.data(),
                    currentbatchsize
                );

                readStorage.gatherSequenceLengths(
                    sequencelengths.data(),
                    currentReadIds.data(),
                    currentbatchsize
                );

                std::vector<kmer_type> tmpkeys(currentbatchsize * numSmallest);
                std::vector<read_number> tmpids(currentbatchsize * numSmallest);

                auto numNewKeys = computeHashesAndReadIds(
                    tmpkeys.data(),
                    tmpids.data(),
                    sequencedata.data(),
                    currentbatchsize,
                    sequencelengths.data(),
                    encodedSequencePitchInInts,
                    currentReadIds.data()
                );

                numKeys += numNewKeys;
                tmpkeys.erase(tmpkeys.begin() + numNewKeys, tmpkeys.end());
                tmpids.erase(tmpids.begin() + numNewKeys, tmpids.end());

                kvtable->secondPassInsert(tmpkeys.data(), tmpids.data(), tmpkeys.size());
            }

            secondPassTimer.print();

            kvtable->secondPassDone();

            setThreadPool(nullptr);
        }
 

        MinhasherHandle makeMinhasherHandle() const override {
            auto data = std::make_unique<QueryData>();

            std::unique_lock<SharedMutex> lock(sharedmutex);
            const int handleid = counter++;
            MinhasherHandle h = constructHandle(handleid);

            tempdataVector.emplace_back(std::move(data));

            return h;
        }

        void destroyHandle(MinhasherHandle& handle) const override{
            std::unique_lock<SharedMutex> lock(sharedmutex);

            const int id = handle.getId();
            assert(id < int(tempdataVector.size()));
            
            tempdataVector[id] = nullptr;
            handle = constructHandle(std::numeric_limits<int>::max());
        }

        void determineNumValues(
            MinhasherHandle& queryHandle,
            const unsigned int* h_sequenceData2Bit,
            std::size_t encodedSequencePitchInInts,
            const int* h_sequenceLengths,
            int numSequences,
            int* h_numValuesPerSequence,
            int& totalNumValues
        ) const override {

            if(numSequences == 0) return;

            QueryData* const queryData = getQueryDataFromHandle(queryHandle);

            queryData->ranges.clear();

            totalNumValues = 0;

            queryData->ranges.resize(numSequences * numSmallest);
            CPUSequenceHasher<kmer_type> sequenceHasher;

            for(int s = 0; s < numSequences; s++){
                const int length = h_sequenceLengths[s];
                const unsigned int* sequence = h_sequenceData2Bit + encodedSequencePitchInInts * s;

                auto hashValues = sequenceHasher.getTopSmallestKmerHashes(
                    sequence, 
                    length, 
                    getKmerSize(), 
                    numSmallest
                );

                int numValues = 0;

                const auto numHashValues = hashValues.size();
                for(std::size_t i = 0; i < numHashValues; i++){
                    const auto mapQueryResult = kvtable->query(hashValues[i]);


                    queryData->ranges[s * numSmallest + i] 
                        = std::make_pair(mapQueryResult.valuesBegin, mapQueryResult.valuesBegin + mapQueryResult.numValues);

                    numValues += mapQueryResult.numValues;
                }

                for(int i = numHashValues; i < numSmallest; i++){
                    queryData->ranges[s * numSmallest + i] = std::make_pair(nullptr, nullptr);
                }

                h_numValuesPerSequence[s] = numValues;
                totalNumValues += numValues;
            }
           
            queryData->previousStage = QueryData::Stage::NumValues;
        }

        void retrieveValues(
            MinhasherHandle& queryHandle,
            const read_number* h_readIds,
            int numSequences,
            int /*totalNumValues*/,
            read_number* h_values,
            int* h_numValuesPerSequence,
            int* h_offsets //numSequences + 1
        ) const override {
            if(numSequences == 0) return;

            QueryData* const queryData = getQueryDataFromHandle(queryHandle);

            assert(queryData->previousStage == QueryData::Stage::NumValues);

            h_offsets[0] = 0;

            for(int s = 0; s < numSequences; s++){
                int numValues = 0;
                for(int i = 0; i < numSmallest; i++){
                    numValues += std::distance(
                        queryData->ranges[s * numSequences + i].first,
                        queryData->ranges[s * numSequences + i].second
                    );
                }

                std::vector<Value_t> valuestmp(numValues);
                auto valueIter = valuestmp.begin();
                for(int i = 0; i < numSmallest; i++){
                    valueIter = std::copy(
                        queryData->ranges[s * numSequences + i].first,
                        queryData->ranges[s * numSequences + i].second,
                        valueIter
                    );
                }

                std::sort(valuestmp.begin(), valuestmp.end());
                auto uniqueEnd = std::unique(valuestmp.begin(), valuestmp.end());

                if(h_readIds != nullptr){
                    auto readIdPos = std::lower_bound(
                        valuestmp.begin(),
                        uniqueEnd,
                        h_readIds[s]
                    );

                    if(readIdPos != uniqueEnd && *readIdPos == h_readIds[s]){
                        //TODO optimization: avoid this copy, instead skip the element when copying to h_values
                        uniqueEnd = std::copy(readIdPos + 1, uniqueEnd, readIdPos);
                    }
                }

                std::copy(valuestmp.begin(), uniqueEnd, h_values + h_offsets[s]);

                h_numValuesPerSequence[s] = std::distance(valuestmp.begin(), uniqueEnd);
                h_offsets[s+1] = h_offsets[s] + std::distance(valuestmp.begin(), uniqueEnd);
            }

            queryData->previousStage = QueryData::Stage::Retrieve;
        }

        MemoryUsage getMemoryInfo() const noexcept override{
            MemoryUsage result;

            result += kvtable->getMemoryInfo();

            return result;
        }

        MemoryUsage getMemoryInfo(const MinhasherHandle& handle) const noexcept override{
            return getQueryDataFromHandle(handle)->getMemoryInfo();
        }

        int getNumResultsPerMapThreshold() const noexcept override{
            return resultsPerMapThreshold;
        }
        
        int getNumberOfMaps() const noexcept override{
            return 1;
        }

        int getKmerSize() const noexcept override{
            return kmerSize;
        }

        void destroy() {
            
        }

        std::uint64_t getKmerMask() const{
            constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;

            return std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - getKmerSize()) * 2);
        }

        std::size_t computeHashesAndReadIds(
            kmer_type* keyoutput,
            read_number* readIdsOutput,
            const unsigned int* h_sequenceData2Bit,
            int numSequences,
            const int* h_sequenceLengths,
            std::size_t encodedSequencePitchInInts,
            const read_number* h_readIds
        ){
            if(numSequences == 0) return 0;

            ThreadPool::ParallelForHandle pforHandle{};

            ForLoopExecutor forLoopExecutor(threadPool, &pforHandle);
            const int numThreads = forLoopExecutor.getNumThreads();

            struct ThreadData{
                std::vector<kmer_type> hashes{};
                std::vector<read_number> ids{};
            };

            std::vector<ThreadData> threadData(numThreads);

            auto hashloopbody = [&](auto begin, auto end, int threadid){
                CPUSequenceHasher<kmer_type> sequenceHasher;

                threadData[threadid].hashes.resize((end - begin) * numSmallest);
                threadData[threadid].ids.resize((end - begin) * numSmallest);

                auto hashIter = threadData[threadid].hashes.begin();
                auto idIter = threadData[threadid].ids.begin();

                for(int s = begin; s < end; s++){
                    const int length = h_sequenceLengths[s];
                    const unsigned int* sequence = h_sequenceData2Bit + encodedSequencePitchInInts * s;

                    auto hashValues = sequenceHasher.getTopSmallestKmerHashes(
                        sequence, 
                        length, 
                        getKmerSize(), 
                        numSmallest
                    );

                    hashIter = std::copy(hashValues.begin(), hashValues.end(), hashIter);
                    std::fill(idIter, idIter + hashValues.size(), h_readIds[s]);
                    idIter = idIter + hashValues.size();
                }

                threadData[threadid].hashes.erase(hashIter, threadData[threadid].hashes.end());
                threadData[threadid].ids.erase(idIter, threadData[threadid].ids.end());

                assert(std::distance(threadData[threadid].hashes.begin(), hashIter) == std::distance(threadData[threadid].ids.begin(), idIter));
            };

            forLoopExecutor(0, numSequences, hashloopbody);

            std::size_t numKeys = 0;

            for(const auto& data : threadData){
                std::copy(data.hashes.begin(), data.hashes.end(), keyoutput + numKeys);
                std::copy(data.ids.begin(), data.ids.end(), readIdsOutput + numKeys);
                numKeys += data.hashes.size();
            }

            return numKeys;
        }

        void setThreadPool(ThreadPool* tp){
            threadPool = tp;
        }

        void setMemoryLimitForConstruction(std::size_t limit){
            memoryLimit = limit;
        }

        void writeToStream(std::ostream& os) const{

            os.write(reinterpret_cast<const char*>(&kmerSize), sizeof(int));
            os.write(reinterpret_cast<const char*>(&numSmallest), sizeof(int));
            os.write(reinterpret_cast<const char*>(&resultsPerMapThreshold), sizeof(int));

            os.write(reinterpret_cast<const char*>(&loadfactor), sizeof(float));

            kvtable->writeToStream(os);
        }

        int loadFromStream(std::ifstream& is, int /*numMapsUpperLimit = std::numeric_limits<int>::max()*/){
            destroy();

            is.read(reinterpret_cast<char*>(&kmerSize), sizeof(int));
            is.read(reinterpret_cast<char*>(&numSmallest), sizeof(int));
            is.read(reinterpret_cast<char*>(&resultsPerMapThreshold), sizeof(int));

            is.read(reinterpret_cast<char*>(&loadfactor), sizeof(float));

            kvtable = std::make_unique<DoublePassMultiValueHashTable<kmer_type, read_number>>(1, loadfactor);
            kvtable->loadFromStream(is);

            return 0;
        }

    private:

        QueryData* getQueryDataFromHandle(const MinhasherHandle& queryHandle) const{
            std::shared_lock<SharedMutex> lock(sharedmutex);

            return tempdataVector[queryHandle.getId()].get();
        }

        mutable int counter = 0;
        mutable SharedMutex sharedmutex{};

        float loadfactor = 0.8f;
        int numSmallest = 1;
        int maxNumKeys{};
        int kmerSize{};
        int resultsPerMapThreshold{};
        std::size_t memoryLimit;
        ThreadPool* threadPool{};
        mutable std::vector<std::unique_ptr<QueryData>> tempdataVector{};
        std::unique_ptr<DoublePassMultiValueHashTable<kmer_type, read_number>> kvtable{};
    };


}

#endif


#endif