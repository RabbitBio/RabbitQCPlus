#ifndef CARE_OrdinaryCpuMinhasher_HPP
#define CARE_OrdinaryCpuMinhasher_HPP


#include <cpuminhasher.hpp>
#include <cpureadstorage.hpp>
#include <groupbykey.hpp>


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

#include <cpusequencehasher.hpp>


#include <cassert>
#include <array>
#include <vector>
#include <memory>
#include <limits>
#include <string>
#include <fstream>
#include <algorithm>

namespace care{


    class OrdinaryCpuMinhasher : public CpuMinhasher{
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

        
    public:

        OrdinaryCpuMinhasher() : OrdinaryCpuMinhasher(0, 50, 16, 0.8f){

        }

        OrdinaryCpuMinhasher(int maxNumKeys_, int maxValuesPerKey, int k, float loadfactor_)
            : loadfactor(loadfactor_), maxNumKeys(maxNumKeys_), kmerSize(k), resultsPerMapThreshold(maxValuesPerKey){

        }

        OrdinaryCpuMinhasher(const OrdinaryCpuMinhasher&) = delete;
        OrdinaryCpuMinhasher(OrdinaryCpuMinhasher&&) = default;
        OrdinaryCpuMinhasher& operator=(const OrdinaryCpuMinhasher&) = delete;
        OrdinaryCpuMinhasher& operator=(OrdinaryCpuMinhasher&&) = default; 

        void setHostMemoryLimitForConstruction(std::size_t bytes) override{
            memoryLimit = bytes;
        }
        
        void setDeviceMemoryLimitsForConstruction(const std::vector<std::size_t>&) override{
            
        }

        void setThreadPool(ThreadPool* tp) override {
            threadPool = tp;
        }

        void constructionIsFinished() override {
            
        }

        int checkInsertionErrors(
            int /*firstHashfunction*/,
            int /*numHashfunctions*/
        ) override{
            return 0;
        }

        bool canWriteToStream() const noexcept override { return true; };
        bool canLoadFromStream() const noexcept override { return true; };


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

            std::vector<kmer_type> allHashValues(numSequences * getNumberOfMaps());

            CPUSequenceHasher<kmer_type> hasher;

            for(int s = 0; s < numSequences; s++){
                const int length = h_sequenceLengths[s];
                const unsigned int* const sequence = h_sequenceData2Bit + encodedSequencePitchInInts * s;

                hasher.hashInto(
                    allHashValues.begin() + getNumberOfMaps() * s,
                    sequence, 
                    length, 
                    getKmerSize(), 
                    getNumberOfMaps(),
                    0
                );
            }

            queryData->ranges.resize(numSequences * getNumberOfMaps());

            std::fill(h_numValuesPerSequence, h_numValuesPerSequence + numSequences, 0);

            for(int map = 0; map < getNumberOfMaps(); ++map){
                for(int s = 0; s < numSequences; s++){
                    const int length = h_sequenceLengths[s];

                    if(length >= getKmerSize()){
                        const kmer_type key = allHashValues[s * getNumberOfMaps() + map];
                        auto entries_range = queryMap(map, key);
                        const int n_entries = std::distance(entries_range.first, entries_range.second);
                        totalNumValues += n_entries;
                        h_numValuesPerSequence[s] += n_entries;
                        //queryData->ranges[s * getNumberOfMaps() + map] = entries_range;
                        queryData->ranges[map * numSequences + s] = entries_range;
                    }
                }
            }

            queryData->previousStage = QueryData::Stage::NumValues;
        }


        void retrieveValues(
            MinhasherHandle& queryHandle,
            int numSequences,
            int /*totalNumValues*/,
            read_number* h_values,
            const int* /*h_numValuesPerSequence*/,
            int* h_offsets //numSequences + 1
        ) const override {
            if(numSequences == 0) return;

            QueryData* const queryData = getQueryDataFromHandle(queryHandle);

            assert(queryData->previousStage == QueryData::Stage::NumValues);

            h_offsets[0] = 0;

            auto iter = h_values;
            for(int s = 0; s < numSequences; s++){
                for(int map = 0; map < getNumberOfMaps(); ++map){
                    const auto& range = queryData->ranges[map * numSequences + s];
                    iter = std::copy(range.first, range.second, iter);
                }
                h_offsets[s+1] = std::distance(h_values, iter);
            }

            queryData->previousStage = QueryData::Stage::Retrieve;
        }

        void compact() override {
            const int num = minhashTables.size();

            auto groupByKey = [&](auto& keys, auto& values, auto& countsPrefixSum){
                constexpr bool valuesOfSameKeyMustBeSorted = true;
                const int maxValuesPerKey = getNumResultsPerMapThreshold();

                //if only 1 value exists, it belongs to the anchor read itself and does not need to be stored.
                constexpr int minValuesPerKey = MINHASHER_MIN_VALUES_PER_KEY;

                care::GroupByKeyCpu<Key_t, Value_t, read_number> groupByKey(valuesOfSameKeyMustBeSorted, maxValuesPerKey, minValuesPerKey);
                groupByKey.execute(keys, values, countsPrefixSum);
            };

            for(int i = 0, l = 0; i < num; i++){
                auto& ptr = minhashTables[i];
            
                if(!ptr->isInitialized()){
                    //after processing 3 tables, available memory should be sufficient for multithreading
                    if(l >= 3){
                        ptr->finalize(groupByKey, threadPool);
                    }else{
                        ptr->finalize(groupByKey, nullptr);
                    }
                    l++;
                }                
            }

            if(threadPool != nullptr){
                threadPool->wait();
            }
        }

        MemoryUsage getMemoryInfo() const noexcept override{
            MemoryUsage result;

            result.host = sizeof(HashTable) * minhashTables.size();
            
            for(const auto& tableptr : minhashTables){
                auto m = tableptr->getMemoryInfo();
                result.host += m.host;

                for(auto pair : m.device){
                    result.device[pair.first] += pair.second;
                }
            }

            return result;
        }

        MemoryUsage getMemoryInfo(const MinhasherHandle& handle) const noexcept override{
            return getQueryDataFromHandle(handle)->getMemoryInfo();
        }

        int getNumResultsPerMapThreshold() const noexcept override{
            return resultsPerMapThreshold;
        }
        
        int getNumberOfMaps() const noexcept override{
            return minhashTables.size();
        }

        int getKmerSize() const noexcept override{
            return kmerSize;
        }

        void destroy() {
            minhashTables.clear();
        }

        void finalize(){
            compact();
        }

        std::uint64_t getKmerMask() const{
            constexpr int maximum_kmer_length = max_k<std::uint64_t>::value;

            return std::numeric_limits<std::uint64_t>::max() >> ((maximum_kmer_length - getKmerSize()) * 2);
        }

        void writeToStream(std::ostream& os) const override{

            os.write(reinterpret_cast<const char*>(&kmerSize), sizeof(int));
            os.write(reinterpret_cast<const char*>(&resultsPerMapThreshold), sizeof(int));

            os.write(reinterpret_cast<const char*>(&loadfactor), sizeof(float));

            const int numTables = getNumberOfMaps();
            os.write(reinterpret_cast<const char*>(&numTables), sizeof(int));

            for(const auto& tableptr : minhashTables){
                tableptr->writeToStream(os);
            }
        }

        int loadFromStream(std::ifstream& is, int numMapsUpperLimit = std::numeric_limits<int>::max()) override{
            destroy();

            is.read(reinterpret_cast<char*>(&kmerSize), sizeof(int));
            is.read(reinterpret_cast<char*>(&resultsPerMapThreshold), sizeof(int));

            is.read(reinterpret_cast<char*>(&loadfactor), sizeof(float));

            int numMaps = 0;

            is.read(reinterpret_cast<char*>(&numMaps), sizeof(int));

            const int mapsToLoad = std::min(numMapsUpperLimit, numMaps);

            for(int i = 0; i < mapsToLoad; i++){
                auto ptr = std::make_unique<HashTable>();
                ptr->loadFromStream(is);
                minhashTables.emplace_back(std::move(ptr));
            }

            return mapsToLoad;
        }
        

        int addHashTables(int numAdditionalTables, const int* /*hashFunctionIds*/) override{
            int added = 0;

            #ifndef NDEBUG
            const int cur = minhashTables.size();
            assert(!(numAdditionalTables + cur > 64));
            #endif

            std::size_t bytesOfCachedConstructedTables = 0;
            for(const auto& ptr : minhashTables){
                auto memusage = ptr->getMemoryInfo();
                bytesOfCachedConstructedTables += memusage.host;
            }

            std::size_t requiredMemPerTable = (sizeof(kmer_type) + sizeof(read_number)) * maxNumKeys;
            int numTablesToConstruct = (memoryLimit - bytesOfCachedConstructedTables) / requiredMemPerTable;
            numTablesToConstruct -= 2; // keep free memory of 2 tables to perform transformation 
            numTablesToConstruct = std::min(numTablesToConstruct, numAdditionalTables);
            //maxNumTablesInIteration = std::min(numTablesToConstruct, 4);

            for(int i = 0; i < numTablesToConstruct; i++){
                try{
                    auto ptr = std::make_unique<HashTable>(maxNumKeys, loadfactor);

                    minhashTables.emplace_back(std::move(ptr));
                    added++;
                }catch(...){

                }
            }

            return added;
        } 

        void insert(            
            const unsigned int* h_sequenceData2Bit,
            int numSequences,
            const int* h_sequenceLengths,
            std::size_t encodedSequencePitchInInts,
            const read_number* h_readIds,
            int firstHashfunction,
            int numHashfunctions,
            const int* /*h_hashFunctionNumbers*/
        ) override {
            if(numSequences == 0) return;

            ThreadPool::ParallelForHandle pforHandle{};

            ForLoopExecutor forLoopExecutor(threadPool, &pforHandle);

            std::vector<kmer_type> allHashValues(numSequences * getNumberOfMaps());

            auto hashloopbody = [&](auto begin, auto end, int /*threadid*/){
                CPUSequenceHasher<kmer_type> hasher;

                for(int s = begin; s < end; s++){
                    const int length = h_sequenceLengths[s];
                    const unsigned int* sequence = h_sequenceData2Bit + encodedSequencePitchInInts * s;

                    auto hashValues = hasher.hash(
                        sequence, 
                        length, 
                        getKmerSize(), 
                        getNumberOfMaps(),
                        0
                    );

                    for(int h = 0; h < getNumberOfMaps(); h++){
                        allHashValues[h * numSequences + s] = hashValues[h];
                    }
                }
            };

            forLoopExecutor(0, numSequences, hashloopbody);

            auto insertloopbody = [&](auto begin, auto end, int /*threadid*/){
                for(int h = begin; h < end; h++){
                    minhashTables[h]->insert(
                        &allHashValues[h * numSequences], h_readIds, numSequences
                    );
                }
            };

            forLoopExecutor(firstHashfunction, firstHashfunction + numHashfunctions, insertloopbody);
        }   


    private:

        QueryData* getQueryDataFromHandle(const MinhasherHandle& queryHandle) const{
            std::shared_lock<SharedMutex> lock(sharedmutex);

            return tempdataVector[queryHandle.getId()].get();
        }

        Range_t queryMap(int mapid, const Key_t& key) const{
            assert(mapid < getNumberOfMaps());

            const int numResultsPerMapQueryThreshold = getNumResultsPerMapThreshold();

            const auto mapQueryResult = minhashTables[mapid]->query(key);

            if(mapQueryResult.numValues > numResultsPerMapQueryThreshold){
                return std::make_pair(nullptr, nullptr); //return empty range
            }

            return std::make_pair(mapQueryResult.valuesBegin, mapQueryResult.valuesBegin + mapQueryResult.numValues);
        }


        mutable int counter = 0;
        mutable SharedMutex sharedmutex{};

        float loadfactor = 0.8f;
        int maxNumKeys{};
        int kmerSize{};
        int resultsPerMapThreshold{};
        ThreadPool* threadPool;
        std::size_t memoryLimit;
        std::vector<std::unique_ptr<HashTable>> minhashTables{};
        mutable std::vector<std::unique_ptr<QueryData>> tempdataVector{};
    };


}

#endif