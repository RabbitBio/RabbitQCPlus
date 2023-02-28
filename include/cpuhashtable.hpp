#ifndef CARE_CPUHASHTABLE_HPP
#define CARE_CPUHASHTABLE_HPP 

#include <config.hpp>
#include <memorymanagement.hpp>
#include <threadpool.hpp>
#include <util.hpp>
#include <hostdevicefunctions.cuh>

#include <map>
#include <vector>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <type_traits>
#include <limits>
#include <iostream>
#include <functional>

#include <hpc_helpers.cuh>


namespace care{

    //computes the new hashtable size from the current hashtable size on automatic rehash
    struct RehashPolicyDouble {
        constexpr std::size_t operator()(const std::size_t x) noexcept { return x * 2; }
    };

    template<class Key, class Value, class RehashPolicy = RehashPolicyDouble>
    class AoSCpuSingleValueHashTable{
        static_assert(std::is_integral<Key>::value, "Key must be integral!");

    public:
        class QueryResult{
        public:
            QueryResult() = default;
            QueryResult(const QueryResult&) = default;
            QueryResult(QueryResult&&) = default;
            QueryResult& operator=(const QueryResult&) = default;
            QueryResult& operator=(QueryResult&&) = default;

            QueryResult(bool b, Value v) : keyFound(b), val(std::move(v)){}

            bool valid() const{
                return keyFound;
            }
            const Value& value() const{
                return val;
            }
            Value& value(){
                return val;
            }
        private:
            bool keyFound;
            Value val;
        };

        AoSCpuSingleValueHashTable(const AoSCpuSingleValueHashTable&) = default;
        AoSCpuSingleValueHashTable(AoSCpuSingleValueHashTable&&) = default;
        AoSCpuSingleValueHashTable& operator=(const AoSCpuSingleValueHashTable&) = default;
        AoSCpuSingleValueHashTable& operator=(AoSCpuSingleValueHashTable&&) = default;

        AoSCpuSingleValueHashTable(std::size_t size = 1, float load = 0.8)
            : load(load), size(size), capacity(size/load)
        {
            storage.resize(capacity, emptySlot);
        }


        bool operator==(const AoSCpuSingleValueHashTable& rhs) const{
            return emptySlot == rhs.emptySlot 
                && feq(load, rhs.load)
                && numKeys == rhs.numKeys
                && maxProbes == rhs.maxProbes
                && size == rhs.size 
                && capacity == rhs.capacity
                && storage == rhs.storage;
        }

        bool operator!=(const AoSCpuSingleValueHashTable& rhs) const{
            return !(operator==(rhs));
        }

        void rehash(std::size_t newsize){
            newsize = std::max(newsize, std::size_t(1));

            std::cerr << "rehash(" << newsize << ")\n";
            AoSCpuSingleValueHashTable newtable(newsize, load);

            forEachKeyValuePair([&](const auto& key, const auto& value){
                newtable.insert(key, value);
            });

            std::swap(newtable, *this);
        }

        //if key already exists, current value is overwritten by passed value
        void insert(const Key& key, const Value& value){
            using hasher = hashers::MurmurHash<std::uint64_t>;

            if(numKeys >= size){
                rehash(RehashPolicy()(size));
            }

            const std::uint64_t key64 = std::uint64_t(key);
            std::size_t probes = 0;
            std::size_t pos = hasher::hash(key64) % capacity;
            while(storage[pos] != emptySlot && storage[pos].first != key){
                pos++;
                //wrap-around
                if(pos == capacity){
                    pos = 0;
                }
                probes++;
            }
            if(storage[pos].first == key){
                storage[pos].second = value;
            }else{
                //found empty slot
                numKeys++;
                storage[pos].first = key;
                storage[pos].second = value;
                maxProbes = std::max(maxProbes, probes);
            }
        }

        template<class Func>
        void insertOrUpdate(const Key& key, const Value& insertValue, Func&& updateFunc){
            using hasher = hashers::MurmurHash<std::uint64_t>;

            if(numKeys >= size){
                rehash(RehashPolicy()(size));
            }

            const std::uint64_t key64 = std::uint64_t(key);
            std::size_t probes = 0;
            std::size_t pos = hasher::hash(key64) % capacity;
            while(storage[pos] != emptySlot && storage[pos].first != key){
                pos++;
                //wrap-around
                if(pos == capacity){
                    pos = 0;
                }
                probes++;
            }
            if(storage[pos].first == key){
                updateFunc(storage[pos].second);
            }else{
                //found empty slot
                numKeys++;
                storage[pos].first = key;
                storage[pos].second = insertValue;
                maxProbes = std::max(maxProbes, probes);
            }
        }

        QueryResult query(const Key& key) const{
            using hasher = hashers::MurmurHash<std::uint64_t>;
            
            const std::uint64_t key64 = std::uint64_t(key);
            std::size_t probes = 0;
            std::size_t pos = hasher::hash(key64) % capacity;
            while(storage[pos].first != key){
                if(storage[pos] == emptySlot){
                    return {false, Value()};
                }
                pos++;
                //wrap-around
                if(pos == capacity){
                    pos = 0;
                }
                probes++;
                if(maxProbes < probes){
                    return {false, Value()};
                }
            }
            return {true, storage[pos].second};
        }

        Value* queryPointer(const Key& key){
            using hasher = hashers::MurmurHash<std::uint64_t>;
            
            const std::uint64_t key64 = std::uint64_t(key);
            std::size_t probes = 0;
            std::size_t pos = hasher::hash(key64) % capacity;
            while(storage[pos].first != key){
                if(storage[pos] == emptySlot){
                    return nullptr;
                }
                pos++;
                //wrap-around
                if(pos == capacity){
                    pos = 0;
                }
                probes++;
                if(maxProbes < probes){
                    return nullptr;
                }
            }
            return &storage[pos].second;
        }

        //Func(key, value)
        template<class Func>
        void forEachKeyValuePair(Func&& func){
            for(std::size_t i = 0; i < capacity; i++){
                if(storage[i] != emptySlot){
                    func(storage[i].first, storage[i].second);
                }
            }
        }

        MemoryUsage getMemoryInfo() const{
            MemoryUsage result;
            result.host = sizeof(Data) * capacity;

            return result;
        }

        void writeToStream(std::ostream& os) const{
            os.write(reinterpret_cast<const char*>(&load), sizeof(float));
            os.write(reinterpret_cast<const char*>(&numKeys), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(&maxProbes), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(&size), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(&capacity), sizeof(std::size_t));

            const std::size_t elements = storage.size();
            const std::size_t bytes = sizeof(Data) * elements;
            os.write(reinterpret_cast<const char*>(&elements), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(storage.data()), bytes);
        }

        void loadFromStream(std::ifstream& is){
            destroy();

            is.read(reinterpret_cast<char*>(&load), sizeof(float));
            is.read(reinterpret_cast<char*>(&numKeys), sizeof(std::size_t));
            is.read(reinterpret_cast<char*>(&maxProbes), sizeof(std::size_t));
            is.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));
            is.read(reinterpret_cast<char*>(&capacity), sizeof(std::size_t));

            std::size_t elements;
            is.read(reinterpret_cast<char*>(&elements), sizeof(std::size_t));
            storage.resize(elements);
            const std::size_t bytes = sizeof(Data) * elements;
            is.read(reinterpret_cast<char*>(storage.data()), bytes);
        }

        void destroy(){
            std::vector<Data> tmp;
            std::swap(storage, tmp);

            numKeys = 0;
            maxProbes = 0;
            size = 0;
            capacity = 0;
        }

        std::size_t getCapacity() const{
            return capacity;
        }

        std::map<int, int> getCountDistribution() const{
            std::map<int, int> map;

            for(std::size_t i = 0; i < capacity; i++){
                if(storage[i] != emptySlot){
                    map[storage[i].second.second]++;
                }
            }
            return map;
        }

        std::size_t getNumKeys() const noexcept{
            return numKeys;
        }

    private:

        using Data = std::pair<Key,Value>;

        Data emptySlot 
            = std::pair<Key,Value>{std::numeric_limits<Key>::max(), std::numeric_limits<Value>::max()};

        float load{};
        std::size_t numKeys{};
        std::size_t maxProbes{};
        std::size_t size{};
        std::size_t capacity{};
        std::vector<Data> storage{};        
    };


    template<class Key, class Value>
    class SoACpuSingleValueHashTable{
        static_assert(std::is_integral<Key>::value, "Key must be integral!");

    public:
        class QueryResult{
        public:
            QueryResult() = default;
            QueryResult(const QueryResult&) = default;
            QueryResult(QueryResult&&) = default;
            QueryResult& operator=(const QueryResult&) = default;
            QueryResult& operator=(QueryResult&&) = default;

            QueryResult(bool b, Value v) : keyFound(b), val(std::move(v)){}

            bool valid() const{
                return keyFound;
            }
            const Value& value() const{
                return val;
            }
            Value& value(){
                return val;
            }
        private:
            bool keyFound;
            Value val;
        };

        SoACpuSingleValueHashTable() = default;
        SoACpuSingleValueHashTable(const SoACpuSingleValueHashTable&) = default;
        SoACpuSingleValueHashTable(SoACpuSingleValueHashTable&&) = default;
        SoACpuSingleValueHashTable& operator=(const SoACpuSingleValueHashTable&) = default;
        SoACpuSingleValueHashTable& operator=(SoACpuSingleValueHashTable&&) = default;

        SoACpuSingleValueHashTable(std::size_t size, float load)
            : load(load), size(size), capacity(size/load), 
              keys(capacity, emptyKey), values(capacity){
        }

        SoACpuSingleValueHashTable(std::size_t size) 
            : SoACpuSingleValueHashTable(size, 0.8f){            
        }

        bool operator==(const SoACpuSingleValueHashTable& rhs) const{
            return emptyKey == rhs.emptyKey 
                && feq(load, rhs.load)
                && maxProbes == rhs.maxProbes
                && size == rhs.size 
                && capacity == rhs.capacity
                && keys == rhs.keys
                && values = rhs.values;
        }

        bool operator!=(const SoACpuSingleValueHashTable& rhs) const{
            return !(operator==(rhs));
        }

        void insert(const Key& key, const Value& value){
            if(key == emptyKey){
                std::cerr << "1\n";
            }
            using hasher = hashers::MurmurHash<std::uint64_t>;

            std::size_t probes = 0;
            std::size_t pos = hasher::hash(key) % capacity;
            while(keys[pos] != emptyKey){
                pos++;
                //wrap-around
                if(pos == capacity){
                    pos = 0;
                }
                probes++;
            }
            keys[pos] = key;
            values[pos] = value;

            maxProbes = std::max(maxProbes, probes);
        }

        QueryResult query(const Key& key) const{
            using hasher = hashers::MurmurHash<std::uint64_t>;
            
            std::size_t probes = 0;
            std::size_t pos = hasher::hash(key) % capacity;
            while(keys[pos] != key){
                if(keys[pos] == emptyKey){
                    return {false, Value{}};
                }
                pos++;
                //wrap-around
                if(pos == capacity){
                    pos = 0;
                }
                probes++;
                if(maxProbes < probes){
                    return {false, Value{}};
                }
            }
            return {true, values[pos]};
        }

        MemoryUsage getMemoryInfo() const{
            MemoryUsage result;
            result.host = sizeof(Key) * capacity;
            result.host += sizeof(Value) * capacity;

            return result;
        }

        void writeToStream(std::ostream& os) const{
            os.write(reinterpret_cast<const char*>(&load), sizeof(float));
            os.write(reinterpret_cast<const char*>(&maxProbes), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(&size), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(&capacity), sizeof(std::size_t));

            const std::size_t elements = keys.size();
            const std::size_t keysbytes = sizeof(Key) * elements;
            const std::size_t valuesbytes = sizeof(Value) * elements;
            os.write(reinterpret_cast<const char*>(&elements), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(keys.data()), keysbytes);
            os.write(reinterpret_cast<const char*>(values.data()), valuesbytes);
        }

        void loadFromStream(std::ifstream& is){
            destroy();

            is.read(reinterpret_cast<char*>(&load), sizeof(float));
            is.read(reinterpret_cast<char*>(&maxProbes), sizeof(std::size_t));
            is.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));
            is.read(reinterpret_cast<char*>(&capacity), sizeof(std::size_t));

            std::size_t elements;
            is.read(reinterpret_cast<char*>(&elements), sizeof(std::size_t));
            keys.resize(elements);
            const std::size_t keysbytes = sizeof(Key) * elements;
            is.read(reinterpret_cast<char*>(keys.data()), keysbytes);

            values.resize(elements);
            const std::size_t valuesbytes = sizeof(Value) * elements;
            is.read(reinterpret_cast<char*>(values.data()), valuesbytes);
        }

        void destroy(){
            std::vector<Key> ktmp;
            std::swap(keys, ktmp);

            std::vector<Value> vtmp;
            std::swap(values, vtmp);

            maxProbes = 0;
            size = 0;
            capacity = 0;
        }

    private:

        Key emptyKey = std::numeric_limits<Key>::max();

        float load{};
        std::size_t maxProbes{};
        std::size_t size{};
        std::size_t capacity{};
        std::vector<Key> keys{};
        std::vector<Value> values{};
    };






    template<class Key, class Value>
    class CpuReadOnlyMultiValueHashTable{
        static_assert(std::is_integral<Key>::value, "Key must be integral!");
    public:

        struct QueryResult{
            int numValues;
            const Value* valuesBegin;
        };

        CpuReadOnlyMultiValueHashTable() = default;
        CpuReadOnlyMultiValueHashTable(const CpuReadOnlyMultiValueHashTable&) = default;
        CpuReadOnlyMultiValueHashTable(CpuReadOnlyMultiValueHashTable&&) = default;
        CpuReadOnlyMultiValueHashTable& operator=(const CpuReadOnlyMultiValueHashTable&) = default;
        CpuReadOnlyMultiValueHashTable& operator=(CpuReadOnlyMultiValueHashTable&&) = default;

        CpuReadOnlyMultiValueHashTable(
            std::uint64_t maxNumValues_,
            float loadfactor_
        ) : loadfactor(loadfactor_), buildMaxNumValues{maxNumValues_}{
            buildkeys.reserve(buildMaxNumValues);
            buildvalues.reserve(buildMaxNumValues);
        }

        bool operator==(const CpuReadOnlyMultiValueHashTable& rhs) const{
            return values == rhs.values && lookup == rhs.lookup;
        }

        bool operator!=(const CpuReadOnlyMultiValueHashTable& rhs) const{
            return !(operator==(rhs));
        }

        void init(
            std::vector<Key> keys, 
            std::vector<Value> vals, 
            int maxValuesPerKey,
            ThreadPool* threadPool,
            bool valuesOfSameKeyMustBeSorted = true
        ){
            init(std::move(keys), std::move(vals), maxValuesPerKey, threadPool, {}, valuesOfSameKeyMustBeSorted);
        }

        template<class GroupByKeyOp>
        void init(
            GroupByKeyOp groupByKey,
            std::vector<Key> keys, 
            std::vector<Value> vals,
            ThreadPool* threadPool
        ){
            assert(keys.size() == vals.size());

            //std::cerr << "init valuesOfSameKeyMustBeSorted = " << valuesOfSameKeyMustBeSorted << "\n";

            if(isInit) return;

            std::vector<read_number> countsPrefixSum;
            values = std::move(vals);

            groupByKey(keys, values, countsPrefixSum);

            std::size_t nonEmtpyKeysCount = 0;
            for(std::size_t i = 0; i < keys.size(); i++){
                const auto count = countsPrefixSum[i+1] - countsPrefixSum[i];
                if(count > 0){
                    nonEmtpyKeysCount++;
                }
            }

            //lookup = std::move(AoSCpuSingleValueHashTable<Key, ValueIndex>(keys.size(), loadfactor));
            lookup = std::move(AoSCpuSingleValueHashTable<Key, ValueIndex>(nonEmtpyKeysCount, loadfactor));

            auto buildKeyLookup = [me=this, keys = std::move(keys), countsPrefixSum = std::move(countsPrefixSum)](){
                for(std::size_t i = 0; i < keys.size(); i++){
                    const auto count = countsPrefixSum[i+1] - countsPrefixSum[i];
                    if(count > 0){
                        me->lookup.insert(
                            keys[i], 
                            ValueIndex{countsPrefixSum[i], count}
                        );
                    }
                }
                me->isInit = true;
            };

            if(threadPool != nullptr){
                threadPool->enqueue(std::move(buildKeyLookup));
            }else{
                buildKeyLookup();
            }
        }

        void insert(const Key* keys, const Value* values, int N){
            assert(keys != nullptr);
            assert(values != nullptr);
            assert(buildMaxNumValues >= buildkeys.size() + N);

            buildkeys.insert(buildkeys.end(), keys, keys + N);
            buildvalues.insert(buildvalues.end(), values, values + N);
        }

        template<class GroupByKeyOp>
        void finalize(GroupByKeyOp groupByKey, ThreadPool* threadPool){
            init(groupByKey, std::move(buildkeys), std::move(buildvalues), threadPool);            
        }

        bool isInitialized() const noexcept{
            return isInit;
        }

        QueryResult query(const Key& key) const{
            assert(isInit);

            auto lookupQueryResult = lookup.query(key);

            if(lookupQueryResult.valid()){
                QueryResult result;

                result.numValues = lookupQueryResult.value().second;
                const auto valuepos = lookupQueryResult.value().first;
                result.valuesBegin = &values[valuepos];

                return result;
            }else{
                QueryResult result;

                result.numValues = 0;
                result.valuesBegin = nullptr;

                return result;
            }
        }

        void query(const Key* keys, std::size_t numKeys, QueryResult* resultsOutput) const{
            assert(isInit);
            for(std::size_t i = 0; i < numKeys; i++){
                resultsOutput[i] = query(keys[i]);
            }
        }

        MemoryUsage getMemoryInfo() const{
            MemoryUsage result;
            result.host = sizeof(Value) * values.capacity();
            result.host += lookup.getMemoryInfo().host;
            result.host += sizeof(Key) * buildkeys.capacity();
            result.host += sizeof(Value) * buildvalues.capacity();

            result.device = lookup.getMemoryInfo().device;

            //std::cerr << lookup.getMemoryInfo().host << " " << result.host << " bytes\n";

            // std::cerr << "readonlytable. keys capacity: " << lookup.getCapacity() << ", values.size() " << values.size() << "\n";

            // auto map = lookup.getCountDistribution();
            // for(auto pair : map){
            //     std::cerr << pair.first << ": " << pair.second << "\n";
            // }

            return result;
        }

        void writeToStream(std::ostream& os) const{
            assert(isInit);

            const std::size_t elements = values.size();
            const std::size_t bytes = sizeof(Value) * elements;
            os.write(reinterpret_cast<const char*>(&elements), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(values.data()), bytes);

            lookup.writeToStream(os);
        }

        void loadFromStream(std::ifstream& is){
            destroy();

            std::size_t elements;
            is.read(reinterpret_cast<char*>(&elements), sizeof(std::size_t));
            values.resize(elements);
            const std::size_t bytes = sizeof(Value) * elements;
            is.read(reinterpret_cast<char*>(values.data()), bytes);

            lookup.loadFromStream(is);
            isInit = true;
        }

        void destroy(){
            std::vector<Value> tmp;
            std::swap(values, tmp);

            lookup.destroy();
            isInit = false;
        }

        static std::size_t estimateGpuMemoryRequiredForInit(std::size_t numElements){

            std::size_t mem = 0;
            mem += sizeof(Key) * numElements; //d_keys
            mem += sizeof(Value) * numElements; //d_values
            mem += sizeof(read_number) * numElements; //d_indices
            mem += std::max(sizeof(read_number), sizeof(Value)) * numElements; //d_indices_tmp for sorting d_indices or d_values_tmp for sorted values

            return mem;
        }

    private:

        using ValueIndex = std::pair<read_number, BucketSize>;
        bool isInit = false;
        float loadfactor = 0.8f;
        std::uint64_t buildMaxNumValues = 0;
        std::vector<Key> buildkeys;
        std::vector<Value> buildvalues;
        // values with the same key are stored in contiguous memory locations
        // a single-value hashmap maps keys to the range of the corresponding values
        std::vector<Value> values; 
        AoSCpuSingleValueHashTable<Key, ValueIndex> lookup;
    };




    template<class Key, class Value>
    class DoublePassMultiValueHashTable{
        static_assert(std::is_integral<Key>::value, "Key must be integral!");
    public:

        struct QueryResult{
            int numValues;
            const Value* valuesBegin;
        };

        DoublePassMultiValueHashTable() = default;
        DoublePassMultiValueHashTable(const DoublePassMultiValueHashTable&) = default;
        DoublePassMultiValueHashTable(DoublePassMultiValueHashTable&&) = default;
        DoublePassMultiValueHashTable& operator=(const DoublePassMultiValueHashTable&) = default;
        DoublePassMultiValueHashTable& operator=(DoublePassMultiValueHashTable&&) = default;

        //dummy value cannot be inserted into hash table
        static constexpr Value dummyValue(){
            return std::numeric_limits<Value>::max();
        }

        DoublePassMultiValueHashTable(
            std::size_t estimatedNumKeys,
            float loadfactor_
        ) : loadfactor(loadfactor_), lookup(estimatedNumKeys, loadfactor){
        }

        bool operator==(const DoublePassMultiValueHashTable& rhs) const{
            return false;
        }

        bool operator!=(const DoublePassMultiValueHashTable& rhs) const{
            return !(operator==(rhs));
        }

        void firstPassInsert(const Key* keys, const Value* values, std::size_t N){
            auto increment = [](MetaData& meta){
                std::size_t currentCount = meta.firstPassCount();
                meta.firstPassCount(currentCount + 1);
            };

            MetaData count1Meta; 
            count1Meta.firstPassCount(1);

            for(std::size_t i = 0; i < N; i++){
                if(values[i] != dummyValue()){
                    lookup.insertOrUpdate(keys[i], count1Meta, increment);
                }
            }
        }

        void firstPassDone(std::size_t minValuesPerKey = 1, std::size_t maxValuesPerKey = std::numeric_limits<std::size_t>::max()){
            helpers::CpuTimer timer1("firstPassDone");

            maxValuesPerKey = std::min(maxValuesPerKey, (1ul << MetaData::numvaluesBits()) - 1);

            std::size_t numKeysWithOkValues = 0;
            std::size_t numValuesOfKeysWithOkValues = 0;

            lookup.forEachKeyValuePair(
                [&](const auto& /*key*/, auto& meta){
                    if(minValuesPerKey <= meta.firstPassCount() && meta.firstPassCount() <= maxValuesPerKey){
                        numKeysWithOkValues++;
                        numValuesOfKeysWithOkValues += meta.firstPassCount();
                    }else{
                        meta.firstPassCount(0);
                    }
                }
            );

            std::cerr << "unique keys: " << lookup.getNumKeys() << "\n";
            std::cerr << "numKeysWithOkValues: " << numKeysWithOkValues << "\n";
            std::cerr << "numValuesOfKeysWithOkValues: " << numValuesOfKeysWithOkValues << "\n";

            //try to create a new lookup which does not contain keys with not ok values
            try{
                helpers::CpuTimer timer2("shrink lookup");

                AoSCpuSingleValueHashTable<Key, MetaData> newLookup(lookup.getNumKeys(), loadfactor);

                lookup.forEachKeyValuePair(
                    [&](const auto& key, const auto& meta){                    
                        if(minValuesPerKey <= meta.firstPassCount() && meta.firstPassCount() <= maxValuesPerKey){
                            newLookup.insert(key, meta);
                        }
                    }
                );

                std::swap(lookup, newLookup);
                timer2.print();
            }catch(...){
                std::cerr << "Could not reduce memory usage of table\n";
            }

            valuestorage.resize(numValuesOfKeysWithOkValues, dummyValue());

            std::size_t position = 0;

            lookup.forEachKeyValuePair(
                [&](const auto& /*key*/, auto& meta){
                    std::size_t count = meta.firstPassCount();
                    meta.offset(position);
                    meta.numValues(count);
                    meta.tmpData(0);
                    position += count;
                }
            );

            timer1.print();
        }

        void secondPassInsert(const Key* keys, const Value* values, std::size_t N){
            for(std::size_t i = 0; i < N; i++){
                if(values[i] != dummyValue()){
                    MetaData* const meta = lookup.queryPointer(keys[i]);
                    if(meta != nullptr){

                        Value* const myValueRangeBegin = valuestorage.data() + meta->offset();
                        const std::size_t pos = meta->tmpData();
                        const std::size_t numValues = meta->numValues();
                        if(pos < numValues){
                            myValueRangeBegin[pos] = values[i];
                            meta->tmpData(pos+1);
                        }
                    }
                }
            }
        }

        void secondPassDone(){
            isInit = true;
        }

        QueryResult query(const Key& key) const{
            assert(isInit);

            auto lookupQueryResult = lookup.query(key);

            if(lookupQueryResult.valid()){
                QueryResult result;

                result.numValues = lookupQueryResult.value().numValues();
                const std::size_t valuepos = lookupQueryResult.value().offset();
                result.valuesBegin = &valuestorage[valuepos];

                return result;
            }else{
                QueryResult result;

                result.numValues = 0;
                result.valuesBegin = nullptr;

                return result;
            }
        }

        void query(const Key* keys, std::size_t numKeys, QueryResult* resultsOutput) const{
            assert(isInit);
            for(std::size_t i = 0; i < numKeys; i++){
                resultsOutput[i] = query(keys[i]);
            }
        }

        MemoryUsage getMemoryInfo() const{
            MemoryUsage result;

            result.host += valuestorage.capacity() * sizeof(Value);
            result += lookup.getMemoryInfo();            

            return result;
        }

        void writeToStream(std::ostream& os) const{
            os.write(reinterpret_cast<const char*>(&isInit), sizeof(bool));
            os.write(reinterpret_cast<const char*>(&loadfactor), sizeof(float));

            const std::size_t elements = valuestorage.size();
            const std::size_t bytes = sizeof(Value) * elements;
            os.write(reinterpret_cast<const char*>(&elements), sizeof(std::size_t));
            os.write(reinterpret_cast<const char*>(valuestorage.data()), bytes);

            lookup.writeToStream(os);
        }

        void loadFromStream(std::ifstream& is){
            destroy();

            is.read(reinterpret_cast<char*>(&isInit), sizeof(bool));
            is.read(reinterpret_cast<char*>(&loadfactor), sizeof(float));

            std::size_t elements;
            is.read(reinterpret_cast<char*>(&elements), sizeof(std::size_t));
            valuestorage.resize(elements);
            const std::size_t bytes = sizeof(Value) * elements;
            is.read(reinterpret_cast<char*>(valuestorage.data()), bytes);

            lookup.loadFromStream(is);
        }

        void destroy(){
            std::vector<Value> tmp;
            std::swap(valuestorage, tmp);

            lookup.destroy();
            isInit = false;
        }


    private:

        struct MetaData{
            static constexpr std::size_t offsetBits(){
                return 48;
            }
            static constexpr std::size_t numvaluesBits(){
                return 8;
            }

            static constexpr std::size_t numTmpBits(){
                return 8;
            }

            static_assert(numvaluesBits() == numTmpBits());

            // offset into value array, num values
            using Data = packed_types::PackedTriple<offsetBits(),numvaluesBits(),numTmpBits()>;

            static constexpr std::size_t maxFirstPassCount(){
                return (1ull << offsetBits()) - 1;
            }

            std::size_t firstPassCount() const noexcept{
                return data.first();
            }

            void firstPassCount(std::size_t newcount) noexcept{
                const std::size_t actualnewcount = std::min(maxFirstPassCount(), newcount);
                data.first(actualnewcount);
            }

            std::size_t numValues() const noexcept{
                return data.second();
            }

            void numValues(std::size_t value) noexcept{
                data.second(value);
            }

            std::size_t tmpData() const noexcept{
                return data.third();
            }

            void tmpData(std::size_t value) noexcept{
                data.third(value);
            }

            std::size_t offset() const noexcept{
                return data.first();
            }

            void offset(std::size_t value) noexcept{
                data.first(value);
            }

            bool operator==(const MetaData& rhs) const noexcept{
                return data == rhs.data;
            }

            bool operator!=(const MetaData& rhs) const noexcept{
                return !operator==(rhs);
            }

            Data data;
        };

        using ValueIndex = std::pair<read_number, BucketSize>;
        bool isInit = false;
        float loadfactor = 0.8f;
        // values with the same key are stored in contiguous memory locations
        // a single-value hashmap maps keys to the range of the corresponding values
        std::vector<Value> valuestorage; 
        AoSCpuSingleValueHashTable<Key, MetaData> lookup;
    };


}

#endif // CARE_CPUHASHTABLE_HPP
