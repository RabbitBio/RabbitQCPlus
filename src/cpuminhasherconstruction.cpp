
#include <cpuminhasherconstruction.hpp>
#include <cpureadstorage.hpp>
#include <cpuminhasher.hpp>
#include <ordinaryminhasher.hpp>
//#include <singlehashminhasher.hpp>

#include <minhasherlimit.hpp>

#include <options.hpp>

#include <memory>
#include <utility>


namespace care{

    std::string to_string(CpuMinhasherType type){
        switch(type){
            case CpuMinhasherType::Ordinary: return "Ordinary";
            //case CpuMinhasherType::OrdinarySingleHash: return "OrdinarySingleHash";
            case CpuMinhasherType::None: return "None";
            default: return "Unknown";
        }
    }

    void constructCpuMinhasherFromReadStorage(
        const ProgramOptions& programOptions,
        const CpuReadStorage& cpuReadStorage,
        CpuMinhasher* cpuMinhasher
    ){
        auto& readStorage = cpuReadStorage;

        const int requestedNumberOfMaps = programOptions.numHashFunctions;

        const read_number numReads = readStorage.getNumberOfReads();
        const int maximumSequenceLength = readStorage.getSequenceLengthUpperBound();
        const std::size_t encodedSequencePitchInInts = SequenceHelpers::getEncodedNumInts2Bit(maximumSequenceLength);
    
        constexpr read_number batchsize = 1000000;
        const int numBatches = SDIV(numReads, batchsize);

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

        const int hashFunctionOffset = 0;

        
        std::vector<int> usedHashFunctionNumbers;

        ThreadPool tpForHashing(programOptions.threads);
        ThreadPool tpForCompacting(std::min(2,programOptions.threads));   
        
        cpuMinhasher->setHostMemoryLimitForConstruction(maxMemoryForTables);
        cpuMinhasher->setDeviceMemoryLimitsForConstruction({0});

        std::vector<read_number> currentReadIds(batchsize);
        std::vector<unsigned int> sequencedata(batchsize * encodedSequencePitchInInts);
        std::vector<int> sequencelengths(batchsize);

        //std::size_t bytesOfCachedConstructedTables = 0;
        int remainingHashFunctions = requestedNumberOfMaps;
        bool keepGoing = true;

        while(remainingHashFunctions > 0 && keepGoing){

            cpuMinhasher->setThreadPool(&tpForHashing);

            const int alreadyExistingHashFunctions = requestedNumberOfMaps - remainingHashFunctions;
            std::vector<int> h_hashfunctionNumbers(remainingHashFunctions);
            std::iota(
                h_hashfunctionNumbers.begin(),
                h_hashfunctionNumbers.end(),
                alreadyExistingHashFunctions + hashFunctionOffset
            );

            int addedHashFunctions = cpuMinhasher->addHashTables(remainingHashFunctions,h_hashfunctionNumbers.data());

            if(addedHashFunctions == 0){
                keepGoing = false;
                break;
            }

            std::cout << "Constructing maps: ";
            for(int i = 0; i < addedHashFunctions; i++){
                std::cout << (alreadyExistingHashFunctions + i) << "(" << (hashFunctionOffset + alreadyExistingHashFunctions + i) << ") ";
            }
            std::cout << '\n';

            usedHashFunctionNumbers.insert(usedHashFunctionNumbers.end(), h_hashfunctionNumbers.begin(), h_hashfunctionNumbers.begin() + addedHashFunctions);

            for (int iter = 0; iter < numBatches; iter++){
                read_number readIdBegin = iter * batchsize;
                read_number readIdEnd = std::min((iter + 1) * batchsize, numReads);

                const std::size_t currentbatchsize = readIdEnd - readIdBegin;

                std::iota(currentReadIds.begin(), currentReadIds.end(), readIdBegin);

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

                cpuMinhasher->insert(
                    sequencedata.data(),
                    currentbatchsize,
                    sequencelengths.data(),
                    encodedSequencePitchInInts,
                    currentReadIds.data(),
                    alreadyExistingHashFunctions,
                    addedHashFunctions,
                    h_hashfunctionNumbers.data()
                );

                int errorcount = cpuMinhasher->checkInsertionErrors(
                    alreadyExistingHashFunctions,
                    addedHashFunctions
                );
                if(errorcount > 0){
                    throw std::runtime_error("An error occurred during hash table construction.");
                }   
            }

            std::cerr << "Compacting\n";
            if(tpForCompacting.getConcurrency() > 1){
                cpuMinhasher->setThreadPool(&tpForCompacting);
            }else{
                cpuMinhasher->setThreadPool(nullptr);
            }
            
            cpuMinhasher->compact();

            remainingHashFunctions -= addedHashFunctions;
        }

        cpuMinhasher->setThreadPool(nullptr);     
        cpuMinhasher->constructionIsFinished();
    }


    std::pair<std::unique_ptr<CpuMinhasher>, CpuMinhasherType>
    constructCpuMinhasherFromCpuReadStorage(
        const ProgramOptions& programOptions,
        const CpuReadStorage& cpuReadStorage,
        CpuMinhasherType requestedType
    ){
        std::unique_ptr<CpuMinhasher> cpuMinhasher;
        CpuMinhasherType cpuMinhasherType = CpuMinhasherType::None;

        auto makeOrdinary = [&](){                
            cpuMinhasher = std::make_unique<OrdinaryCpuMinhasher>(
                cpuReadStorage.getNumberOfReads(),
                calculateResultsPerMapThreshold(programOptions.estimatedCoverage),
                programOptions.kmerlength,
                programOptions.hashtableLoadfactor
            );

            cpuMinhasherType = CpuMinhasherType::Ordinary;
        };

        
        if(requestedType == CpuMinhasherType::Ordinary){
            makeOrdinary();
        }else{
            makeOrdinary();
        }

        if(programOptions.load_hashtables_from != "" && cpuMinhasher->canLoadFromStream()){

            std::ifstream is(programOptions.load_hashtables_from);
            assert((bool)is);

            const int loadedMaps = cpuMinhasher->loadFromStream(is, programOptions.numHashFunctions);

            std::cout << "Loaded " << loadedMaps << " hash tables from " << programOptions.load_hashtables_from << std::endl;
        }else{
            constructCpuMinhasherFromReadStorage(
                programOptions,
                cpuReadStorage,
                cpuMinhasher.get()
            );
        }

        return {std::move(cpuMinhasher), cpuMinhasherType};
    }
    
    
}
