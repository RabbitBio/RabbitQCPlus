#ifndef CARE_READSTORAGECONSTRUCTION2_HPP 
#define CARE_READSTORAGECONSTRUCTION2_HPP

#include <util.hpp>
#include <config.hpp>
#include <threadpool.hpp>
#include <readlibraryio.hpp>
#include <sequencehelpers.hpp>
#include <concurrencyhelpers.hpp>
#include <lengthstorage.hpp>
#include <chunkedreadstorage.hpp>
#include <options.hpp>
#include <qualityscorecompression.hpp>

#include <vector>
#include <array>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <string>
#include <iostream>
#include <future>
#include <limits>
#include <map>
#include <set>

namespace care{

std::unique_ptr<ChunkedReadStorage> constructChunkedReadStorageFromFiles(
    const ProgramOptions& programOptions
){
    const bool useQualityScores = programOptions.useQualityScores;
    const int numQualityBits = programOptions.qualityScoreBits;

    if(programOptions.load_binary_reads_from != ""){
        const bool pairedEnd = programOptions.pairType == SequencePairType::PairedEnd;
        std::unique_ptr<ChunkedReadStorage> readStorage = std::make_unique<ChunkedReadStorage>(pairedEnd, useQualityScores, numQualityBits);

        readStorage->loadFromFile(programOptions.load_binary_reads_from);
        std::cerr << "Loaded readstorage from file " << programOptions.load_binary_reads_from << "\n";

        return readStorage;
    }else{
        const bool pairedEnd = programOptions.pairType == SequencePairType::PairedEnd;
        std::unique_ptr<ChunkedReadStorage> readStorage = std::make_unique<ChunkedReadStorage>(pairedEnd, useQualityScores, numQualityBits);

        const bool showProgress = programOptions.showProgress;

        auto showProgressFunc = [showProgress](auto totalCount, auto seconds){
            if(showProgress){
                std::cout << "Processed " << totalCount << " reads in file. Elapsed time: " 
                                << seconds << " seconds." << std::endl;
            }
        };

        auto updateShowProgressInterval = [](auto duration){
            return duration * 2;
        };

        ProgressThread<std::size_t> progressThread(
            std::numeric_limits<std::size_t>::max(), 
            showProgressFunc, 
            updateShowProgressInterval
        );
                

        auto preprocessSequence = [&](std::string& sequence, int& Ncount){

            auto isValidBase = [](char c){
                constexpr std::array<char, 10> validBases{'A','C','G','T','a','c','g','t'};
                return validBases.end() != std::find(validBases.begin(), validBases.end(), c);
            };

            const int numundeterminedBasesInSequence = std::count_if(sequence.begin(), sequence.end(), [&](char c){
                return !isValidBase(c);
            });

            constexpr std::array<char, 4> bases = {'A', 'C', 'G', 'T'};

            for(auto& c : sequence){
                if(c == 'a') c = 'A';
                else if(c == 'c') c = 'C';
                else if(c == 'g') c = 'G';
                else if(c == 't') c = 'T';
                else if(!isValidBase(c)){
                    c = bases[Ncount];
                    Ncount = (Ncount + 1) % 4;
                }
            }

            return numundeterminedBasesInSequence > 0;
        };


        struct BatchFromFile{
            int validItems = 0;
            read_number firstReadId = 0;
            std::vector<std::string> sequences{};
            std::vector<std::string> qualities{};
        };

    
        SimpleConcurrentQueue<BatchFromFile*> freeBatchFromFile;
        SimpleConcurrentQueue<BatchFromFile*> unprocessedBatchFromFile;

        constexpr int fileParserMaxBatchsize = 65536;

        auto fileParserThreadFunction = [&](
            const std::string& filename, 
            std::size_t readIdOffset
        ){
            int batchId = 0;

            BatchFromFile* sbatch = nullptr;

            auto initbatch = [&](){
                sbatch = freeBatchFromFile.pop();
                sbatch->validItems = 0;
                sbatch->firstReadId = readIdOffset;
                sbatch->sequences.resize(fileParserMaxBatchsize);
                if(useQualityScores){
                    sbatch->qualities.resize(fileParserMaxBatchsize);
                }
            };

            initbatch();

            batchId++;

            std::size_t totalNumberOfReads = 0;

            forEachReadInFile(
                filename,
                [&](auto /*readnum*/, auto& read){

                    std::swap(sbatch->sequences[sbatch->validItems], read.sequence);
                    if(useQualityScores){
                        std::swap(sbatch->qualities[sbatch->validItems], read.quality);
                    }
                    sbatch->validItems++;

                    if(sbatch->validItems >= fileParserMaxBatchsize){
                        unprocessedBatchFromFile.push(sbatch);

                        readIdOffset += sbatch->validItems;

                        initbatch();

                        batchId++;
                    }

                    totalNumberOfReads++;
                }
            );        

            sbatch->sequences.resize(sbatch->validItems);
            sbatch->qualities.resize(sbatch->validItems);
            unprocessedBatchFromFile.push(sbatch);

            return totalNumberOfReads;
        };

        auto pairedfileParserThreadFunction = [&](
            const std::string& filename1, 
            const std::string& filename2,
            std::size_t readIdOffset
        ){
            int batchId = 0;

            BatchFromFile* sbatch = nullptr;

            auto initbatch = [&](){
                sbatch = freeBatchFromFile.pop();
                sbatch->validItems = 0;
                sbatch->firstReadId = readIdOffset;
                sbatch->sequences.resize(fileParserMaxBatchsize);
                if(useQualityScores){
                    sbatch->qualities.resize(fileParserMaxBatchsize);
                }
            };

            initbatch();

            batchId++;

            std::size_t totalNumberOfReads = 0;

            forEachReadInPairedFiles(
                filename1, filename2,
                [&](auto /*readnum*/, auto& read){

                    std::swap(sbatch->sequences[sbatch->validItems], read.sequence);
                    if(useQualityScores){
                        std::swap(sbatch->qualities[sbatch->validItems], read.quality);
                    }
                    sbatch->validItems++;

                    if(sbatch->validItems >= fileParserMaxBatchsize){
                        unprocessedBatchFromFile.push(sbatch);

                        readIdOffset += sbatch->validItems;

                        initbatch();

                        batchId++;
                    }

                    totalNumberOfReads++;
                }
            );        

            sbatch->sequences.resize(sbatch->validItems);
            sbatch->qualities.resize(sbatch->validItems);
            unprocessedBatchFromFile.push(sbatch);

            return totalNumberOfReads;
        };


        struct EncodedBatch{
            int validItems = 0;
            read_number firstReadId = 0;
            std::size_t encodedSequencePitchInInts = 0;
            //std::size_t qualityPitchInBytes = 0;
            std::size_t encodedQualityPitchInInts = 0;
            std::vector<int> sequenceLengths{};
            std::vector<unsigned int> encodedSequences{};
            //std::vector<char> qualities{};
            std::vector<unsigned int> encodedQualities{};
            std::vector<read_number> ambiguousReadIds{};
        };

        SimpleConcurrentQueue<EncodedBatch*> freeEncodedBatches;
        SimpleConcurrentQueue<EncodedBatch*> unprocessedEncodedBatches;

        auto encoderThreadFunction = [&](){
            BatchFromFile* sbatch = unprocessedBatchFromFile.pop();
            EncodedBatch* encbatch = nullptr;

            QualityCompressorWrapper qualityCompressor(numQualityBits);

            auto initEncBatch = [&](auto sequencepitchInInts, auto qualityPitchInInts){
                encbatch = freeEncodedBatches.pop();
                assert(encbatch != nullptr);

                encbatch->validItems = sbatch->validItems;
                encbatch->firstReadId = sbatch->firstReadId;
                encbatch->encodedSequencePitchInInts = sequencepitchInInts;
                encbatch->sequenceLengths.resize(encbatch->validItems);
                encbatch->encodedSequences.resize(encbatch->validItems * sequencepitchInInts);
                if(useQualityScores){
                    encbatch->encodedQualityPitchInInts = qualityPitchInInts;
                    encbatch->encodedQualities.resize(encbatch->validItems * qualityPitchInInts);
                }
                encbatch->ambiguousReadIds.clear();
            };

            while(sbatch != nullptr){
                int maxLength = 0;
                int Ncount = 0;

                for(int i = 0; i < sbatch->validItems; i++){
                    const int length = sbatch->sequences[i].length();
                    maxLength = std::max(maxLength, length);
                }

                const std::size_t sequencepitchInInts = SequenceHelpers::getEncodedNumInts2Bit(maxLength);
                const std::size_t qualityPitchInInts = QualityCompressionHelper::getNumInts(maxLength, numQualityBits);

                initEncBatch(sequencepitchInInts, qualityPitchInInts);

                for(int i = 0; i < sbatch->validItems; i++){
                    const int length = sbatch->sequences[i].length();
                    encbatch->sequenceLengths[i] = length;

                    bool isAmbig = preprocessSequence(sbatch->sequences[i], Ncount);
                    if(isAmbig){
                        const read_number readId = sbatch->firstReadId + i;
                        encbatch->ambiguousReadIds.emplace_back(readId);
                    }

                    SequenceHelpers::encodeSequence2Bit(
                        encbatch->encodedSequences.data() + i * sequencepitchInInts,
                        sbatch->sequences[i].c_str(),
                        length
                    );

                    if(useQualityScores){
                        qualityCompressor.encodeQualityString(
                            encbatch->encodedQualities.data() + i * qualityPitchInInts,
                            sbatch->qualities[i].data(),
                            length
                        );
                    }
                    
                    encbatch->sequenceLengths[i] = length;
                }

                freeBatchFromFile.push(sbatch);
                unprocessedEncodedBatches.push(encbatch);           

                sbatch = unprocessedBatchFromFile.pop();
            }
        };

        auto inserterThreadFunction = [&](){
            EncodedBatch* sbatch = unprocessedEncodedBatches.pop();

            while(sbatch != nullptr){
                const int numSequences = sbatch->validItems;

                if(sbatch->ambiguousReadIds.size() > 0){
                    readStorage->appendAmbiguousReadIds(
                        sbatch->ambiguousReadIds
                    );
                }

                readStorage->appendConsecutiveReads(
                    sbatch->firstReadId,
                    sbatch->validItems,
                    std::move(sbatch->sequenceLengths),
                    std::move(sbatch->encodedSequences),
                    sbatch->encodedSequencePitchInInts,
                    std::move(sbatch->encodedQualities),
                    sbatch->encodedQualityPitchInInts
                );

                progressThread.addProgress(numSequences);

                freeEncodedBatches.push(sbatch);
                sbatch = unprocessedEncodedBatches.pop();
            }
                    
        };



        constexpr int numParsers = 1;
        int numEncoders = 4;
        int numInserters = 1;

        const int numFilebatches = numParsers + numEncoders;

        std::vector<BatchFromFile> batchesFromFile(numFilebatches);

        for(int i = 0; i < numFilebatches; i++){
            freeBatchFromFile.push(&batchesFromFile[i]);
        }

        const int numEncodedBatches = numEncoders + numInserters;

        std::vector<EncodedBatch> encodedBatches(numEncodedBatches);

        for(int i = 0; i < numEncodedBatches; i++){
            freeEncodedBatches.push(&encodedBatches[i]);
        }

        std::vector<std::future<void>> encoderFutures;
        std::vector<std::future<void>> inserterFutures;
        
        for(int i = 0; i < numEncoders; i++){
            encoderFutures.emplace_back(
                std::async(
                    std::launch::async,
                    encoderThreadFunction
                )
            );
        }

        for(int i = 0; i < numInserters; i++){
            inserterFutures.emplace_back(
                std::async(
                    std::launch::async,
                    inserterThreadFunction
                )
            );
        }

        std::size_t totalNumReads = 0;

        if(programOptions.pairType == SequencePairType::SingleEnd){

            const int numInputFiles = programOptions.inputfiles.size();

            for(int i = 0; i < numInputFiles; i++){
                std::string inputfile = programOptions.inputfiles[i];

                std::future<std::size_t> future = std::async(
                    std::launch::async,
                    fileParserThreadFunction,
                    std::move(inputfile), totalNumReads
                );

                std::size_t numReads = future.get();
                totalNumReads += numReads;
            }

        }else if(programOptions.pairType == SequencePairType::PairedEnd){
            //paired end may be one of the following. 1 file with interleaved reads,
            //or two files with separate reads

            const int numInputFiles = programOptions.inputfiles.size();

            assert(numInputFiles > 0);
            assert(numInputFiles <= 2);

            if(numInputFiles == 1){

                std::string inputfile = programOptions.inputfiles[0];

                std::cout << "Converting paired reads of file " 
                    << inputfile  << "\n";

                std::future<std::size_t> future = std::async(
                    std::launch::async,
                    fileParserThreadFunction,
                    std::move(inputfile), totalNumReads
                );

                std::size_t numReads = future.get();
                totalNumReads += numReads;
            
            }else{
                assert(numInputFiles == 2);

                std::string filename1 = programOptions.inputfiles[0];
                std::string filename2 = programOptions.inputfiles[1];

                std::cout << "Converting paired reads of files " 
                    << filename1 << " and " << filename2 << ", storing them in memory\n";

                std::future<std::size_t> future = std::async(
                    std::launch::async,
                    pairedfileParserThreadFunction,
                    std::move(filename1), std::move(filename2), totalNumReads
                );

                std::size_t numReads = future.get();
                totalNumReads += numReads;
            
            }
        }else{
            throw std::runtime_error("invalid pair mode / number of input files");
        }

        //parsing done. flush queues to exit other threads

        for(int i = 0; i < numEncoders; i++){
            unprocessedBatchFromFile.push(nullptr);      
        }

        for(auto& f : encoderFutures){
            f.wait();
        }

        for(int i = 0; i < numInserters; i++){
            unprocessedEncodedBatches.push(nullptr);
        }
        
        for(auto& f : inserterFutures){
            f.wait();
        }
        
        progressThread.finished();
        if(showProgress){
            std::cout << "\n";
        }

        helpers::CpuTimer footimer("init readstorage after construction");
        
        readStorage->appendingFinished(
            programOptions.memoryTotalLimit - readStorage->getMemoryInfo().host
        );

        footimer.print();

        return readStorage;
    }
    
}

} //namespace care


#endif
