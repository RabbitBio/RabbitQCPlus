

#include <config.hpp>

#include <options.hpp>
#include <cpureadstorage.hpp>
#include <cpu_alignment.hpp>
#include <alignmentorientation.hpp>
#include <msa.hpp>
#include <qualityscoreweights.hpp>
#include <rangegenerator.hpp>
#include <correctedsequence.hpp>

#include <threadpool.hpp>
#include <serializedobjectstorage.hpp>
#include <util.hpp>
#include <filehelpers.hpp>
#include <hostdevicefunctions.cuh>

#include <classification.hpp>


#define ENABLE_CPU_CORRECTOR_TIMING
#include <corrector.hpp>
#include <cpuminhasher.hpp>

#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>

#include <vector>

#include <omp.h>

#include <thrust/iterator/counting_iterator.h>


namespace care{
namespace cpu{


SerializedObjectStorage correct_cpu(
    const ProgramOptions& programOptions,
    CpuMinhasher& minhasher,
    CpuReadStorage& readStorage
){

    omp_set_num_threads(programOptions.threads);

    const auto rsMemInfo = readStorage.getMemoryInfo();
    const auto mhMemInfo = minhasher.getMemoryInfo();

    std::size_t memoryAvailableBytesHost = programOptions.memoryTotalLimit;
    if(memoryAvailableBytesHost > rsMemInfo.host){
        memoryAvailableBytesHost -= rsMemInfo.host;
    }else{
        memoryAvailableBytesHost = 0;
    }
    if(memoryAvailableBytesHost > mhMemInfo.host){
        memoryAvailableBytesHost -= mhMemInfo.host;
    }else{
        memoryAvailableBytesHost = 0;
    }


    ReadCorrectionFlags correctionFlags(readStorage.getNumberOfReads());


    std::cerr << "correctionStatusFlagsPerRead bytes: " << correctionFlags.sizeInBytes() / 1024. / 1024. << " MB\n";

    if(memoryAvailableBytesHost > correctionFlags.sizeInBytes()){
        memoryAvailableBytesHost -= correctionFlags.sizeInBytes();
    }else{
        memoryAvailableBytesHost = 0;
    }

    const std::size_t availableMemoryInBytes = memoryAvailableBytesHost; //getAvailableMemoryInKB() * 1024;
    std::size_t memoryForPartialResultsInBytes = 0;

    if(availableMemoryInBytes > 2*(std::size_t(1) << 30)){
        memoryForPartialResultsInBytes = availableMemoryInBytes - 2*(std::size_t(1) << 30);
    }

    std::cerr << "Partial results may occupy " << (memoryForPartialResultsInBytes /1024. / 1024. / 1024.) 
        << " GB in memory. Remaining partial results will be stored in temp directory. \n";

    const std::size_t memoryLimitData = memoryForPartialResultsInBytes * 0.75;
    const std::size_t memoryLimitOffsets = memoryForPartialResultsInBytes * 0.25;

    SerializedObjectStorage partialResults(memoryLimitData, memoryLimitOffsets, programOptions.tempdirectory + "/");

    const std::size_t numReadsToProcess = getNumReadsToProcess(&readStorage, programOptions);

    IteratorRangeTraversal<thrust::counting_iterator<read_number>> readIdGenerator(
        thrust::make_counting_iterator<read_number>(0),
        thrust::make_counting_iterator<read_number>(0) + numReadsToProcess
    );
    
    BackgroundThread outputThread(true);
    outputThread.setMaximumQueueSize(programOptions.threads);

    CpuErrorCorrector::TimeMeasurements timingsOfAllThreads;

    ClfAgent clfAgent_(programOptions);
    
    auto showProgress = [&](auto totalCount, auto seconds){
        if(programOptions.showProgress){
            std::size_t totalNumReads = numReadsToProcess;

            printf("Processed %10u of %10lu reads (Runtime: %03d:%02d:%02d)\r",
                    totalCount, totalNumReads,
                    int(seconds / 3600),
                    int(seconds / 60) % 60,
                    int(seconds) % 60);
            std::cout.flush();

            if(totalCount == totalNumReads){
                std::cout << '\n';
            }
        }        
    };

    auto updateShowProgressInterval = [](auto duration){
        return duration;
    };

    ProgressThread<read_number> progressThread(numReadsToProcess, showProgress, updateShowProgressInterval);

    #pragma omp parallel
    {
        //const int threadId = omp_get_thread_num();

        ClfAgent clfAgent = clfAgent_;
        const std::size_t encodedSequencePitchInInts2Bit = SequenceHelpers::getEncodedNumInts2Bit(readStorage.getSequenceLengthUpperBound());
        const std::size_t decodedSequencePitchInBytes = readStorage.getSequenceLengthUpperBound();
        const std::size_t qualityPitchInBytes = readStorage.getSequenceLengthUpperBound();

        const int myBatchsize = std::max((readStorage.isPairedEnd() ? 2 : 1), programOptions.batchsize);

        CpuErrorCorrector errorCorrector(
            encodedSequencePitchInInts2Bit,
            decodedSequencePitchInBytes,
            qualityPitchInBytes,
            programOptions,
            minhasher,
            readStorage,
            correctionFlags,
            clfAgent
        );

        std::vector<read_number> batchReadIds(myBatchsize);
        std::vector<unsigned int> batchEncodedData(myBatchsize * encodedSequencePitchInInts2Bit);
        std::vector<char> batchQualities(myBatchsize * qualityPitchInBytes);
        std::vector<int> batchReadLengths(myBatchsize);    

        while(!(readIdGenerator.empty())){

            batchReadIds.resize(myBatchsize);

            readIdGenerator.process_next_n(
                myBatchsize, 
                [&](auto begin, auto end){
                    auto readIdsEnd = std::copy(begin, end, batchReadIds.begin());
                    batchReadIds.erase(readIdsEnd, batchReadIds.end());
                }
            );

            if(batchReadIds.empty()){
                continue;
            }

            //collect input data of all reads in batch

            readStorage.gatherSequenceLengths(
                batchReadLengths.data(),
                batchReadIds.data(),
                batchReadIds.size()
            );

            readStorage.gatherSequences(
                batchEncodedData.data(),
                encodedSequencePitchInInts2Bit,
                batchReadIds.data(),
                batchReadIds.size()
            );

            if(programOptions.useQualityScores){
                readStorage.gatherQualities(
                    batchQualities.data(),
                    qualityPitchInBytes,
                    batchReadIds.data(),
                    batchReadIds.size()
                );
            }

            CorrectionOutput correctionOutput;

            if(readStorage.isPairedEnd()){
                assert(batchReadIds.size() % 2 == 0);

                CpuErrorCorrectorMultiInput input{};
                input.anchorLengths.resize(batchReadIds.size());
                input.anchorReadIds.resize(batchReadIds.size());
                input.encodedAnchors.resize(batchReadIds.size());
                input.anchorQualityscores.resize(batchReadIds.size());

                for(size_t i = 0; i < batchReadIds.size(); i++){
                    input.anchorReadIds[i] = batchReadIds[i];
                    input.encodedAnchors[i] = batchEncodedData.data() + i * encodedSequencePitchInInts2Bit;
                    input.anchorQualityscores[i] = batchQualities.data() + i * qualityPitchInBytes;
                    input.anchorLengths[i] = batchReadLengths[i];
                }

                auto outputs = errorCorrector.processMulti(input);

                for(auto& output : outputs){

                    if(output.hasAnchorCorrection){
                        correctionOutput.anchorCorrections.emplace_back(std::move(output.anchorCorrection));
                    }

                    for(auto& tmp : output.candidateCorrections){
                        correctionOutput.candidateCorrections.emplace_back(std::move(tmp));
                    }
                }

            }else{

                for(size_t i = 0; i < batchReadIds.size(); i++){
                    const read_number readId = batchReadIds[i];

                    CpuErrorCorrectorInput input;
                    input.anchorReadId = readId;
                    input.encodedAnchor = batchEncodedData.data() + i * encodedSequencePitchInInts2Bit;
                    input.anchorQualityscores = batchQualities.data() + i * qualityPitchInBytes;
                    input.anchorLength = batchReadLengths[i];

                    auto output = errorCorrector.process(input);

                    if(output.hasAnchorCorrection){
                        correctionOutput.anchorCorrections.emplace_back(std::move(output.anchorCorrection));
                    }

                    for(auto& tmp : output.candidateCorrections){
                        correctionOutput.candidateCorrections.emplace_back(std::move(tmp));
                    }
                }
            }

            EncodedCorrectionOutput encodedCorrectionOutput = correctionOutput;

            auto outputfunction = [
                &, 
                encodedCorrectionOutput = std::move(encodedCorrectionOutput)
            ](){
                std::vector<std::uint8_t> tempbuffer(256);

                auto saveEncodedCorrectedSequence = [&](const EncodedTempCorrectedSequence* encoded){
                    if(!(encoded->isHQ() && encoded->useEdits() && encoded->getNumEdits() == 0)){
                        const std::size_t serializedSize = encoded->getSerializedNumBytes();
                        tempbuffer.resize(serializedSize);

                        auto end = encoded->copyToContiguousMemory(tempbuffer.data(), tempbuffer.data() + tempbuffer.size());
                        assert(end != nullptr);

                        partialResults.insert(tempbuffer.data(), end);
                    }
                };

                const int numA = encodedCorrectionOutput.encodedAnchorCorrections.size();
                const int numC = encodedCorrectionOutput.encodedCandidateCorrections.size();

                for(int i = 0; i < numA; i++){
                    saveEncodedCorrectedSequence(
                        &encodedCorrectionOutput.encodedAnchorCorrections[i]
                    );
                }

                for(int i = 0; i < numC; i++){
                    saveEncodedCorrectedSequence(
                        &encodedCorrectionOutput.encodedCandidateCorrections[i]
                    );
                }
            };

            outputThread.enqueue(std::move(outputfunction));

            clfAgent.flush();

            progressThread.addProgress(batchReadIds.size()); 

            // if(omp_get_thread_num() == 0){
            //     std::cerr << "getCandidatesTimeTotal: " << errorCorrector.getTimings().getCandidatesTimeTotal.count() << "\n";
            // }
            
        } //while unprocessed reads exist loop end   

        #pragma omp critical
        {
            timingsOfAllThreads += errorCorrector.getTimings();            
        }

    } // parallel end

    progressThread.finished();

    outputThread.stopThread(BackgroundThread::StopType::FinishAndStop);

    #ifdef ENABLE_CPU_CORRECTOR_TIMING

    auto totalDurationOfThreads = timingsOfAllThreads.getSumOfDurations();
    const int numThreads = programOptions.threads;

    auto printDuration = [&](const auto& name, const auto& duration){
        std::cout << "# average time per thread ("<< name << "): "
                  << duration.count() / numThreads  << " s. "
                  << (100.0 * duration.count() / totalDurationOfThreads.count()) << " %."<< std::endl;
    };

    #define printme(x) printDuration((#x), timingsOfAllThreads.x);

    printme(getAnchorSequenceDataTimeTotal);
    printme(getCandidatesTimeTotal);
    printme(copyCandidateDataToBufferTimeTotal);
    printme(getAlignmentsTimeTotal);
    printme(findBestAlignmentDirectionTimeTotal);
    printme(gatherBestAlignmentDataTimeTotal);
    printme(mismatchRatioFilteringTimeTotal);
    printme(compactBestAlignmentDataTimeTotal);
    printme(fetchQualitiesTimeTotal);
    printme(makeCandidateStringsTimeTotal);
    printme(msaAddSequencesTimeTotal);
    printme(msaFindConsensusTimeTotal);
    printme(msaMinimizationTimeTotal);
    printme(msaCorrectAnchorTimeTotal);
    printme(msaCorrectCandidatesTimeTotal);

    #undef printme

    #endif

    return partialResults;
}

} //namespace cpu

} //namespace care


