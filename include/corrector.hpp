#ifndef CARE_CORRECTOR_HPP
#define CARE_CORRECTOR_HPP

#include <config.hpp>

#include <options.hpp>

#include <cpuminhasher.hpp>

#include <cpureadstorage.hpp>
#include <cpu_alignment.hpp>
#include <alignmentorientation.hpp>
#include <msa.hpp>
#include <qualityscoreweights.hpp>
#include <classification.hpp>
#include <hostdevicefunctions.cuh>
#include <cpucorrectortask.hpp>
#include <correctedsequence.hpp>
#include <hostdevicefunctions.cuh>
#include <corrector_common.hpp>
#include <cpucorrectortask.hpp>
#include <util.hpp>

#include <cstddef>
#include <memory>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>

#include <omp.h>

namespace care{

    struct CpuErrorCorrectorOutput {
        bool hasAnchorCorrection{};
        TempCorrectedSequence anchorCorrection{};
        std::vector<TempCorrectedSequence> candidateCorrections{};
    };    


class CpuErrorCorrector{
public:

    struct MultiCandidateIds{
        std::vector<read_number> candidateReadIds;
        std::vector<int> numCandidatesPerAnchor;
        std::vector<int> numCandidatesPerAnchorPS;
        std::vector<bool> isPairedCandidate;
    };

    struct MultiCandidateData{
        std::vector<int> candidateLengths;
        std::vector<unsigned int> encodedCandidates;
        std::vector<char> candidateQualities;
    };

    struct TimeMeasurements{
        std::chrono::duration<double> getAnchorSequenceDataTimeTotal{0};
        std::chrono::duration<double> getCandidatesTimeTotal{0};
        std::chrono::duration<double> copyCandidateDataToBufferTimeTotal{0};
        std::chrono::duration<double> getAlignmentsTimeTotal{0};
        std::chrono::duration<double> findBestAlignmentDirectionTimeTotal{0};
        std::chrono::duration<double> gatherBestAlignmentDataTimeTotal{0};
        std::chrono::duration<double> mismatchRatioFilteringTimeTotal{0};
        std::chrono::duration<double> compactBestAlignmentDataTimeTotal{0};
        std::chrono::duration<double> fetchQualitiesTimeTotal{0};
        std::chrono::duration<double> makeCandidateStringsTimeTotal{0};
        std::chrono::duration<double> msaAddSequencesTimeTotal{0};
        std::chrono::duration<double> msaFindConsensusTimeTotal{0};
        std::chrono::duration<double> msaMinimizationTimeTotal{0};
        std::chrono::duration<double> msaCorrectAnchorTimeTotal{0};
        std::chrono::duration<double> msaCorrectCandidatesTimeTotal{0};

        TimeMeasurements& operator+=(const TimeMeasurements& rhs) noexcept{
            getAnchorSequenceDataTimeTotal += rhs.getAnchorSequenceDataTimeTotal;
            getCandidatesTimeTotal += rhs.getCandidatesTimeTotal;
            copyCandidateDataToBufferTimeTotal += rhs.copyCandidateDataToBufferTimeTotal;
            getAlignmentsTimeTotal += rhs.getAlignmentsTimeTotal;
            findBestAlignmentDirectionTimeTotal += rhs.findBestAlignmentDirectionTimeTotal;
            gatherBestAlignmentDataTimeTotal += rhs.gatherBestAlignmentDataTimeTotal;
            mismatchRatioFilteringTimeTotal += rhs.mismatchRatioFilteringTimeTotal;
            compactBestAlignmentDataTimeTotal += rhs.compactBestAlignmentDataTimeTotal;
            fetchQualitiesTimeTotal += rhs.fetchQualitiesTimeTotal;
            makeCandidateStringsTimeTotal += rhs.makeCandidateStringsTimeTotal;
            msaAddSequencesTimeTotal += rhs.msaAddSequencesTimeTotal;
            msaFindConsensusTimeTotal += rhs.msaFindConsensusTimeTotal;
            msaMinimizationTimeTotal += rhs.msaMinimizationTimeTotal;
            msaCorrectAnchorTimeTotal += rhs.msaCorrectAnchorTimeTotal;
            msaCorrectCandidatesTimeTotal += rhs.msaCorrectCandidatesTimeTotal;

            return *this;
        }

        std::chrono::duration<double> getSumOfDurations() const noexcept{
            std::chrono::duration<double> sum = getAnchorSequenceDataTimeTotal
                                            + getCandidatesTimeTotal
                                            + copyCandidateDataToBufferTimeTotal
                                            + getAlignmentsTimeTotal
                                            + findBestAlignmentDirectionTimeTotal
                                            + gatherBestAlignmentDataTimeTotal
                                            + mismatchRatioFilteringTimeTotal
                                            + compactBestAlignmentDataTimeTotal
                                            + fetchQualitiesTimeTotal
                                            + makeCandidateStringsTimeTotal
                                            + msaAddSequencesTimeTotal
                                            + msaFindConsensusTimeTotal
                                            + msaMinimizationTimeTotal
                                            + msaCorrectAnchorTimeTotal
                                            + msaCorrectCandidatesTimeTotal;
            return sum;
        }
    };   



    CpuErrorCorrector() = default;
    CpuErrorCorrector(
        std::size_t encodedSequencePitchInInts_,
        std::size_t decodedSequencePitchInBytes_,
        std::size_t qualityPitchInBytes_,
        const ProgramOptions& programOptions_,
        const CpuMinhasher& minhasher_,
        const CpuReadStorage& readStorage_,
        ReadCorrectionFlags& correctionFlags_,
        ClfAgent& clfAgent_
   ) : encodedSequencePitchInInts(encodedSequencePitchInInts_),
        decodedSequencePitchInBytes(decodedSequencePitchInBytes_),
        qualityPitchInBytes(qualityPitchInBytes_),
        programOptions(&programOptions_),
        minhasher{&minhasher_},
        minhashHandle{minhasher->makeMinhasherHandle()},
        readStorage{&readStorage_},
        correctionFlags(&correctionFlags_),
        clfAgent(&clfAgent_),
        qualityCoversion(std::make_unique<cpu::QualityScoreConversion>())
    {

    }

    ~CpuErrorCorrector(){
        minhasher->destroyHandle(minhashHandle);
    }

    CpuErrorCorrectorOutput process(const CpuErrorCorrectorInput input){
        CpuErrorCorrectorTask task = makeTask(input);
    
        TimeMeasurements timings;

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        auto tpa = std::chrono::system_clock::now();
        #endif

        determineCandidateReadIds(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.getCandidatesTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif

        computePairFlags(task);


        if(task.candidateReadIds.size() == 0){
            //return uncorrected anchor
            return CpuErrorCorrectorOutput{};
        }

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif

        getCandidateSequenceData(task);

        computeReverseComplementCandidates(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.copyCandidateDataToBufferTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif


        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif

        getCandidateAlignments(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.getAlignmentsTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif


        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif

        filterCandidatesByAlignmentFlag(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.gatherBestAlignmentDataTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif


        if(task.candidateReadIds.size() == 0){
            //return uncorrected anchor
            return CpuErrorCorrectorOutput{};
        }

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif

        filterCandidatesByAlignmentMismatchRatio(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.mismatchRatioFilteringTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif


        if(task.candidateReadIds.size() == 0){
            //return uncorrected anchor
            return CpuErrorCorrectorOutput{};
        }

        if(programOptions->useQualityScores){

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            getCandidateQualities(task);
            reverseQualitiesOfRCAlignments(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.fetchQualitiesTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif

        }

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif

        makeCandidateStrings(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.makeCandidateStringsTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif


        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif

        alignmentsComputeWeightsAndAoStoSoA(task);

        buildMultipleSequenceAlignment(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.msaFindConsensusTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif


        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif

        refineMSA(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.msaMinimizationTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif


        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif

        correctAnchor(task);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.msaCorrectAnchorTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif

        if(task.anchorCorrection.isCorrected){
            if(task.msaProperties.isHQ){
                correctionFlags->setCorrectedAsHqAnchor(task.input.anchorReadId);
            }
        }else{
            correctionFlags->setCouldNotBeCorrectedAsAnchor(task.input.anchorReadId);
        }


        if(task.msaProperties.isHQ && programOptions->correctCandidates){

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            correctCandidates(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.msaCorrectCandidatesTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif

        }

        CpuErrorCorrectorOutput correctionOutput = makeOutputOfTask(task);

        totalTime += timings;

        return correctionOutput;
    }

    const TimeMeasurements& getTimings() const noexcept{
        return totalTime;
    }

    std::stringstream& getMlStreamAnchor(){
        return ml_stream_anchor;
    }

    std::stringstream& getMlStreamCandidates(){
        return ml_stream_cands;
    }
    std::vector<CpuErrorCorrectorOutput> processMulti(const CpuErrorCorrectorMultiInput input){
        //if paired end, both reads of a pair must be processed simultaneously
        //paired anchors are stored consecutively
        assert(!readStorage->isPairedEnd() || input.anchorReadIds.size() % 2 == 0);

        const int numAnchors = input.anchorReadIds.size();
        if(numAnchors == 0){
            return {};
        }

        std::vector<CpuErrorCorrectorOutput> resultVector;
        resultVector.reserve(numAnchors);

        TimeMeasurements timings;

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        auto tpa = std::chrono::system_clock::now();
        #endif

        MultiCandidateIds multiIds = determineCandidateReadIds(input);

        // if(omp_get_thread_num() == 0){
        //     std::cerr << "get candidates: " << (std::chrono::system_clock::now() - tpa).count() << ", num candidates: " << multiIds.candidateReadIds.size() << "\n";
        // }

        computePairFlags(multiIds);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.getCandidatesTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        tpa = std::chrono::system_clock::now();
        #endif
        MultiCandidateData multiCandidates = getCandidateSequencesData(multiIds);

        #ifdef ENABLE_CPU_CORRECTOR_TIMING
        timings.copyCandidateDataToBufferTimeTotal += std::chrono::system_clock::now() - tpa;
        #endif

        for(int anchorIndex = 0; anchorIndex < numAnchors; anchorIndex++){

            CpuErrorCorrectorTask task = makeTask(input, multiIds, multiCandidates, anchorIndex);

            if(task.candidateReadIds.size() == 0){
                //return uncorrected anchor
                resultVector.emplace_back();
                continue;
            }

            computeReverseComplementCandidates(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            getCandidateAlignments(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.getAlignmentsTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif


            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            filterCandidatesByAlignmentFlag(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.gatherBestAlignmentDataTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif


            if(task.candidateReadIds.size() == 0){
                //return uncorrected anchor
                resultVector.emplace_back();
                continue;
            }

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            filterCandidatesByAlignmentMismatchRatio(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.mismatchRatioFilteringTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif


            if(task.candidateReadIds.size() == 0){
                //return uncorrected anchor
                resultVector.emplace_back();
                continue;
            }

            if(programOptions->useQualityScores){

                #ifdef ENABLE_CPU_CORRECTOR_TIMING
                tpa = std::chrono::system_clock::now();
                #endif

                getQualitiesFromMultiCandidates(task, multiIds, multiCandidates, anchorIndex);

                reverseQualitiesOfRCAlignments(task);

                #ifdef ENABLE_CPU_CORRECTOR_TIMING
                timings.fetchQualitiesTimeTotal += std::chrono::system_clock::now() - tpa;
                #endif

            }

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            makeCandidateStrings(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.makeCandidateStringsTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif


            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            alignmentsComputeWeightsAndAoStoSoA(task);

            buildMultipleSequenceAlignment(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.msaFindConsensusTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif


            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            refineMSA(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.msaMinimizationTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif


            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            tpa = std::chrono::system_clock::now();
            #endif

            correctAnchor(task);

            #ifdef ENABLE_CPU_CORRECTOR_TIMING
            timings.msaCorrectAnchorTimeTotal += std::chrono::system_clock::now() - tpa;
            #endif

            if(task.anchorCorrection.isCorrected){
                if(task.msaProperties.isHQ){
                    correctionFlags->setCorrectedAsHqAnchor(task.input.anchorReadId);
                }
            }else{
                correctionFlags->setCouldNotBeCorrectedAsAnchor(task.input.anchorReadId);
            }


            if(task.msaProperties.isHQ && programOptions->correctCandidates){

                #ifdef ENABLE_CPU_CORRECTOR_TIMING
                tpa = std::chrono::system_clock::now();
                #endif

                correctCandidates(task);

                #ifdef ENABLE_CPU_CORRECTOR_TIMING
                timings.msaCorrectCandidatesTimeTotal += std::chrono::system_clock::now() - tpa;
                #endif

            }

            resultVector.emplace_back(makeOutputOfTask(task));

        }

        totalTime += timings;

        return resultVector;
    }   

private:

    CpuErrorCorrectorTask makeTask(const CpuErrorCorrectorInput& input){
        CpuErrorCorrectorTask task;
        task.active = true;
        task.input = input;
        task.multipleSequenceAlignment.setQualityConversion(qualityCoversion.get());

        const int length = input.anchorLength;

        //decode anchor
        task.decodedAnchor.resize(length);
        SequenceHelpers::decode2BitSequence(task.decodedAnchor.data(), input.encodedAnchor, length);

        return task;
    }

    CpuErrorCorrectorTask makeTask(const CpuErrorCorrectorMultiInput& multiinput, const MultiCandidateIds& multiids, const MultiCandidateData& multicandidateData, int index){
        CpuErrorCorrectorInput input;
        input.anchorLength = multiinput.anchorLengths[index];
        input.anchorReadId = multiinput.anchorReadIds[index];
        input.encodedAnchor = multiinput.encodedAnchors[index];
        input.anchorQualityscores = multiinput.anchorQualityscores[index];

        CpuErrorCorrectorTask task = makeTask(input);

        const int offsetBegin = multiids.numCandidatesPerAnchorPS[index];
        const int offsetEnd = multiids.numCandidatesPerAnchorPS[index + 1];

        task.candidateReadIds.insert(
            task.candidateReadIds.end(),
            multiids.candidateReadIds.begin() + offsetBegin,
            multiids.candidateReadIds.begin() + offsetEnd
        );
        task.candidateSequencesLengths.insert(
            task.candidateSequencesLengths.end(),
            multicandidateData.candidateLengths.begin() + offsetBegin,
            multicandidateData.candidateLengths.begin() + offsetEnd
        );
        task.candidateSequencesData.insert(
            task.candidateSequencesData.end(),
            multicandidateData.encodedCandidates.begin() + encodedSequencePitchInInts * offsetBegin,
            multicandidateData.encodedCandidates.begin() + encodedSequencePitchInInts * offsetEnd
        );
        task.isPairedCandidate.insert(
            task.isPairedCandidate.end(),
            multiids.isPairedCandidate.begin() + offsetBegin,
            multiids.isPairedCandidate.begin() + offsetEnd
        );

        assert(task.isPairedCandidate.size() == task.candidateReadIds.size());
        
        return task;
    }


    void determineCandidateReadIds(CpuErrorCorrectorTask& task) const{

        task.candidateReadIds.clear();

        const read_number readId = task.input.anchorReadId;

        //const bool containsN = readProvider->readContainsN(readId);
        bool containsN = false;
        readStorage->areSequencesAmbiguous(&containsN, &readId, 1);

        //exclude anchors with ambiguous bases
        if(!(programOptions->excludeAmbiguousReads && containsN)){

            assert(task.input.anchorLength == int(task.decodedAnchor.size()));

            // candidateIdsProvider->getCandidates(
            //     task.candidateReadIds,
            //     task.decodedAnchor.data(),
            //     task.decodedAnchor.size()
            // );
            int numPerSequence = 0;
            int maxnumcandidates = 0;
            std::array<int, 2> offsets;

            minhasher->determineNumValues(
                minhashHandle,
                task.input.encodedAnchor,
                encodedSequencePitchInInts,
                &task.input.anchorLength,
                1,
                &numPerSequence,
                maxnumcandidates
            );

            task.candidateReadIds.resize(maxnumcandidates);

            minhasher->retrieveValues(
                minhashHandle,
                1,
                maxnumcandidates,
                task.candidateReadIds.data(),
                &numPerSequence,
                offsets.data()
            );

            std::sort(task.candidateReadIds.begin(), task.candidateReadIds.end());
            task.candidateReadIds.erase(
                std::unique(task.candidateReadIds.begin(), task.candidateReadIds.end()),
                task.candidateReadIds.end()
            );

            //remove self
            auto readIdPos = std::lower_bound(
                task.candidateReadIds.begin(),
                task.candidateReadIds.end(),
                readId
            );

            if(readIdPos != task.candidateReadIds.end() && *readIdPos == readId){
                task.candidateReadIds.erase(readIdPos);
            }

            auto resultsEnd = task.candidateReadIds.end();
            //exclude candidates with ambiguous bases

            if(programOptions->excludeAmbiguousReads){
                resultsEnd = std::remove_if(
                    task.candidateReadIds.begin(),
                    task.candidateReadIds.end(),
                    [&](read_number readId){
                        bool isAmbig = false;
                        readStorage->areSequencesAmbiguous(&isAmbig, &readId, 1);
                        return isAmbig;
                    }
                );
            }

            task.candidateReadIds.erase(resultsEnd, task.candidateReadIds.end());
        }
    }



    MultiCandidateIds determineCandidateReadIds(const CpuErrorCorrectorMultiInput& multiInput) const{
        const int numAnchors = multiInput.anchorReadIds.size();

        MultiCandidateIds multiCandidateIds;
        multiCandidateIds.numCandidatesPerAnchor.resize(numAnchors, 0);
        multiCandidateIds.numCandidatesPerAnchorPS.resize(numAnchors + 1, 0);

        for(int i = 0; i < numAnchors; i++){
            const read_number readId = multiInput.anchorReadIds[i];
            const int readlength = multiInput.anchorLengths[i];

            std::vector<char> decodedAnchor(readlength);
            SequenceHelpers::decode2BitSequence(decodedAnchor.data(), multiInput.encodedAnchors[i], readlength);

            bool containsN = false;
            readStorage->areSequencesAmbiguous(&containsN, &readId, 1);

            std::vector<read_number> candidateIds;

            //exclude anchors with ambiguous bases
            if(!(programOptions->excludeAmbiguousReads && containsN)){

                // candidateIdsProvider->getCandidates(
                //     candidateIds,
                //     decodedAnchor.data(),
                //     decodedAnchor.size()
                // );

                int numPerSequence = 0;
                int maxnumcandidates = 0;
                std::array<int, 2> offsets;

                minhasher->determineNumValues(
                    minhashHandle,
                    multiInput.encodedAnchors[i],
                    encodedSequencePitchInInts,
                    &readlength,
                    1,
                    &numPerSequence,
                    maxnumcandidates
                );

                candidateIds.resize(maxnumcandidates);

                minhasher->retrieveValues(
                    minhashHandle,
                    1,
                    maxnumcandidates,
                    candidateIds.data(),
                    &numPerSequence,
                    offsets.data()
                );

                std::sort(candidateIds.begin(), candidateIds.end());
                candidateIds.erase(
                    std::unique(candidateIds.begin(), candidateIds.end()),
                    candidateIds.end()
                );

                //remove self
                auto readIdPos = std::lower_bound(
                    candidateIds.begin(),
                    candidateIds.end(),
                    readId
                );

                if(readIdPos != candidateIds.end() && *readIdPos == readId){
                    candidateIds.erase(readIdPos);
                }

                auto resultsEnd = candidateIds.end();
                //exclude candidates with ambiguous bases

                if(programOptions->excludeAmbiguousReads){
                    resultsEnd = std::remove_if(
                        candidateIds.begin(),
                        candidateIds.end(),
                        [&](read_number readId){
                            bool isAmbig = false;
                            readStorage->areSequencesAmbiguous(&isAmbig, &readId, 1);
                            return isAmbig;
                        }
                    );
                }

                candidateIds.erase(resultsEnd, candidateIds.end());                
            }
        
            multiCandidateIds.numCandidatesPerAnchor[i] = candidateIds.size();
            multiCandidateIds.numCandidatesPerAnchorPS[i + 1] = 
                multiCandidateIds.numCandidatesPerAnchorPS[i] + candidateIds.size();
            multiCandidateIds.candidateReadIds.insert(multiCandidateIds.candidateReadIds.end(), candidateIds.begin(), candidateIds.end());
        
            multiCandidateIds.isPairedCandidate.resize(multiCandidateIds.candidateReadIds.size(), false);

            assert(multiCandidateIds.isPairedCandidate.size() == multiCandidateIds.candidateReadIds.size());
        }

        return multiCandidateIds;
    }


    void computePairFlags(CpuErrorCorrectorTask& task) const{
        //pair flags cannot be computed for a single task

        task.isPairedCandidate.resize(task.candidateReadIds.size());
        std::fill(task.isPairedCandidate.begin(), task.isPairedCandidate.end(), false);
    }

    void computePairFlags(MultiCandidateIds& multiIds) const{
        if(readStorage->isPairedEnd()){
            const int numAnchors = multiIds.numCandidatesPerAnchor.size();
            assert(numAnchors % 2 == 0);

            const int numAnchorPairs = numAnchors / 2;

            assert(multiIds.isPairedCandidate.size() == multiIds.candidateReadIds.size());
            assert(multiIds.numCandidatesPerAnchorPS.size() == multiIds.numCandidatesPerAnchor.size() + 1);

            for(int ap = 0; ap < numAnchorPairs; ap++){
                const int begin1 = multiIds.numCandidatesPerAnchorPS[2*ap + 0];
                const int end1 = multiIds.numCandidatesPerAnchorPS[2*ap + 1];
                const int begin2 = multiIds.numCandidatesPerAnchorPS[2*ap + 1];
                const int end2 = multiIds.numCandidatesPerAnchorPS[2*ap + 2];

                // assert(std::is_sorted(pairIds + begin1, pairIds + end1));
                // assert(std::is_sorted(pairIds + begin2, pairIds + end2));

                std::vector<int> pairedPositions(std::min(end1-begin1, end2-begin2));
                std::vector<int> pairedPositions2(std::min(end1-begin1, end2-begin2));

                auto endIters = findPositionsOfPairedReadIds(
                    multiIds.candidateReadIds.data() + begin1,
                    multiIds.candidateReadIds.data() + end1,
                    multiIds.candidateReadIds.data() + begin2,
                    multiIds.candidateReadIds.data() + end2,
                    pairedPositions.begin(),
                    pairedPositions2.begin()
                );

                pairedPositions.erase(endIters.first, pairedPositions.end());
                pairedPositions2.erase(endIters.second, pairedPositions2.end());
                for(auto i : pairedPositions){
                    multiIds.isPairedCandidate[begin1 + i] = true;
                }
                for(auto i : pairedPositions2){
                    multiIds.isPairedCandidate[begin2 + i] = true;
                }            
            }
        }else{
            assert(multiIds.isPairedCandidate.size() == multiIds.candidateReadIds.size());
            std::fill(multiIds.isPairedCandidate.begin(), multiIds.isPairedCandidate.end(), false);
        }
    }



    MultiCandidateData getCandidateSequencesData(const MultiCandidateIds& multiIds) const{
        const int numIds = multiIds.candidateReadIds.size();
        if(numIds == 0){
            return MultiCandidateData{};
        }

        MultiCandidateData multiData;
        multiData.candidateLengths.resize(numIds);
        multiData.encodedCandidates.resize(numIds * encodedSequencePitchInInts);
        multiData.candidateQualities.resize(numIds * qualityPitchInBytes);

        readStorage->gatherSequenceLengths(
            multiData.candidateLengths.data(),
            multiIds.candidateReadIds.data(),
            numIds
        );

        readStorage->gatherSequences(
            multiData.encodedCandidates.data(),
            encodedSequencePitchInInts,
            multiIds.candidateReadIds.data(),
            numIds
        ); 

        if(programOptions->useQualityScores){

            readStorage->gatherQualities(
                multiData.candidateQualities.data(),
                qualityPitchInBytes,
                multiIds.candidateReadIds.data(),
                numIds
            ); 
        }

        return multiData;
    }


    void getQualitiesFromMultiCandidates(CpuErrorCorrectorTask& task, const MultiCandidateIds& multiids, const MultiCandidateData& multicandidateData, int index) const{
        const int offsetBegin = multiids.numCandidatesPerAnchorPS[index];
        const int offsetEnd = multiids.numCandidatesPerAnchorPS[index + 1];

        const int numCandidates = task.candidateReadIds.size();
        task.candidateQualities.resize(qualityPitchInBytes * numCandidates);

        auto first1 = task.candidateReadIds.cbegin();
        auto last1 = task.candidateReadIds.cend();
        auto first2 = multiids.candidateReadIds.cbegin() + offsetBegin;
        auto last2 = multiids.candidateReadIds.cbegin() + offsetEnd;

        auto qualOutputIter = task.candidateQualities.begin();
        auto qualInputIter = multicandidateData.candidateQualities.cbegin() + offsetBegin * qualityPitchInBytes;

        //copy quality scores of candidates which have not been removed (set_intersection, range 1 is a subset of range2)
        while (first1 != last1 && first2 != last2) {
            if (*first1 < *first2) {
                ++first1;
            } else  {
                if (!(*first2 < *first1)) {
                    qualOutputIter = std::copy_n(qualInputIter, qualityPitchInBytes, qualOutputIter);
                }
                ++first2;
                qualInputIter += qualityPitchInBytes;
            }
        }
    }


    //Gets forward sequence and reverse complement sequence of each candidate
    void getCandidateSequenceData(CpuErrorCorrectorTask& task) const{

        const int numCandidates = task.candidateReadIds.size();
        
        if(numCandidates == 0) return;

        task.candidateSequencesLengths.resize(numCandidates);
        task.candidateSequencesData.clear();
        task.candidateSequencesData.resize(size_t(encodedSequencePitchInInts) * numCandidates, 0);
        task.candidateSequencesRevcData.clear();
        task.candidateSequencesRevcData.resize(size_t(encodedSequencePitchInInts) * numCandidates, 0);

        readStorage->gatherSequenceLengths(
            task.candidateSequencesLengths.data(),
            task.candidateReadIds.data(),
            numCandidates
        );

        readStorage->gatherSequences(
            task.candidateSequencesData.data(),
            encodedSequencePitchInInts,
            task.candidateReadIds.data(),
            numCandidates
        );        
    }

    void computeReverseComplementCandidates(CpuErrorCorrectorTask& task){
        const int numCandidates = task.candidateReadIds.size();
        task.candidateSequencesRevcData.resize(task.candidateSequencesData.size());

        for(int i = 0; i < numCandidates; i++){
            const unsigned int* const seqPtr = task.candidateSequencesData.data() 
                                                + std::size_t(encodedSequencePitchInInts) * i;
            unsigned int* const seqrevcPtr = task.candidateSequencesRevcData.data() 
                                                + std::size_t(encodedSequencePitchInInts) * i;

            SequenceHelpers::reverseComplementSequence2Bit(
                seqrevcPtr,  
                seqPtr,
                task.candidateSequencesLengths[i]
            );
        }
    } 

    //compute alignments between anchor sequence and candidate sequences
    void getCandidateAlignments(CpuErrorCorrectorTask& task) const{
        const int numCandidates = task.candidateReadIds.size();

        task.alignments.resize(numCandidates);
        task.revcAlignments.resize(numCandidates);
        task.alignmentFlags.resize(numCandidates);

        cpu::shd::cpuShiftedHammingDistancePopcount2Bit(
            alignmentHandle,
            task.alignments.begin(),
            task.input.encodedAnchor,
            task.input.anchorLength,
            task.candidateSequencesData.data(),
            encodedSequencePitchInInts,
            task.candidateSequencesLengths.data(),
            numCandidates,
            programOptions->min_overlap,
            programOptions->maxErrorRate,
            programOptions->min_overlap_ratio
        );

        cpu::shd::cpuShiftedHammingDistancePopcount2Bit(
            alignmentHandle,
            task.revcAlignments.begin(),
            task.input.encodedAnchor,
            task.input.anchorLength,
            task.candidateSequencesRevcData.data(),
            encodedSequencePitchInInts,
            task.candidateSequencesLengths.data(),
            numCandidates,
            programOptions->min_overlap,
            programOptions->maxErrorRate,
            programOptions->min_overlap_ratio
        );

        //decide whether to keep forward or reverse complement

        for(int i = 0; i < numCandidates; i++){
            const auto& forwardAlignment = task.alignments[i];
            const auto& revcAlignment = task.revcAlignments[i];
            const int candidateLength = task.candidateSequencesLengths[i];

            AlignmentOrientation bestAlignmentFlag = care::chooseBestAlignmentOrientation(
                forwardAlignment,
                revcAlignment,
                task.input.anchorLength,
                candidateLength,
                programOptions->min_overlap_ratio,
                programOptions->min_overlap,
                programOptions->estimatedErrorrate
            );

            task.alignmentFlags[i] = bestAlignmentFlag;
        }
    }

    //remove candidates with alignment flag None
    void filterCandidatesByAlignmentFlag(CpuErrorCorrectorTask& task) const{

        const int numCandidates = task.candidateReadIds.size();

        int insertpos = 0;

        for(int i = 0; i < numCandidates; i++){

            const AlignmentOrientation flag = task.alignmentFlags[i];

            if(flag == AlignmentOrientation::Forward){
                task.candidateReadIds[insertpos] = task.candidateReadIds[i];
                std::copy_n(
                    task.candidateSequencesData.data() + i * size_t(encodedSequencePitchInInts),
                    encodedSequencePitchInInts,
                    task.candidateSequencesData.data() + insertpos * size_t(encodedSequencePitchInInts)
                );
                task.candidateSequencesLengths[insertpos] = task.candidateSequencesLengths[i];
                task.alignmentFlags[insertpos] = task.alignmentFlags[i];
                task.alignments[insertpos] = task.alignments[i];   
                task.isPairedCandidate[insertpos] = task.isPairedCandidate[i];             

                insertpos++;
            }else if(flag == AlignmentOrientation::ReverseComplement){
                task.candidateReadIds[insertpos] = task.candidateReadIds[i];
                std::copy_n(
                    task.candidateSequencesRevcData.data() + i * size_t(encodedSequencePitchInInts),
                    encodedSequencePitchInInts,
                    task.candidateSequencesData.data() + insertpos * size_t(encodedSequencePitchInInts)
                );
                task.candidateSequencesLengths[insertpos] = task.candidateSequencesLengths[i];
                task.alignmentFlags[insertpos] = task.alignmentFlags[i];
                task.alignments[insertpos] = task.revcAlignments[i];
                task.isPairedCandidate[insertpos] = task.isPairedCandidate[i];               

                insertpos++;
            }else{
                ;//AlignmentOrientation::None discard alignment
            }
        }

        task.candidateReadIds.erase(
            task.candidateReadIds.begin() + insertpos, 
            task.candidateReadIds.end()
        );
        task.candidateSequencesData.erase(
            task.candidateSequencesData.begin() + encodedSequencePitchInInts * insertpos, 
            task.candidateSequencesData.end()
        );
        task.candidateSequencesLengths.erase(
            task.candidateSequencesLengths.begin() + insertpos, 
            task.candidateSequencesLengths.end()
        );
        task.alignmentFlags.erase(
            task.alignmentFlags.begin() + insertpos, 
            task.alignmentFlags.end()
        );
        task.alignments.erase(
            task.alignments.begin() + insertpos, 
            task.alignments.end()
        );
        task.isPairedCandidate.erase(
            task.isPairedCandidate.begin() + insertpos, 
            task.isPairedCandidate.end()
        );

        task.revcAlignments.clear();
        task.candidateSequencesRevcData.clear();
    }

    //remove candidates with bad alignment mismatch ratio
    void filterCandidatesByAlignmentMismatchRatio(CpuErrorCorrectorTask& task) const{

        if(readStorage->isPairedEnd()){
            filterCandidatesByAlignmentMismatchRatioWithPairFlags(task);
            return;
        }
        
        auto lastResortFunc = [](){
            return false;
        };

        const float mismatchratioBaseFactor = programOptions->estimatedErrorrate * 1.0f;
        const float goodAlignmentsCountThreshold = programOptions->estimatedCoverage * programOptions->m_coverage;
        const int numCandidates = task.candidateReadIds.size();

        std::array<int, 3> counts({0,0,0});

        for(int i = 0; i < numCandidates; i++){
            const auto& alignment = task.alignments[i];
            const float mismatchratio = float(alignment.nOps) / float(alignment.overlap);

            if (mismatchratio < 2 * mismatchratioBaseFactor) {
                counts[0] += 1;
            }
            if (mismatchratio < 3 * mismatchratioBaseFactor) {
                counts[1] += 1;
            }
            if (mismatchratio < 4 * mismatchratioBaseFactor) {
                counts[2] += 1;
            }
        }

        float mismatchratioThreshold = std::numeric_limits<float>::lowest();
        if( std::any_of(counts.begin(), counts.end(), [](auto c){return c > 0;}) ){
            
            if (counts[0] >= goodAlignmentsCountThreshold) {
                mismatchratioThreshold = 2 * mismatchratioBaseFactor;
            } else if (counts[1] >= goodAlignmentsCountThreshold) {
                mismatchratioThreshold = 3 * mismatchratioBaseFactor;
            } else if (counts[2] >= goodAlignmentsCountThreshold) {
                mismatchratioThreshold = 4 * mismatchratioBaseFactor;
            } else {
                if(lastResortFunc()){
                    mismatchratioThreshold = 4 * mismatchratioBaseFactor;
                }else{
                    ; //mismatchratioThreshold = std::numeric_limits<float>::lowest();
                }
            }
        }

        int insertpos = 0;
        for(int i = 0; i < numCandidates; i++){
            const auto& alignment = task.alignments[i];
            const float mismatchratio = float(alignment.nOps) / float(alignment.overlap);
            const bool notremoved = mismatchratio < mismatchratioThreshold;

            if(notremoved){
                task.candidateReadIds[insertpos] = task.candidateReadIds[i];
                std::copy_n(
                    task.candidateSequencesData.data() + i * size_t(encodedSequencePitchInInts),
                    encodedSequencePitchInInts,
                    task.candidateSequencesData.data() + insertpos * size_t(encodedSequencePitchInInts)
                );
                task.candidateSequencesLengths[insertpos] = task.candidateSequencesLengths[i];
                task.alignmentFlags[insertpos] = task.alignmentFlags[i];
                task.alignments[insertpos] = task.alignments[i];
                task.isPairedCandidate[insertpos] = task.isPairedCandidate[i];

                insertpos++;
            }
        }

        task.candidateReadIds.erase(
            task.candidateReadIds.begin() + insertpos, 
            task.candidateReadIds.end()
        );
        task.candidateSequencesData.erase(
            task.candidateSequencesData.begin() + encodedSequencePitchInInts * insertpos, 
            task.candidateSequencesData.end()
        );
        task.candidateSequencesLengths.erase(
            task.candidateSequencesLengths.begin() + insertpos, 
            task.candidateSequencesLengths.end()
        );
        task.alignmentFlags.erase(
            task.alignmentFlags.begin() + insertpos, 
            task.alignmentFlags.end()
        );
        task.alignments.erase(
            task.alignments.begin() + insertpos, 
            task.alignments.end()
        );
        task.isPairedCandidate.erase(
            task.isPairedCandidate.begin() + insertpos, 
            task.isPairedCandidate.end()
        );
    }

    void filterCandidatesByAlignmentMismatchRatioWithPairFlags(CpuErrorCorrectorTask& task) const{
        const float threshold = programOptions->pairedFilterThreshold;
        const int numCandidates = task.candidateReadIds.size();

        int insertpos = 0;
        auto keep = [&](int i){
            task.candidateReadIds[insertpos] = task.candidateReadIds[i];
            std::copy_n(
                task.candidateSequencesData.data() + i * size_t(encodedSequencePitchInInts),
                encodedSequencePitchInInts,
                task.candidateSequencesData.data() + insertpos * size_t(encodedSequencePitchInInts)
            );
            task.candidateSequencesLengths[insertpos] = task.candidateSequencesLengths[i];
            task.alignmentFlags[insertpos] = task.alignmentFlags[i];
            task.alignments[insertpos] = task.alignments[i];
            task.isPairedCandidate[insertpos] = task.isPairedCandidate[i];

            insertpos++;
        };

        //paired candidates are kept unconditionally. others are filtered by mismatch ratio
        assert(task.isPairedCandidate.size() == task.candidateReadIds.size());

        for(int candidate_index = 0; candidate_index < numCandidates; candidate_index += 1) {
            if(!task.isPairedCandidate[candidate_index]){
                const auto& alignment = task.alignments[candidate_index];
                const float mismatchratio = float(alignment.nOps) / float(alignment.overlap);

                if(mismatchratio < threshold) {
                    keep(candidate_index);
                }
            }else{                
                keep(candidate_index);
            }
        }

        task.candidateReadIds.erase(
            task.candidateReadIds.begin() + insertpos, 
            task.candidateReadIds.end()
        );
        task.candidateSequencesData.erase(
            task.candidateSequencesData.begin() + encodedSequencePitchInInts * insertpos, 
            task.candidateSequencesData.end()
        );
        task.candidateSequencesLengths.erase(
            task.candidateSequencesLengths.begin() + insertpos, 
            task.candidateSequencesLengths.end()
        );
        task.alignmentFlags.erase(
            task.alignmentFlags.begin() + insertpos, 
            task.alignmentFlags.end()
        );
        task.alignments.erase(
            task.alignments.begin() + insertpos, 
            task.alignments.end()
        );
        task.isPairedCandidate.erase(
            task.isPairedCandidate.begin() + insertpos, 
            task.isPairedCandidate.end()
        );
    }

    //get quality scores of candidates with respect to alignment direction
    void getCandidateQualities(CpuErrorCorrectorTask& task) const{
        const int numCandidates = task.candidateReadIds.size();

        task.candidateQualities.resize(qualityPitchInBytes * numCandidates);

        readStorage->gatherQualities(
            task.candidateQualities.data(),
            qualityPitchInBytes,
            task.candidateReadIds.data(),
            numCandidates
        );         
    }

    void reverseQualitiesOfRCAlignments(CpuErrorCorrectorTask& task){
        const int numCandidates = task.candidateReadIds.size();

        //reverse quality scores of candidates with reverse complement alignment
        for(int c = 0; c < numCandidates; c++){
            if(task.alignmentFlags[c] == AlignmentOrientation::ReverseComplement){
                std::reverse(
                    task.candidateQualities.data() + c * size_t(qualityPitchInBytes),
                    task.candidateQualities.data() + (c+1) * size_t(qualityPitchInBytes)
                );
            }
        }             
    }

    //compute decoded candidate strings with respect to alignment direction
    void makeCandidateStrings(CpuErrorCorrectorTask& task) const{
        const int numCandidates = task.candidateReadIds.size();

        task.decodedCandidateSequences.resize(decodedSequencePitchInBytes * numCandidates);
        
        for(int i = 0; i < numCandidates; i++){
            const unsigned int* const srcptr = task.candidateSequencesData.data() + i * encodedSequencePitchInInts;
            char* const destptr = task.decodedCandidateSequences.data() + i * decodedSequencePitchInBytes;
            const int length = task.candidateSequencesLengths[i];

            SequenceHelpers::decode2BitSequence(
                destptr,
                srcptr,
                length
            );
        }
    }

    void alignmentsComputeWeightsAndAoStoSoA(CpuErrorCorrectorTask& task) const{
        const int numCandidates = task.candidateReadIds.size();

        task.alignmentShifts.resize(numCandidates);
        task.alignmentOps.resize(numCandidates);
        task.alignmentOverlaps.resize(numCandidates);
        task.alignmentWeights.resize(numCandidates);

        for(int i = 0; i < numCandidates; i++){
            task.alignmentShifts[i] = task.alignments[i].shift;
            task.alignmentOps[i] = task.alignments[i].nOps;
            task.alignmentOverlaps[i] = task.alignments[i].overlap;

            task.alignmentWeights[i] = calculateOverlapWeight(
                task.input.anchorLength,
                task.alignments[i].nOps, 
                task.alignments[i].overlap,
                programOptions->maxErrorRate
            );
        }
    }

    void buildMultipleSequenceAlignment(CpuErrorCorrectorTask& task) const{

        const int numCandidates = task.candidateReadIds.size();

        const char* const candidateQualityPtr = programOptions->useQualityScores ?
                                                task.candidateQualities.data()
                                                : nullptr;

        MultipleSequenceAlignment::InputData buildArgs;
        buildArgs.useQualityScores = programOptions->useQualityScores;
        buildArgs.anchorLength = task.input.anchorLength;
        buildArgs.nCandidates = numCandidates;
        buildArgs.candidatesPitch = decodedSequencePitchInBytes;
        buildArgs.candidateQualitiesPitch = qualityPitchInBytes;
        buildArgs.anchor = task.decodedAnchor.data();
        buildArgs.candidates = task.decodedCandidateSequences.data();
        buildArgs.anchorQualities = task.input.anchorQualityscores;
        buildArgs.candidateQualities = candidateQualityPtr;
        buildArgs.candidateLengths = task.candidateSequencesLengths.data();
        buildArgs.candidateShifts = task.alignmentShifts.data();
        buildArgs.candidateDefaultWeightFactors = task.alignmentWeights.data();
    
        task.multipleSequenceAlignment.build(buildArgs);
    }

    void refineMSA(CpuErrorCorrectorTask& task) const{

        constexpr int max_num_minimizations = 5;

        auto removeCandidatesOfDifferentRegion = [&](const auto& minimizationResult){
            const int numCandidates = task.candidateReadIds.size();

            int insertpos = 0;
            for(int i = 0; i < numCandidates; i++){
                if(!minimizationResult.differentRegionCandidate[i]){                        
                    //keep candidate

                    task.candidateReadIds[insertpos] = task.candidateReadIds[i];

                    std::copy_n(
                        task.candidateSequencesData.data() + i * size_t(encodedSequencePitchInInts),
                        encodedSequencePitchInInts,
                        task.candidateSequencesData.data() + insertpos * size_t(encodedSequencePitchInInts)
                    );

                    task.candidateSequencesLengths[insertpos] = task.candidateSequencesLengths[i];
                    task.alignmentFlags[insertpos] = task.alignmentFlags[i];
                    task.alignments[insertpos] = task.alignments[i]; 
                    task.alignmentOps[insertpos] = task.alignmentOps[i];
                    task.alignmentShifts[insertpos] = task.alignmentShifts[i];
                    task.alignmentOverlaps[insertpos] = task.alignmentOverlaps[i];
                    task.alignmentWeights[insertpos] = task.alignmentWeights[i];

                    if(programOptions->useQualityScores){
                        std::copy_n(
                            task.candidateQualities.data() + i * size_t(qualityPitchInBytes),
                            qualityPitchInBytes,
                            task.candidateQualities.data() + insertpos * size_t(qualityPitchInBytes)
                        );
                    }

                    std::copy_n(
                        task.decodedCandidateSequences.data() + i * size_t(decodedSequencePitchInBytes),
                        decodedSequencePitchInBytes,
                        task.decodedCandidateSequences.data() + insertpos * size_t(decodedSequencePitchInBytes)
                    );

                    insertpos++;
                }
            }

            task.candidateReadIds.erase(
                task.candidateReadIds.begin() + insertpos, 
                task.candidateReadIds.end()
            );
            task.candidateSequencesData.erase(
                task.candidateSequencesData.begin() + encodedSequencePitchInInts * insertpos, 
                task.candidateSequencesData.end()
            );
            task.candidateSequencesLengths.erase(
                task.candidateSequencesLengths.begin() + insertpos, 
                task.candidateSequencesLengths.end()
            );
            task.alignmentFlags.erase(
                task.alignmentFlags.begin() + insertpos, 
                task.alignmentFlags.end()
            );
            task.alignments.erase(
                task.alignments.begin() + insertpos, 
                task.alignments.end()
            );

            if(programOptions->useQualityScores){
                task.candidateQualities.erase(
                    task.candidateQualities.begin() + qualityPitchInBytes * insertpos, 
                    task.candidateQualities.end()
                );
            }
            
            task.decodedCandidateSequences.erase(
                task.decodedCandidateSequences.begin() + decodedSequencePitchInBytes * insertpos, 
                task.decodedCandidateSequences.end()
            );
            task.alignmentOps.erase(
                task.alignmentOps.begin() + insertpos, 
                task.alignmentOps.end()
            );
            task.alignmentShifts.erase(
                task.alignmentShifts.begin() + insertpos, 
                task.alignmentShifts.end()
            );
            task.alignmentOverlaps.erase(
                task.alignmentOverlaps.begin() + insertpos, 
                task.alignmentOverlaps.end()
            );
            task.alignmentWeights.erase(
                task.alignmentWeights.begin() + insertpos, 
                task.alignmentWeights.end()
            );
        };

        if(max_num_minimizations > 0){                

            for(int numIterations = 0; numIterations < max_num_minimizations; numIterations++){
                const auto minimizationResult = task.multipleSequenceAlignment.findCandidatesOfDifferentRegion(
                    programOptions->estimatedCoverage
                );

                if(minimizationResult.performedMinimization){
                    removeCandidatesOfDifferentRegion(minimizationResult);

                    //build minimized multiple sequence alignment
                    buildMultipleSequenceAlignment(task);
                }else{
                    break;
                }               
                
            }
        }      
    }

    void correctAnchorClassic(CpuErrorCorrectorTask& task) const{

        assert(programOptions->correctionType == CorrectionType::Classic);

        const int anchorColumnsBegin_incl = task.multipleSequenceAlignment.anchorColumnsBegin_incl;
        const int anchorColumnsEnd_excl = task.multipleSequenceAlignment.anchorColumnsEnd_excl;

        task.msaProperties = task.multipleSequenceAlignment.getMSAProperties(
            anchorColumnsBegin_incl,
            anchorColumnsEnd_excl,
            programOptions->estimatedErrorrate,
            programOptions->estimatedCoverage,
            programOptions->m_coverage
        );

        task.anchorCorrection = task.multipleSequenceAlignment.getCorrectedAnchor(
            task.msaProperties,
            programOptions->estimatedErrorrate,
            programOptions->estimatedCoverage,
            programOptions->m_coverage,
            programOptions->kmerlength,
            task.input.anchorReadId
        );        
    }       

    void correctAnchorClf(CpuErrorCorrectorTask& task) const
    {
        auto& msa = task.multipleSequenceAlignment;
        int anchor_b = msa.anchorColumnsBegin_incl;
        int anchor_e = msa.anchorColumnsEnd_excl;
        auto& cons = msa.consensus;
        auto& orig = task.decodedAnchor;
        auto& corr = task.anchorCorrection.correctedSequence;

        task.msaProperties = msa.getMSAProperties(
            anchor_b, anchor_e, programOptions->estimatedErrorrate, programOptions->estimatedCoverage,
            programOptions->m_coverage);

        corr.insert(0, cons.data()+anchor_b, task.input.anchorLength);
        if (!task.msaProperties.isHQ) {
            for (int i = 0; i < task.input.anchorLength; ++i) {
                if (orig[i] != cons[anchor_b+i] && !clfAgent->decide_anchor(task, i, *programOptions))
                {
                    corr[i] = orig[i];
                }
            }
        }

        task.anchorCorrection.isCorrected = true;
    }

    void correctAnchorPrint(CpuErrorCorrectorTask& task) const{
        auto& msa = task.multipleSequenceAlignment;
        int anchor_b = msa.anchorColumnsBegin_incl;
        int anchor_e = msa.anchorColumnsEnd_excl;
        auto& cons = msa.consensus;
        auto& orig = task.decodedAnchor;

        task.msaProperties = msa.getMSAProperties(
            anchor_b, anchor_e, programOptions->estimatedErrorrate, programOptions->estimatedCoverage,
            programOptions->m_coverage);

        if (!task.msaProperties.isHQ) {
            for (int i = 0; i < task.input.anchorLength; ++i) {
                if (orig[i] != cons[anchor_b+i]) {
                    clfAgent->print_anchor(task, i, *programOptions);
                }
            }
        }
        task.anchorCorrection.isCorrected = false;
    }

    void correctAnchor(CpuErrorCorrectorTask& task) const{
        if(programOptions->correctionType == CorrectionType::Classic)
            correctAnchorClassic(task);
        else if(programOptions->correctionType == CorrectionType::Forest)
            correctAnchorClf(task);
        else if (programOptions->correctionType == CorrectionType::Print)
            correctAnchorPrint(task);
    }

    void correctCandidatesClassic(CpuErrorCorrectorTask& task) const{

        task.candidateCorrections = task.multipleSequenceAlignment.getCorrectedCandidates(
            programOptions->estimatedErrorrate,
            programOptions->estimatedCoverage,
            programOptions->m_coverage,
            programOptions->new_columns_to_correct
        );
    }

    void correctCandidatesPrint(CpuErrorCorrectorTask& task) const{

        const auto& msa = task.multipleSequenceAlignment;
        const int anchor_begin = msa.anchorColumnsBegin_incl;
        const int anchor_end = msa.anchorColumnsEnd_excl;

        for(int cand = 0; cand < msa.nCandidates; ++cand) {
            const int cand_begin = msa.anchorColumnsBegin_incl + task.alignmentShifts[cand];
            const int cand_length = task.candidateSequencesLengths[cand];
            const int cand_end = cand_begin + cand_length;
            const int offset = cand * decodedSequencePitchInBytes;
            
            if(cand_begin >= anchor_begin - programOptions->new_columns_to_correct
                && cand_end <= anchor_end + programOptions->new_columns_to_correct)
            {
                for (int i = 0; i < cand_length; ++i) {
                    if (task.decodedCandidateSequences[offset+i] != msa.consensus[cand_begin+i]) {
                        clfAgent->print_cand(task, i, *programOptions, cand, offset);
                    }
                }
            }
        }
        task.candidateCorrections = std::vector<CorrectedCandidate>{};
    }

    void correctCandidatesClf(CpuErrorCorrectorTask& task) const 
    {
        const auto& msa = task.multipleSequenceAlignment;

        task.candidateCorrections = std::vector<CorrectedCandidate>{};

        // auto it = std::find(task.candidateReadIds.begin(), task.candidateReadIds.end(), 37);
        // if(it != task.candidateReadIds.end()){
        //     std::cerr << "found 37 at local position " << std::distance(task.candidateReadIds.begin(), it) << "\n";
        // }

        const int anchor_begin = msa.anchorColumnsBegin_incl;
        const int anchor_end = msa.anchorColumnsEnd_excl;
        for(int cand = 0; cand < msa.nCandidates; ++cand) {
            read_number candidateReadId = task.candidateReadIds[cand];

            //if this read has already been corrected as a high quality anchor, the candidate correction will not be used for output construction.
            //-> don't compute it.
            if(correctionFlags->isCorrectedAsHQAnchor(candidateReadId)){
                continue;
            }
            const int cand_begin = msa.anchorColumnsBegin_incl + task.alignmentShifts[cand];
            const int cand_length = task.candidateSequencesLengths[cand];
            const int cand_end = cand_begin + cand_length;
            const int offset = cand * decodedSequencePitchInBytes;

            // if(task.candidateReadIds[cand] == 37){
            //     std::cerr <<  "candidateIndex " << cand << "\n";
            //     std::cerr <<  "cand_begin " << cand_begin << "\n";
            //     std::cerr <<  "cand_end " << cand_end << "\n";
            //     std::cerr <<  "cand_length " << cand_length << "\n";
            //     std::cerr <<  "candidateReadId " << task.candidateReadIds[cand] << "\n";
            //     std::cerr <<  "anchorReadId " << task.input.anchorReadId << "\n";
            // }

            
            if(cand_begin >= anchor_begin - programOptions->new_columns_to_correct
                && cand_end <= anchor_end + programOptions->new_columns_to_correct)
            {

                task.candidateCorrections.emplace_back(cand, task.alignmentShifts[cand],
                    std::string(&msa.consensus[cand_begin], cand_length));

                // if(task.candidateReadIds[cand] == 37){
                //     std::cerr << "in range\n";
                //     std::cerr <<  "candidateIndex " << cand << "\n";
                //     std::cerr <<  "cand_begin " << cand_begin << "\n";
                //     std::cerr <<  "cand_end " << cand_end << "\n";
                //     std::cerr <<  "cand_length " << cand_length << "\n";
                //     std::cerr <<  "candidateReadId " << task.candidateReadIds[cand] << "\n";

                //     std::cerr << "decodedCandidate:\n";
                //     for(int k = 0; k < cand_length; k++){
                //         std::cerr << task.decodedCandidateSequences[offset+k];
                //     }
                //     std::cerr << "\n";

                //     std::cerr << "consensusCandidate:\n";
                //     for(int k = 0; k < cand_length; k++){
                //         std::cerr << task.candidateCorrections.back().sequence[k];
                //     }
                //     std::cerr << "\n";

                // }


                for (int i = 0; i < cand_length; ++i) {
                    if (task.decodedCandidateSequences[offset+i] != msa.consensus[cand_begin+i]
                        && !clfAgent->decide_cand(task, i, *programOptions, cand, offset))
                    {
                        task.candidateCorrections.back().sequence[i] = task.decodedCandidateSequences[offset+i];
                        //if(task.input.anchorReadId == 5 && task.candidateReadIds[cand] == 633) std::cerr << "checking position " << i << ". revert consensus\n";
                    }else{
                        // if (task.decodedCandidateSequences[offset+i] != msa.consensus[cand_begin+i]){
                        //     if(task.input.anchorReadId == 5 && task.candidateReadIds[cand] == 633) std::cerr << "checking position " << i << ". keep consensus\n";
                        // }
                    }
                }
            }else{
                // if(task.candidateReadIds[cand] == 37){
                //     std::cerr << "not in range with shift " << task.alignmentShifts[cand] << "\n";
                // }
            }
        }
    }

    void correctCandidates(CpuErrorCorrectorTask& task) const{
        switch (programOptions->correctionTypeCands) {
            case CorrectionType::Print:
                correctCandidatesPrint(task);
                break;
            case CorrectionType::Forest:
                correctCandidatesClf(task);
                break;
            case CorrectionType::Classic:
            default:
                correctCandidatesClassic(task);
        }
    }

    CpuErrorCorrectorOutput makeOutputOfTask(CpuErrorCorrectorTask& task) const{
        CpuErrorCorrectorOutput result;

        result.hasAnchorCorrection = task.anchorCorrection.isCorrected;

        if(result.hasAnchorCorrection){
            auto& correctedSequenceString = task.anchorCorrection.correctedSequence;
            const int correctedlength = correctedSequenceString.length();
            bool originalReadContainsN = false;
            readStorage->areSequencesAmbiguous(&originalReadContainsN, &task.input.anchorReadId, 1);
            
            TempCorrectedSequence tmp;
            
            if(!originalReadContainsN){
                const int maxEdits = correctedlength / 7;
                int edits = 0;
                for(int i = 0; i < correctedlength && edits <= maxEdits; i++){
                    if(correctedSequenceString[i] != task.decodedAnchor[i]){
                        tmp.edits.emplace_back(i, correctedSequenceString[i]);
                        edits++;
                    }
                }
                tmp.useEdits = edits <= maxEdits;
            }else{
                tmp.useEdits = false;
            }
            
            tmp.hq = task.msaProperties.isHQ;
            tmp.type = TempCorrectedSequenceType::Anchor;
            tmp.readId = task.input.anchorReadId;
            tmp.sequence = std::move(correctedSequenceString); 
            
            result.anchorCorrection = std::move(tmp);
        }
        
        
        
        for(const auto& correctedCandidate : task.candidateCorrections){
            const read_number candidateId = task.candidateReadIds[correctedCandidate.index];
            
            //don't save candidate if hq anchor correction exists
            bool savingIsOk = !correctionFlags->isCorrectedAsHQAnchor(candidateId);
                        
            if (savingIsOk) {                            
                
                TempCorrectedSequence tmp;
                
                tmp.type = TempCorrectedSequenceType::Candidate;
                tmp.readId = candidateId;
                tmp.shift = correctedCandidate.shift;

                const bool candidateIsForward = task.alignmentFlags[correctedCandidate.index] == AlignmentOrientation::Forward;

                if(candidateIsForward){
                    tmp.sequence = std::move(correctedCandidate.sequence);
                }else{
                    //if input candidate for correction is reverse complement, corrected candidate is also reverse complement
                    //get forward sequence
                    std::string fwd;
                    fwd.resize(correctedCandidate.sequence.length());
                    SequenceHelpers::reverseComplementSequenceDecoded(
                        &fwd[0], 
                        correctedCandidate.sequence.c_str(), 
                                            correctedCandidate.sequence.length()
                    );
                    tmp.sequence = std::move(fwd);
                }

                // if(candidateId == 1180257){
                //     std::cerr << "processing 1180257\n";
                //     std::cerr << "tmp.sequence = " << tmp.sequence << "\n";
                // }
                
                bool originalCandidateReadContainsN = false;
                readStorage->areSequencesAmbiguous(&originalCandidateReadContainsN, &candidateId, 1);
                
                if(!originalCandidateReadContainsN){
                    const std::size_t offset = correctedCandidate.index * decodedSequencePitchInBytes;
                    const char* const uncorrectedCandidate = &task.decodedCandidateSequences[offset];
                    const int uncorrectedCandidateLength = task.candidateSequencesLengths[correctedCandidate.index];
                    const int correctedCandidateLength = tmp.sequence.length();
                    
                    assert(uncorrectedCandidateLength == correctedCandidateLength);
                    
                    const int maxEdits = correctedCandidateLength / 7;
                    int edits = 0;
                    if(candidateIsForward){
                        for(int pos = 0; pos < correctedCandidateLength && edits <= maxEdits; pos++){
                            if(tmp.sequence[pos] != uncorrectedCandidate[pos]){
                                tmp.edits.emplace_back(pos, tmp.sequence[pos]);
                                edits++;
                            }
                        }
                    }else{
                        //tmp.sequence is forward sequence, but uncorrectedCandidate is reverse complement
                        std::string fwduncorrected;
                        fwduncorrected.resize(uncorrectedCandidateLength);
                        SequenceHelpers::reverseComplementSequenceDecoded(
                            &fwduncorrected[0], 
                            uncorrectedCandidate, 
                            uncorrectedCandidateLength
                        );

                        for(int pos = 0; pos < correctedCandidateLength && edits <= maxEdits; pos++){
                            if(tmp.sequence[pos] != fwduncorrected[pos]){
                                tmp.edits.emplace_back(pos, tmp.sequence[pos]);
                                edits++;
                            }
                        }
                    }
                    

                    tmp.useEdits = edits <= maxEdits;
                }else{
                    tmp.useEdits = false;
                }
                
                result.candidateCorrections.emplace_back(std::move(tmp));
            }
        }

        return result;
    }


private:

    std::size_t encodedSequencePitchInInts{};
    std::size_t decodedSequencePitchInBytes{};
    std::size_t qualityPitchInBytes{};

    const ProgramOptions* programOptions{};
    const CpuMinhasher* minhasher{};
    mutable MinhasherHandle minhashHandle;
    const CpuReadStorage* readStorage{};

    ReadCorrectionFlags* correctionFlags{};
    ClfAgent* clfAgent{};

    mutable cpu::shd::CpuAlignmentHandle alignmentHandle{};

    mutable std::stringstream ml_stream_anchor;
    mutable std::stringstream ml_stream_cands;

    std::unique_ptr<cpu::QualityScoreConversion> qualityCoversion;

    TimeMeasurements totalTime{};
};



}



#endif
