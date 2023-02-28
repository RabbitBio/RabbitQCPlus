#ifndef CARE_CPUCORRECTOR_TASK_HPP
#define CARE_CPUCORRECTOR_TASK_HPP

#include <config.hpp>

#include <cpu_alignment.hpp>
#include <alignmentorientation.hpp>
#include <msa.hpp>

#include <vector>

namespace care{

    struct CpuErrorCorrectorInput{
        int anchorLength{};
        read_number anchorReadId{};
        const unsigned int* encodedAnchor{};
        const char* anchorQualityscores{};
    };

    struct CpuErrorCorrectorMultiInput{
        std::vector<int> anchorLengths;
        std::vector<read_number> anchorReadIds;
        std::vector<const unsigned int*> encodedAnchors;
        std::vector<const char*> anchorQualityscores;
    };

    struct CpuErrorCorrectorTask{
        bool active{};

        std::vector<read_number> candidateReadIds{};
        std::vector<read_number> filteredReadIds{};
        std::vector<unsigned int> candidateSequencesData{};
        std::vector<unsigned int> candidateSequencesRevcData{};
        std::vector<int> candidateSequencesLengths{};
        std::vector<int> alignmentShifts{};
        std::vector<int> alignmentOps{};
        std::vector<int> alignmentOverlaps{};
        std::vector<float> alignmentWeights{};
        std::vector<char> candidateQualities{};
        std::vector<char> decodedAnchor{};
        std::vector<char> decodedCandidateSequences{};
        std::vector<cpu::SHDResult> alignments{};
        std::vector<cpu::SHDResult> revcAlignments{};
        std::vector<AlignmentOrientation> alignmentFlags{};
        std::vector<bool> isPairedCandidate{};

        CpuErrorCorrectorInput input{};

        CorrectionResult anchorCorrection;
        std::vector<CorrectedCandidate> candidateCorrections;
        MSAProperties msaProperties;
        MultipleSequenceAlignment multipleSequenceAlignment;
    };

}

#endif