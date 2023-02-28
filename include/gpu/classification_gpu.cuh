#ifndef CARE_CLASSIFICATION_GPU_CUH
#define CARE_CLASSIFICATION_GPU_CUH

#include <classification.hpp>
#include <gpu/gpumsa.cuh>
#include <type_traits>
namespace care{

namespace gpu{

    struct ExtractAnchorInputData{
        char origBase;
        char consensusBase;
        float estimatedCoverage;
        int msaPos;
        int anchorColumnsBegin_incl;
        int anchorColumnsEnd_excl;
        GpuMSAProperties msaProperties;
        GpuSingleMSA msa;
    };

    struct ExtractCandidateInputData{
        char origBase;
        char consensusBase;
        float estimatedCoverage;
        int msaPos;
        int anchorColumnsBegin_incl;
        int anchorColumnsEnd_excl;
        int queryColumnsBegin_incl;
        int queryColumnsEnd_excl;
        GpuMSAProperties msaProperties;
        GpuSingleMSA msa;
    };

    namespace detail{

        struct extract_anchor_transformed{

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 37;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractAnchorInputData& input) const noexcept{
                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];

                *(features + 0) = float(input.origBase == 'A');
                *(features + 1) = float(input.origBase == 'C');
                *(features + 2) = float(input.origBase == 'G');
                *(features + 3) = float(input.origBase == 'T');
                *(features + 4) = float(input.consensusBase == 'A');
                *(features + 5) = float(input.consensusBase == 'C');
                *(features + 6) = float(input.consensusBase == 'G');
                *(features + 7) = float(input.consensusBase == 'T');
                *(features + 8) = input.origBase == 'A' ? countsA[msaPos] / countsACGT : 0;
                *(features + 9) = input.origBase == 'C' ? countsC[msaPos] / countsACGT : 0;
                *(features + 10) = input.origBase == 'G' ? countsG[msaPos] / countsACGT : 0;
                *(features + 11) = input.origBase == 'T' ? countsT[msaPos] / countsACGT : 0;
                *(features + 12) = input.origBase == 'A' ? weightsA[msaPos]:0;
                *(features + 13) = input.origBase == 'C' ? weightsC[msaPos]:0;
                *(features + 14) = input.origBase == 'G' ? weightsG[msaPos]:0;
                *(features + 15) = input.origBase == 'T' ? weightsT[msaPos]:0;
                *(features + 16) = input.consensusBase == 'A' ? countsA[msaPos] / countsACGT : 0;
                *(features + 17) = input.consensusBase == 'C' ? countsC[msaPos] / countsACGT : 0;
                *(features + 18) = input.consensusBase == 'G' ? countsG[msaPos] / countsACGT : 0;
                *(features + 19) = input.consensusBase == 'T' ? countsT[msaPos] / countsACGT : 0;
                *(features + 20) = input.consensusBase == 'A' ? weightsA[msaPos]:0;
                *(features + 21) = input.consensusBase == 'C' ? weightsC[msaPos]:0;
                *(features + 22) = input.consensusBase == 'G' ? weightsG[msaPos]:0;
                *(features + 23) = input.consensusBase == 'T' ? weightsT[msaPos]:0;
                *(features + 24) = weightsA[msaPos];
                *(features + 25) = weightsC[msaPos];
                *(features + 26) = weightsG[msaPos];
                *(features + 27) = weightsT[msaPos];
                *(features + 28) = countsA[msaPos] / countsACGT;
                *(features + 29) = countsC[msaPos] / countsACGT;
                *(features + 30) = countsG[msaPos] / countsACGT;
                *(features + 31) = countsT[msaPos] / countsACGT;
                *(features + 32) = input.msaProperties.avg_support;
                *(features + 33) = input.msaProperties.min_support;
                *(features + 34) = float(input.msaProperties.max_coverage) / input.estimatedCoverage;
                *(features + 35) = float(input.msaProperties.min_coverage) / input.estimatedCoverage;
                *(features + 36) = float(std::max(input.anchorColumnsBegin_incl - msaPos, msaPos - input.anchorColumnsEnd_excl)) / (input.anchorColumnsEnd_excl-input.anchorColumnsBegin_incl);
            }
        };

        struct extract_cands_transformed{

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 42;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractCandidateInputData& input) const noexcept{
                const int a_begin = input.anchorColumnsBegin_incl;
                const int a_end = input.anchorColumnsEnd_excl;
                const int c_begin = input.queryColumnsBegin_incl;
                const int c_end = input.queryColumnsEnd_excl;

                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];

                *(features + 0) = float(input.origBase == 'A');
                *(features + 1) = float(input.origBase == 'C');
                *(features + 2) = float(input.origBase == 'G');
                *(features + 3) = float(input.origBase == 'T');
                *(features + 4) = float(input.consensusBase == 'A');
                *(features + 5) = float(input.consensusBase == 'C');
                *(features + 6) = float(input.consensusBase == 'G');
                *(features + 7) = float(input.consensusBase == 'T');
                *(features + 8) = input.origBase == 'A'? countsA[msaPos] / countsACGT : 0;
                *(features + 9) = input.origBase == 'C'? countsC[msaPos] / countsACGT : 0;
                *(features + 10) = input.origBase == 'G'? countsG[msaPos] / countsACGT : 0;
                *(features + 11) = input.origBase == 'T'? countsT[msaPos] / countsACGT : 0;
                *(features + 12) = input.origBase == 'A'? weightsA[msaPos]:0;
                *(features + 13) = input.origBase == 'C'? weightsC[msaPos]:0;
                *(features + 14) = input.origBase == 'G'? weightsG[msaPos]:0;
                *(features + 15) = input.origBase == 'T'? weightsT[msaPos]:0;
                *(features + 16) = input.consensusBase == 'A'? countsA[msaPos] / countsACGT : 0;
                *(features + 17) = input.consensusBase == 'C'? countsC[msaPos] / countsACGT : 0;
                *(features + 18) = input.consensusBase == 'G'? countsG[msaPos] / countsACGT : 0;
                *(features + 19) = input.consensusBase == 'T'? countsT[msaPos] / countsACGT : 0;
                *(features + 20) = input.consensusBase == 'A'? weightsA[msaPos]:0;
                *(features + 21) = input.consensusBase == 'C'? weightsC[msaPos]:0;
                *(features + 22) = input.consensusBase == 'G'? weightsG[msaPos]:0;
                *(features + 23) = input.consensusBase == 'T'? weightsT[msaPos]:0;
                *(features + 24) = weightsA[msaPos];
                *(features + 25) = weightsC[msaPos];
                *(features + 26) = weightsG[msaPos];
                *(features + 27) = weightsT[msaPos];
                *(features + 28) = countsA[msaPos] / countsACGT;
                *(features + 29) = countsC[msaPos] / countsACGT;
                *(features + 30) = countsG[msaPos] / countsACGT;
                *(features + 31) = countsT[msaPos] / countsACGT;
                *(features + 32) = input.msaProperties.avg_support;
                *(features + 33) = input.msaProperties.min_support;
                *(features + 34) = float(input.msaProperties.max_coverage)/input.estimatedCoverage;
                *(features + 35) = float(input.msaProperties.min_coverage)/input.estimatedCoverage;
                *(features + 36) = float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(c_end-c_begin); // absolute shift (compatible with differing read lengths)
                *(features + 37) = float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(a_end-a_begin);
                *(features + 38) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(a_end-a_begin); // relative overlap (ratio of a or c length in case of diff. read len)
                *(features + 39) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(c_end-c_begin);
                *(features + 40) = float(std::max(a_begin-msaPos, msaPos-a_end))/(a_end-a_begin);
                *(features + 41) = float(std::max(a_begin-msaPos, msaPos-a_end))/(c_end-c_begin);
            }
        };

        struct extract_anchor {

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 21;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractAnchorInputData& input) const noexcept{
                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];

                *(features + 0) = float(input.origBase == 'A');
                *(features + 1) = float(input.origBase == 'C');
                *(features + 2) = float(input.origBase == 'G');
                *(features + 3) = float(input.origBase == 'T');
                *(features + 4) = float(input.consensusBase == 'A');
                *(features + 5) = float(input.consensusBase == 'C');
                *(features + 6) = float(input.consensusBase == 'G');
                *(features + 7) = float(input.consensusBase == 'T');
                *(features + 8) = weightsA[msaPos];
                *(features + 9) = weightsC[msaPos];
                *(features + 10) = weightsG[msaPos];
                *(features + 11) = weightsT[msaPos];
                *(features + 12) = countsA[msaPos]/countsACGT;
                *(features + 13) = countsC[msaPos]/countsACGT;
                *(features + 14) = countsG[msaPos]/countsACGT;
                *(features + 15) = countsT[msaPos]/countsACGT;
                *(features + 16) = input.msaProperties.avg_support;
                *(features + 17) = input.msaProperties.min_support;
                *(features + 18) = float(input.msaProperties.max_coverage)/input.estimatedCoverage;
                *(features + 19) = float(input.msaProperties.min_coverage)/input.estimatedCoverage;
                *(features + 20) = float(std::max(input.anchorColumnsBegin_incl-msaPos, msaPos-input.anchorColumnsEnd_excl))/(input.anchorColumnsEnd_excl-input.anchorColumnsBegin_incl);
            }
        };

        struct extract_cands {

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 26;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractCandidateInputData& input) const noexcept{  
                const int a_begin = input.anchorColumnsBegin_incl;
                const int a_end = input.anchorColumnsEnd_excl;
                const int c_begin = input.queryColumnsBegin_incl;
                const int c_end = input.queryColumnsEnd_excl;

                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];
                
                *(features + 0) = float(input.origBase == 'A');
                *(features + 1) = float(input.origBase == 'C');
                *(features + 2) = float(input.origBase == 'G');
                *(features + 3) = float(input.origBase == 'T');
                *(features + 4) = float(input.consensusBase == 'A');
                *(features + 5) = float(input.consensusBase == 'C');
                *(features + 6) = float(input.consensusBase == 'G');
                *(features + 7) = float(input.consensusBase == 'T');
                *(features + 8) = weightsA[msaPos];
                *(features + 9) = weightsC[msaPos];
                *(features + 10) = weightsG[msaPos];
                *(features + 11) = weightsT[msaPos];
                *(features + 12) = countsA[msaPos]/countsACGT;
                *(features + 13) = countsC[msaPos]/countsACGT;
                *(features + 14) = countsG[msaPos]/countsACGT;
                *(features + 15) = countsT[msaPos]/countsACGT;
                *(features + 16) = input.msaProperties.avg_support;
                *(features + 17) = input.msaProperties.min_support;
                *(features + 18) = float(input.msaProperties.max_coverage)/input.estimatedCoverage;
                *(features + 19) = float(input.msaProperties.min_coverage)/input.estimatedCoverage;
                *(features + 20) = float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(c_end-c_begin); // absolute shift (compatible with differing read lengths)
                *(features + 21) = float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(a_end-a_begin);
                *(features + 22) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(a_end-a_begin); // relative overlap (ratio of a or c length in case of diff. read len)
                *(features + 23) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(c_end-c_begin);
                *(features + 24) = float(std::max(a_begin-msaPos, msaPos-a_end))/(a_end-a_begin);
                *(features + 25) = float(std::max(a_begin-msaPos, msaPos-a_end))/(c_end-c_begin);
            }
        };

        struct extract_anchor_normed_weights {
            
            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 21;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractAnchorInputData& input) const noexcept{  
                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];
                const float weightsACGT = weightsA[msaPos] + weightsC[msaPos] + weightsG[msaPos] + weightsT[msaPos];

                const int a_begin = input.anchorColumnsBegin_incl;
                const int a_end = input.anchorColumnsEnd_excl;

                *(features + 0) = float(input.origBase == 'A');
                *(features + 1) = float(input.origBase == 'C');
                *(features + 2) = float(input.origBase == 'G');
                *(features + 3) = float(input.origBase == 'T');
                *(features + 4) = float(input.consensusBase == 'A');
                *(features + 5) = float(input.consensusBase == 'C');
                *(features + 6) = float(input.consensusBase == 'G');
                *(features + 7) = float(input.consensusBase == 'T');
                *(features + 8) = weightsA[msaPos]/weightsACGT;
                *(features + 9) = weightsC[msaPos]/weightsACGT;
                *(features + 10) = weightsG[msaPos]/weightsACGT;
                *(features + 11) = weightsT[msaPos]/weightsACGT;
                *(features + 12) = countsA[msaPos]/countsACGT;
                *(features + 13) = countsC[msaPos]/countsACGT;
                *(features + 14) = countsG[msaPos]/countsACGT;
                *(features + 15) = countsT[msaPos]/countsACGT;
                *(features + 16) = input.msaProperties.avg_support;
                *(features + 17) = input.msaProperties.min_support;
                *(features + 18) = float(input.msaProperties.max_coverage)/input.estimatedCoverage;
                *(features + 19) = float(input.msaProperties.min_coverage)/input.estimatedCoverage;
                *(features + 20) = float(std::max(a_begin-msaPos, msaPos-a_end))/(a_end-a_begin);
            }
        };

        struct extract_cands_normed_weights {
            
            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 26;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractCandidateInputData& input) const noexcept{  
                const int a_begin = input.anchorColumnsBegin_incl;
                const int a_end = input.anchorColumnsEnd_excl;
                const int c_begin = input.queryColumnsBegin_incl;
                const int c_end = input.queryColumnsEnd_excl;

                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];
                const float weightsACGT = weightsA[msaPos] + weightsC[msaPos] + weightsG[msaPos] + weightsT[msaPos];

                *(features + 0) = float(input.origBase == 'A');
                *(features + 1) = float(input.origBase == 'C');
                *(features + 2) = float(input.origBase == 'G');
                *(features + 3) = float(input.origBase == 'T');
                *(features + 4) = float(input.consensusBase == 'A');
                *(features + 5) = float(input.consensusBase == 'C');
                *(features + 6) = float(input.consensusBase == 'G');
                *(features + 7) = float(input.consensusBase == 'T');
                *(features + 8) = weightsA[msaPos]/weightsACGT;
                *(features + 9) = weightsC[msaPos]/weightsACGT;
                *(features + 10) = weightsG[msaPos]/weightsACGT;
                *(features + 11) = weightsT[msaPos]/weightsACGT;
                *(features + 12) = countsA[msaPos]/countsACGT;
                *(features + 13) = countsC[msaPos]/countsACGT;
                *(features + 14) = countsG[msaPos]/countsACGT;
                *(features + 15) = countsT[msaPos]/countsACGT;
                *(features + 16) = input.msaProperties.avg_support;
                *(features + 17) = input.msaProperties.min_support;
                *(features + 18) = float(input.msaProperties.max_coverage)/input.estimatedCoverage;
                *(features + 19) = float(input.msaProperties.min_coverage)/input.estimatedCoverage;
                *(features + 20) = float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(c_end-c_begin); // absolute shift (compatible with differing read lengths)
                *(features + 21) = float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(a_end-a_begin);
                *(features + 22) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(a_end-a_begin); // relative overlap (ratio of a or c length in case of diff. read len)
                *(features + 23) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(c_end-c_begin);
                *(features + 24) = float(std::max(a_begin-msaPos, msaPos-a_end))/(a_end-a_begin);
                *(features + 25) = float(std::max(a_begin-msaPos, msaPos-a_end))/(c_end-c_begin);
            }
        };

        struct extract_anchor_transformed_normed_weights{

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 37;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractAnchorInputData& input) const noexcept{
                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];
                const float weightsACGT = weightsA[msaPos] + weightsC[msaPos] + weightsG[msaPos] + weightsT[msaPos];

                *(features + 0) = float(input.origBase == 'A');
                *(features + 1) = float(input.origBase == 'C');
                *(features + 2) = float(input.origBase == 'G');
                *(features + 3) = float(input.origBase == 'T');
                *(features + 4) = float(input.consensusBase == 'A');
                *(features + 5) = float(input.consensusBase == 'C');
                *(features + 6) = float(input.consensusBase == 'G');
                *(features + 7) = float(input.consensusBase == 'T');
                *(features + 8) = input.origBase == 'A' ? countsA[msaPos] / countsACGT : 0;
                *(features + 9) = input.origBase == 'C' ? countsC[msaPos] / countsACGT : 0;
                *(features + 10) = input.origBase == 'G' ? countsG[msaPos] / countsACGT : 0;
                *(features + 11) = input.origBase == 'T' ? countsT[msaPos] / countsACGT : 0;
                *(features + 12) = input.origBase == 'A' ? weightsA[msaPos]/weightsACGT:0;
                *(features + 13) = input.origBase == 'C' ? weightsC[msaPos]/weightsACGT:0;
                *(features + 14) = input.origBase == 'G' ? weightsG[msaPos]/weightsACGT:0;
                *(features + 15) = input.origBase == 'T' ? weightsT[msaPos]/weightsACGT:0;
                *(features + 16) = input.consensusBase == 'A' ? countsA[msaPos] / countsACGT : 0;
                *(features + 17) = input.consensusBase == 'C' ? countsC[msaPos] / countsACGT : 0;
                *(features + 18) = input.consensusBase == 'G' ? countsG[msaPos] / countsACGT : 0;
                *(features + 19) = input.consensusBase == 'T' ? countsT[msaPos] / countsACGT : 0;
                *(features + 20) = input.consensusBase == 'A' ? weightsA[msaPos]/weightsACGT:0;
                *(features + 21) = input.consensusBase == 'C' ? weightsC[msaPos]/weightsACGT:0;
                *(features + 22) = input.consensusBase == 'G' ? weightsG[msaPos]/weightsACGT:0;
                *(features + 23) = input.consensusBase == 'T' ? weightsT[msaPos]/weightsACGT:0;
                *(features + 24) = weightsA[msaPos]/weightsACGT;
                *(features + 25) = weightsC[msaPos]/weightsACGT;
                *(features + 26) = weightsG[msaPos]/weightsACGT;
                *(features + 27) = weightsT[msaPos]/weightsACGT;
                *(features + 28) = countsA[msaPos] / countsACGT;
                *(features + 29) = countsC[msaPos] / countsACGT;
                *(features + 30) = countsG[msaPos] / countsACGT;
                *(features + 31) = countsT[msaPos] / countsACGT;
                *(features + 32) = input.msaProperties.avg_support;
                *(features + 33) = input.msaProperties.min_support;
                *(features + 34) = float(input.msaProperties.max_coverage) / input.estimatedCoverage;
                *(features + 35) = float(input.msaProperties.min_coverage) / input.estimatedCoverage;
                *(features + 36) = float(std::max(input.anchorColumnsBegin_incl - msaPos, msaPos - input.anchorColumnsEnd_excl)) / (input.anchorColumnsEnd_excl-input.anchorColumnsBegin_incl);
            }
        };

        struct extract_cands_transformed_normed_weights{

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 42;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractCandidateInputData& input) const noexcept{
                const int a_begin = input.anchorColumnsBegin_incl;
                const int a_end = input.anchorColumnsEnd_excl;
                const int c_begin = input.queryColumnsBegin_incl;
                const int c_end = input.queryColumnsEnd_excl;

                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];
                const float weightsACGT = weightsA[msaPos] + weightsC[msaPos] + weightsG[msaPos] + weightsT[msaPos];

                *(features + 0) = float(input.origBase == 'A');
                *(features + 1) = float(input.origBase == 'C');
                *(features + 2) = float(input.origBase == 'G');
                *(features + 3) = float(input.origBase == 'T');
                *(features + 4) = float(input.consensusBase == 'A');
                *(features + 5) = float(input.consensusBase == 'C');
                *(features + 6) = float(input.consensusBase == 'G');
                *(features + 7) = float(input.consensusBase == 'T');
                *(features + 8) = input.origBase == 'A'? countsA[msaPos] / countsACGT : 0;
                *(features + 9) = input.origBase == 'C'? countsC[msaPos] / countsACGT : 0;
                *(features + 10) = input.origBase == 'G'? countsG[msaPos] / countsACGT : 0;
                *(features + 11) = input.origBase == 'T'? countsT[msaPos] / countsACGT : 0;
                *(features + 12) = input.origBase == 'A'? weightsA[msaPos] / weightsACGT :0;
                *(features + 13) = input.origBase == 'C'? weightsC[msaPos] / weightsACGT :0;
                *(features + 14) = input.origBase == 'G'? weightsG[msaPos] / weightsACGT :0;
                *(features + 15) = input.origBase == 'T'? weightsT[msaPos] / weightsACGT :0;
                *(features + 16) = input.consensusBase == 'A'? countsA[msaPos] / countsACGT : 0;
                *(features + 17) = input.consensusBase == 'C'? countsC[msaPos] / countsACGT : 0;
                *(features + 18) = input.consensusBase == 'G'? countsG[msaPos] / countsACGT : 0;
                *(features + 19) = input.consensusBase == 'T'? countsT[msaPos] / countsACGT : 0;
                *(features + 20) = input.consensusBase == 'A'? weightsA[msaPos] / weightsACGT :0;
                *(features + 21) = input.consensusBase == 'C'? weightsC[msaPos] / weightsACGT :0;
                *(features + 22) = input.consensusBase == 'G'? weightsG[msaPos] / weightsACGT :0;
                *(features + 23) = input.consensusBase == 'T'? weightsT[msaPos] / weightsACGT :0;
                *(features + 24) = weightsA[msaPos] / weightsACGT ;
                *(features + 25) = weightsC[msaPos] / weightsACGT ;
                *(features + 26) = weightsG[msaPos] / weightsACGT ;
                *(features + 27) = weightsT[msaPos] / weightsACGT ;
                *(features + 28) = countsA[msaPos] / countsACGT;
                *(features + 29) = countsC[msaPos] / countsACGT;
                *(features + 30) = countsG[msaPos] / countsACGT;
                *(features + 31) = countsT[msaPos] / countsACGT;
                *(features + 32) = input.msaProperties.avg_support;
                *(features + 33) = input.msaProperties.min_support;
                *(features + 34) = float(input.msaProperties.max_coverage)/input.estimatedCoverage;
                *(features + 35) = float(input.msaProperties.min_coverage)/input.estimatedCoverage;
                *(features + 36) = float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(c_end-c_begin); // absolute shift (compatible with differing read lengths)
                *(features + 37) = float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(a_end-a_begin);
                *(features + 38) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(a_end-a_begin); // relative overlap (ratio of a or c length in case of diff. read len)
                *(features + 39) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(c_end-c_begin);
                *(features + 40) = float(std::max(a_begin-msaPos, msaPos-a_end))/(a_end-a_begin);
                *(features + 41) = float(std::max(a_begin-msaPos, msaPos-a_end))/(c_end-c_begin);
            }
        };


        struct extract_anchor_v2{

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 11;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractAnchorInputData& input) const noexcept{
                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];

                const float countsACGT = input.msa.coverages[msaPos];
                const float weightsACGT = weightsA[msaPos] + weightsC[msaPos] + weightsG[msaPos] + weightsT[msaPos];

                float weightCons = 0, countCons = 0;

                switch (input.consensusBase) {
                    case 'A':
                        weightCons = weightsA[msaPos];
                        countCons = countsA[msaPos];
                        break;
                    case 'C':
                        weightCons = weightsC[msaPos];
                        countCons = countsC[msaPos];
                        break;
                    case 'G':
                        weightCons = weightsG[msaPos];
                        countCons = countsG[msaPos];
                        break;
                    case 'T':
                        weightCons = weightsT[msaPos];
                        countCons = countsT[msaPos];
                        break;
                    default: break;
                }

                *(features + 0) = float(input.msa.origCoverages[msaPos]) / countsACGT;
                *(features + 1) = input.msa.origWeights[msaPos] / weightsACGT;
                *(features + 2) = input.msa.origWeights[msaPos] / float(input.msa.origCoverages[msaPos]);
                *(features + 3) = countCons / countsACGT;
                *(features + 4) = input.msa.support[msaPos];
                *(features + 5) = weightCons / countCons;
                *(features + 6) = weightsACGT / countsACGT;
                *(features + 7) = input.msaProperties.avg_support;
                *(features + 8) = input.msaProperties.min_support;
                *(features + 9) = float(input.msaProperties.max_coverage) / input.estimatedCoverage;
                *(features + 10) = float(input.msaProperties.min_coverage) / input.estimatedCoverage;
             }
        };

        struct extract_cands_v2{

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 12;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractCandidateInputData& input) const noexcept{
                const int a_begin = input.anchorColumnsBegin_incl;
                const int a_end = input.anchorColumnsEnd_excl;
                const int c_begin = input.queryColumnsBegin_incl;
                const int c_end = input.queryColumnsEnd_excl;

                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];
                const float weightsACGT = weightsA[msaPos] + weightsC[msaPos] + weightsG[msaPos] + weightsT[msaPos];

                float weightCons = 0, countCons = 0, weightOrig = 0, countOrig = 0;

                switch (input.consensusBase) {
                    case 'A':
                        weightCons = weightsA[msaPos];
                        countCons = countsA[msaPos];
                        break;
                    case 'C':
                        weightCons = weightsC[msaPos];
                        countCons = countsC[msaPos];
                        break;
                    case 'G':
                        weightCons = weightsG[msaPos];
                        countCons = countsG[msaPos];
                        break;
                    case 'T':
                        weightCons = weightsT[msaPos];
                        countCons = countsT[msaPos];
                        break;
                    default: break;
                }

                switch (input.origBase) {
                    case 'A':
                        weightOrig = weightsA[msaPos];
                        countOrig = countsA[msaPos];
                        break;
                    case 'C':
                        weightOrig = weightsC[msaPos];
                        countOrig = countsC[msaPos];
                        break;
                    case 'G':
                        weightOrig = weightsG[msaPos];
                        countOrig = countsG[msaPos];
                        break;
                    case 'T':
                        weightOrig = weightsT[msaPos];
                        countOrig = countsT[msaPos];
                        break;
                    default: break;
                }

                *(features + 0) = countOrig / countsACGT;
                *(features + 1) = weightOrig / weightsACGT;
                *(features + 2) = weightOrig / countOrig;
                *(features + 3) = countCons / countsACGT;
                *(features + 4) = input.msa.support[msaPos];
                *(features + 5) = weightCons / countCons;
                *(features + 6) = weightsACGT / countsACGT;
                *(features + 7) = input.msaProperties.avg_support;
                *(features + 8) = input.msaProperties.min_support;
                *(features + 9) = float(input.msaProperties.max_coverage) / input.estimatedCoverage;
                *(features + 10) = float(input.msaProperties.min_coverage) / input.estimatedCoverage;
                *(features + 11) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))
                    / (std::max(a_end, c_end)-std::min(a_begin, c_begin)); // jaccard
            }
        };

        struct extract_anchor_v3 {

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 13;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractAnchorInputData& input) const noexcept{
                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];

                const float countsACGT = input.msa.coverages[msaPos];
                const float weightsACGT = weightsA[msaPos] + weightsC[msaPos] + weightsG[msaPos] + weightsT[msaPos];

                float weightCons = 0, countCons = 0;

                switch (input.consensusBase) {
                    case 'A':
                        weightCons = weightsA[msaPos];
                        countCons = countsA[msaPos];
                        break;
                    case 'C':
                        weightCons = weightsC[msaPos];
                        countCons = countsC[msaPos];
                        break;
                    case 'G':
                        weightCons = weightsG[msaPos];
                        countCons = countsG[msaPos];
                        break;
                    case 'T':
                        weightCons = weightsT[msaPos];
                        countCons = countsT[msaPos];
                        break;
                    default: break;
                }

                *(features + 0) = float(input.msa.origCoverages[msaPos]) / countsACGT;
                *(features + 1) = input.msa.origWeights[msaPos] / weightsACGT;
                *(features + 2) = input.msa.origWeights[msaPos] / float(input.msa.origCoverages[msaPos]);
                *(features + 3) = countCons / countsACGT;
                *(features + 4) = input.msa.support[msaPos];
                *(features + 5) = weightCons / countCons;
                *(features + 6) = weightsACGT / countsACGT;
                *(features + 7) = countsACGT / input.estimatedCoverage;
                *(features + 8) = weightsACGT / input.estimatedCoverage;
                *(features + 9) = input.msaProperties.avg_support;
                *(features + 10) = input.msaProperties.min_support;
                *(features + 11) = float(input.msaProperties.max_coverage) / input.estimatedCoverage;
                *(features + 12) = float(input.msaProperties.min_coverage) / input.estimatedCoverage;
             }
        };

        struct extract_cands_v3 {

            HOSTDEVICEQUALIFIER
            static constexpr int numFeatures() noexcept{
                return 14;
            }

            template<class OutIter>
            HOSTDEVICEQUALIFIER
            void operator()(OutIter features, const ExtractCandidateInputData& input) const noexcept{
                const int a_begin = input.anchorColumnsBegin_incl;
                const int a_end = input.anchorColumnsEnd_excl;
                const int c_begin = input.queryColumnsBegin_incl;
                const int c_end = input.queryColumnsEnd_excl;

                const int msaPos = input.msaPos;
                const int pitch = input.msa.columnPitchInElements;
                const float countsACGT = input.msa.coverages[msaPos];
                const int* const countsA = &input.msa.counts[0 * pitch];
                const int* const countsC = &input.msa.counts[1 * pitch];
                const int* const countsG = &input.msa.counts[2 * pitch];
                const int* const countsT = &input.msa.counts[3 * pitch];

                const float* const weightsA = &input.msa.weights[0 * pitch];
                const float* const weightsC = &input.msa.weights[1 * pitch];
                const float* const weightsG = &input.msa.weights[2 * pitch];
                const float* const weightsT = &input.msa.weights[3 * pitch];
                const float weightsACGT = weightsA[msaPos] + weightsC[msaPos] + weightsG[msaPos] + weightsT[msaPos];

                float weightCons = 0, countCons = 0, weightOrig = 0, countOrig = 0;

                switch (input.consensusBase) {
                    case 'A':
                        weightCons = weightsA[msaPos];
                        countCons = countsA[msaPos];
                        break;
                    case 'C':
                        weightCons = weightsC[msaPos];
                        countCons = countsC[msaPos];
                        break;
                    case 'G':
                        weightCons = weightsG[msaPos];
                        countCons = countsG[msaPos];
                        break;
                    case 'T':
                        weightCons = weightsT[msaPos];
                        countCons = countsT[msaPos];
                        break;
                    default: break;
                }

                switch (input.origBase) {
                    case 'A':
                        weightOrig = weightsA[msaPos];
                        countOrig = countsA[msaPos];
                        break;
                    case 'C':
                        weightOrig = weightsC[msaPos];
                        countOrig = countsC[msaPos];
                        break;
                    case 'G':
                        weightOrig = weightsG[msaPos];
                        countOrig = countsG[msaPos];
                        break;
                    case 'T':
                        weightOrig = weightsT[msaPos];
                        countOrig = countsT[msaPos];
                        break;
                    default: break;
                }

                *(features + 0) = countOrig / countsACGT;
                *(features + 1) = weightOrig / weightsACGT;
                *(features + 2) = weightOrig / countOrig;
                *(features + 3) = countCons / countsACGT;
                *(features + 4) = input.msa.support[msaPos];
                *(features + 5) = weightCons / countCons;
                *(features + 6) = weightsACGT / countsACGT;
                *(features + 7) = countsACGT / input.estimatedCoverage;
                *(features + 8) = weightsACGT / input.estimatedCoverage;
                *(features + 9) = input.msaProperties.avg_support;
                *(features + 10) = input.msaProperties.min_support;
                *(features + 11) = float(input.msaProperties.max_coverage) / input.estimatedCoverage;
                *(features + 12) = float(input.msaProperties.min_coverage) / input.estimatedCoverage;
                *(features + 13) = float(std::min(a_end, c_end)-std::max(a_begin, c_begin))
                    / (std::max(a_end, c_end)-std::min(a_begin, c_begin)); // jaccard
            }
        };


        //map cpu extractors to corresponding gpu extractors

        template<class CpuExtractor> struct GpuExtractorSelection;
        template<> struct GpuExtractorSelection<care::detail::extract_anchor>{using type = care::gpu::detail::extract_anchor;};
        template<> struct GpuExtractorSelection<care::detail::extract_cands>{using type = care::gpu::detail::extract_cands;};
        template<> struct GpuExtractorSelection<care::detail::extract_anchor_transformed>{using type = care::gpu::detail::extract_anchor_transformed;};
        template<> struct GpuExtractorSelection<care::detail::extract_cands_transformed>{using type = care::gpu::detail::extract_cands_transformed;};
        template<> struct GpuExtractorSelection<care::detail::extract_anchor_normed_weights>{using type = care::gpu::detail::extract_anchor_normed_weights;};
        template<> struct GpuExtractorSelection<care::detail::extract_cands_normed_weights>{using type = care::gpu::detail::extract_cands_normed_weights;};
        template<> struct GpuExtractorSelection<care::detail::extract_anchor_transformed_normed_weights>{using type = care::gpu::detail::extract_anchor_transformed_normed_weights;};
        template<> struct GpuExtractorSelection<care::detail::extract_cands_transformed_normed_weights>{using type = care::gpu::detail::extract_cands_transformed_normed_weights;};
        template<> struct GpuExtractorSelection<care::detail::extract_anchor_v2>{using type = care::gpu::detail::extract_anchor_v2;};
        template<> struct GpuExtractorSelection<care::detail::extract_cands_v2>{using type = care::gpu::detail::extract_cands_v2;};
        template<> struct GpuExtractorSelection<care::detail::extract_anchor_support>{using type = care::gpu::detail::extract_anchor_v2;};
        template<> struct GpuExtractorSelection<care::detail::extract_cands_support>{using type = care::gpu::detail::extract_cands_v2;};
        template<> struct GpuExtractorSelection<care::detail::extract_anchor_v3>{using type = care::gpu::detail::extract_anchor_v3;};
        template<> struct GpuExtractorSelection<care::detail::extract_cands_v3>{using type = care::gpu::detail::extract_cands_v3;};

    } //namespace detail

    using anchor_extractor = detail::GpuExtractorSelection<care::anchor_extractor>::type;
    using cands_extractor = detail::GpuExtractorSelection<care::cands_extractor>::type;

} //namespace gpu

} //namespace care



#endif