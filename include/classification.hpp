#ifndef CARE_CLASSIFICATION_HPP
#define CARE_CLASSIFICATION_HPP

#include <array>
#include <random>
#include <forest.hpp>
#include <logreg.hpp>
#include <msa.hpp>
#include <cpucorrectortask.hpp>
#include <options.hpp>
#include <chrono>
#include <thread>

// This header allows toggling of feature transformations and classifiers,
// and seperates classification logic from the main corrector code.
// SEE BOTTOM OF FILE FOR TOGGLES! 

// TODO: implement logistic regression and investigate peformance
// Possibly same accuracy with VASTLY superior performance

// The current features were designed with the possibility of using logistic regression in mind
// but are highly redundant for any decision tree.
// However, since sklearn only supports float-type features, we might aswell do the one-hot-vector-encoding here for now.


namespace care {

struct ClfAgentDecisionInputData{
    const char* decodedAnchor{};
    int anchorColumnsBegin_incl{};
    int anchorColumnsEnd_excl{};
    const int* alignmentShifts{};
    const int* candidateSequencesLengths{};
    const char* decodedCandidateSequences{};
    const int* countsA{};
    const int* countsC{};
    const int* countsG{};
    const int* countsT{};
    const float* weightsA{};
    const float* weightsC{};
    const float* weightsG{};
    const float* weightsT{};
    const int* coverages{};
    const char* consensus{};

    std::function<MSAProperties(int,int,float,float,float)> getMSAProperties{};

    MSAProperties anchorMsaProperties{};
};

template<typename AnchorClf,
         typename CandClf,
         typename AnchorExtractor,
         typename CandsExtractor>
struct clf_agent
{

    static uint64_t get_seed() {
        return std::chrono::high_resolution_clock::now().time_since_epoch().count() + std::hash<std::thread::id>()(std::this_thread::get_id());
    }

    //TODO: access permission
    std::shared_ptr<AnchorClf> classifier_anchor;
    std::shared_ptr<CandClf> classifier_cands;
    std::stringstream anchor_stream, cands_stream;
    std::shared_ptr<std::ofstream> anchor_print_file, cands_print_file;
    std::ranlux48 rng;
    std::bernoulli_distribution coinflip_anchor, coinflip_cands;
    AnchorExtractor extract_anchor;
    CandsExtractor extract_cands;

    clf_agent(const ProgramOptions& opts) :
        classifier_anchor(opts.correctionType == CorrectionType::Forest ? std::make_shared<AnchorClf>(opts.mlForestfileAnchor, opts.maxForestTreesAnchor, opts.thresholdAnchor) : nullptr),
        classifier_cands(opts.correctionTypeCands == CorrectionType::Forest ? std::make_shared<CandClf>(opts.mlForestfileCands, opts.maxForestTreesCands, opts.thresholdCands) : nullptr),
        anchor_print_file(opts.correctionType == CorrectionType::Print ? std::make_shared<std::ofstream>(opts.mlForestfilePrintAnchor) : nullptr),
        cands_print_file(opts.correctionTypeCands == CorrectionType::Print ? std::make_shared<std::ofstream>(opts.mlForestfilePrintCands) : nullptr),
        rng(get_seed()),
        coinflip_anchor(opts.sampleRateAnchor),
        coinflip_cands(opts.sampleRateCands)
    {
        if (opts.correctionType == CorrectionType::Print) {
            *anchor_print_file << extract_anchor << std::endl;
        }

        if (opts.correctionTypeCands == CorrectionType::Print) {
            *cands_print_file << extract_cands << std::endl;
        }
    }

    clf_agent(const clf_agent& other) :
        classifier_anchor(other.classifier_anchor),
        classifier_cands(other.classifier_cands),
        anchor_print_file(other.anchor_print_file),
        cands_print_file(other.cands_print_file),
        rng(get_seed()),
        coinflip_anchor(other.coinflip_anchor),
        coinflip_cands(other.coinflip_cands)
    {}

    void print_anchor(const CpuErrorCorrectorTask& task, size_t i, const ProgramOptions& opt) {       
        if(!coinflip_anchor(rng)) return;

        anchor_stream << task.input.anchorReadId << ' ' << i << ' ' << task.multipleSequenceAlignment.consensus[task.multipleSequenceAlignment.anchorColumnsBegin_incl+i] << ' ';
        for (float j: extract_anchor(task, i, opt))
            anchor_stream << j << ' ';
        anchor_stream << '\n';
    }

    void print_cand(const CpuErrorCorrectorTask& task, int i, const ProgramOptions& opt, size_t cand, size_t offset) {       
        if(!coinflip_cands(rng)) return;

        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int c_begin = a_begin + task.alignmentShifts[cand];
        int pos = c_begin + i;

        cands_stream << task.candidateReadIds[cand] << ' ' << (task.alignmentFlags[cand]==AlignmentOrientation::ReverseComplement?-i-1:i) << ' ' << task.multipleSequenceAlignment.consensus[pos] << ' ';
        for (float j: extract_cands(task, i, opt, cand, offset))
            cands_stream << j << ' ';
        cands_stream << '\n';
    }

    template<typename... Args>
    bool decide_anchor(Args&&...args) {
        auto feats = extract_anchor(std::forward<Args>(args)...);
        // #pragma omp critical
        // {
        //     for (const auto& f: feats)
        //         std::cerr << f << ' ';
        //     std::cerr << "=>  " << classifier_anchor->prob_debug(feats) << std::endl;
        // }
        return classifier_anchor->decide(feats);
    }

    template<typename... Args>
    bool decide_cand(Args&&...args) {
        return classifier_cands->decide(extract_cands(std::forward<Args>(args)...));
    }

    void flush() {
        if (anchor_print_file && anchor_stream.peek() != decltype(anchor_stream)::traits_type::eof()) {
            #pragma omp critical
            {
                *anchor_print_file << anchor_stream.rdbuf();
            }
        }
        anchor_stream = std::stringstream();

        if (cands_print_file && cands_stream.peek() != decltype(cands_stream)::traits_type::eof()) {
            #pragma omp critical
            {
                *cands_print_file << cands_stream.rdbuf();
            }
        }
        cands_stream = std::stringstream();
    }
};

namespace detail {

struct extract_anchor {

    constexpr operator auto() {
        return u8"21 extract_anchor";
    }

    using features_t = std::array<float, 21>;

    features_t operator()(const CpuErrorCorrectorTask& task, int i, const ProgramOptions& opt) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int pos = a_begin + i;
        char orig = task.decodedAnchor[i];
        float countsACGT = msa.countsA[pos] + msa.countsC[pos] + msa.countsG[pos] + msa.countsT[pos];
        return {
            float(orig == 'A'),
            float(orig == 'C'),
            float(orig == 'G'),
            float(orig == 'T'),
            float(msa.consensus[pos] == 'A'),
            float(msa.consensus[pos] == 'C'),
            float(msa.consensus[pos] == 'G'),
            float(msa.consensus[pos] == 'T'),
            msa.weightsA[pos],
            msa.weightsC[pos],
            msa.weightsG[pos],
            msa.weightsT[pos],
            msa.countsA[pos]/countsACGT,
            msa.countsC[pos]/countsACGT,
            msa.countsG[pos]/countsACGT,
            msa.countsT[pos]/countsACGT,
            task.msaProperties.avg_support,
            task.msaProperties.min_support,
            float(task.msaProperties.max_coverage)/opt.estimatedCoverage,
            float(task.msaProperties.min_coverage)/opt.estimatedCoverage,
            float(std::max(a_begin-pos, pos-a_end))/(a_end-a_begin)
        };
    }
};

struct extract_cands {

    constexpr operator auto() {
        return u8"26 extract_cands";
    }

    using features_t = std::array<float, 26>;

    features_t operator()(const CpuErrorCorrectorTask& task, size_t i, const ProgramOptions& opt, size_t cand, size_t offset) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int c_begin = a_begin + task.alignmentShifts[cand];
        int c_end = c_begin + task.candidateSequencesLengths[cand];
        int pos = c_begin + i;
        char orig = task.decodedCandidateSequences[offset+i];
        float countsACGT = msa.countsA[pos] + msa.countsC[pos] + msa.countsG[pos] + msa.countsT[pos];
        MSAProperties props = msa.getMSAProperties(c_begin, c_end, opt.estimatedErrorrate, opt.estimatedCoverage, opt.m_coverage);
        return {
            float(orig == 'A'),
            float(orig == 'C'),
            float(orig == 'G'),
            float(orig == 'T'),
            float(msa.consensus[pos] == 'A'),
            float(msa.consensus[pos] == 'C'),
            float(msa.consensus[pos] == 'G'),
            float(msa.consensus[pos] == 'T'),
            msa.weightsA[pos],
            msa.weightsC[pos],
            msa.weightsG[pos],
            msa.weightsT[pos],
            msa.countsA[pos]/countsACGT,
            msa.countsC[pos]/countsACGT,
            msa.countsG[pos]/countsACGT,
            msa.countsT[pos]/countsACGT,
            props.avg_support,
            props.min_support,
            float(props.max_coverage)/opt.estimatedCoverage,
            float(props.min_coverage)/opt.estimatedCoverage,
            float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(c_end-c_begin), // absolute shift (compatible with differing read lengths)
            float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(a_end-a_begin),
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(a_end-a_begin), // relative overlap (ratio of a or c length in case of diff. read len)
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(c_end-c_begin),
            float(std::max(a_begin-pos, pos-a_end))/(a_end-a_begin),
            float(std::max(a_begin-pos, pos-a_end))/(c_end-c_begin)
        };
    }
};

struct extract_anchor_transformed {
    using features_t = std::array<float, 37>;

    constexpr operator auto() {
        return u8"37 extract_anchor_transformed";
    }

    features_t operator()(const CpuErrorCorrectorTask& task, int i, const ProgramOptions& opt) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int pos = a_begin + i;
        char orig = task.decodedAnchor[i];
        float countsACGT = msa.countsA[pos] + msa.countsC[pos] + msa.countsG[pos] + msa.countsT[pos];
        return {
            float(orig == 'A'),
            float(orig == 'C'),
            float(orig == 'G'),
            float(orig == 'T'),
            float(msa.consensus[pos] == 'A'),
            float(msa.consensus[pos] == 'C'),
            float(msa.consensus[pos] == 'G'),
            float(msa.consensus[pos] == 'T'),
            orig == 'A'?msa.countsA[pos]/countsACGT:0,
            orig == 'C'?msa.countsC[pos]/countsACGT:0,
            orig == 'G'?msa.countsG[pos]/countsACGT:0,
            orig == 'T'?msa.countsT[pos]/countsACGT:0,
            orig == 'A'?msa.weightsA[pos]:0,
            orig == 'C'?msa.weightsC[pos]:0,
            orig == 'G'?msa.weightsG[pos]:0,
            orig == 'T'?msa.weightsT[pos]:0,
            msa.consensus[pos] == 'A'?msa.countsA[pos]/countsACGT:0,
            msa.consensus[pos] == 'C'?msa.countsC[pos]/countsACGT:0,
            msa.consensus[pos] == 'G'?msa.countsG[pos]/countsACGT:0,
            msa.consensus[pos] == 'T'?msa.countsT[pos]/countsACGT:0,
            msa.consensus[pos] == 'A'?msa.weightsA[pos]:0,
            msa.consensus[pos] == 'C'?msa.weightsC[pos]:0,
            msa.consensus[pos] == 'G'?msa.weightsG[pos]:0,
            msa.consensus[pos] == 'T'?msa.weightsT[pos]:0,
            msa.weightsA[pos],
            msa.weightsC[pos],
            msa.weightsG[pos],
            msa.weightsT[pos],
            msa.countsA[pos]/countsACGT,
            msa.countsC[pos]/countsACGT,
            msa.countsG[pos]/countsACGT,
            msa.countsT[pos]/countsACGT,
            task.msaProperties.avg_support,
            task.msaProperties.min_support,
            float(task.msaProperties.max_coverage)/opt.estimatedCoverage,
            float(task.msaProperties.min_coverage)/opt.estimatedCoverage,
            float(std::max(a_begin-pos, pos-a_end))/(a_end-a_begin)
        };
    }
};

struct extract_cands_transformed {
    using features_t = std::array<float, 42>;

    constexpr operator auto() {
        return u8"42 extract_cands_transformed";
    }

    features_t operator()(const CpuErrorCorrectorTask& task, size_t i, const ProgramOptions& opt, size_t cand, size_t offset) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int c_begin = a_begin + task.alignmentShifts[cand];
        int c_end = c_begin + task.candidateSequencesLengths[cand];
        int pos = c_begin + i;
        char orig = task.decodedCandidateSequences[offset+i];
        float countsACGT = msa.countsA[pos] + msa.countsC[pos] + msa.countsG[pos] + msa.countsT[pos];
        MSAProperties props = msa.getMSAProperties(c_begin, c_end, opt.estimatedErrorrate, opt.estimatedCoverage, opt.m_coverage);
        return {
            float(orig == 'A'),
            float(orig == 'C'),
            float(orig == 'G'),
            float(orig == 'T'),
            float(msa.consensus[pos] == 'A'),
            float(msa.consensus[pos] == 'C'),
            float(msa.consensus[pos] == 'G'),
            float(msa.consensus[pos] == 'T'),
            orig == 'A'?msa.countsA[pos]/countsACGT:0,
            orig == 'C'?msa.countsC[pos]/countsACGT:0,
            orig == 'G'?msa.countsG[pos]/countsACGT:0,
            orig == 'T'?msa.countsT[pos]/countsACGT:0,
            orig == 'A'?msa.weightsA[pos]:0,
            orig == 'C'?msa.weightsC[pos]:0,
            orig == 'G'?msa.weightsG[pos]:0,
            orig == 'T'?msa.weightsT[pos]:0,
            msa.consensus[pos] == 'A'?msa.countsA[pos]/countsACGT:0,
            msa.consensus[pos] == 'C'?msa.countsC[pos]/countsACGT:0,
            msa.consensus[pos] == 'G'?msa.countsG[pos]/countsACGT:0,
            msa.consensus[pos] == 'T'?msa.countsT[pos]/countsACGT:0,
            msa.consensus[pos] == 'A'?msa.weightsA[pos]:0,
            msa.consensus[pos] == 'C'?msa.weightsC[pos]:0,
            msa.consensus[pos] == 'G'?msa.weightsG[pos]:0,
            msa.consensus[pos] == 'T'?msa.weightsT[pos]:0,
            msa.weightsA[pos],
            msa.weightsC[pos],
            msa.weightsG[pos],
            msa.weightsT[pos],
            msa.countsA[pos]/countsACGT,
            msa.countsC[pos]/countsACGT,
            msa.countsG[pos]/countsACGT,
            msa.countsT[pos]/countsACGT,
            props.avg_support,
            props.min_support,
            float(props.max_coverage)/opt.estimatedCoverage,
            float(props.min_coverage)/opt.estimatedCoverage,
            float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(c_end-c_begin), // absolute shift (compatible with differing read lengths)
            float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(a_end-a_begin),
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(a_end-a_begin), // relative overlap (ratio of a or c length in case of diff. read len)
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(c_end-c_begin),
            float(std::max(a_begin-pos, pos-a_end))/(a_end-a_begin),
            float(std::max(a_begin-pos, pos-a_end))/(c_end-c_begin)
        };
    }
};

struct extract_anchor_normed_weights {
    using features_t = std::array<float, 21>;

    constexpr operator auto() {
        return u8"21 extract_anchor_normed_weights";
    }

    features_t operator()(const CpuErrorCorrectorTask& task, int i, const ProgramOptions& opt) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int pos = a_begin + i;
        char orig = task.decodedAnchor[i];
        float countsACGT = msa.countsA[pos] + msa.countsC[pos] + msa.countsG[pos] + msa.countsT[pos];
        float weightsACGT = msa.weightsA[pos] + msa.weightsC[pos] + msa.weightsG[pos] + msa.weightsT[pos];
        return {
            float(orig == 'A'),
            float(orig == 'C'),
            float(orig == 'G'),
            float(orig == 'T'),
            float(msa.consensus[pos] == 'A'),
            float(msa.consensus[pos] == 'C'),
            float(msa.consensus[pos] == 'G'),
            float(msa.consensus[pos] == 'T'),
            msa.weightsA[pos]/weightsACGT,
            msa.weightsC[pos]/weightsACGT,
            msa.weightsG[pos]/weightsACGT,
            msa.weightsT[pos]/weightsACGT,
            msa.countsA[pos]/countsACGT,
            msa.countsC[pos]/countsACGT,
            msa.countsG[pos]/countsACGT,
            msa.countsT[pos]/countsACGT,
            task.msaProperties.avg_support,
            task.msaProperties.min_support,
            float(task.msaProperties.max_coverage)/opt.estimatedCoverage,
            float(task.msaProperties.min_coverage)/opt.estimatedCoverage,
            float(std::max(a_begin-pos, pos-a_end))/(a_end-a_begin)
        };
    }
};

struct extract_cands_normed_weights {
    using features_t = std::array<float, 26>;

    constexpr operator auto() {
        return u8"26 extract_cands_normed_weights";
    }

    features_t operator()(const CpuErrorCorrectorTask& task, size_t i, const ProgramOptions& opt, size_t cand, size_t offset) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int c_begin = a_begin + task.alignmentShifts[cand];
        int c_end = c_begin + task.candidateSequencesLengths[cand];
        int pos = c_begin + i;
        char orig = task.decodedCandidateSequences[offset+i];
        float countsACGT = msa.countsA[pos] + msa.countsC[pos] + msa.countsG[pos] + msa.countsT[pos];
        float weightsACGT = msa.weightsA[pos] + msa.weightsC[pos] + msa.weightsG[pos] + msa.weightsT[pos];
        MSAProperties props = msa.getMSAProperties(c_begin, c_end, opt.estimatedErrorrate, opt.estimatedCoverage, opt.m_coverage);
        return {
            float(orig == 'A'),
            float(orig == 'C'),
            float(orig == 'G'),
            float(orig == 'T'),
            float(msa.consensus[pos] == 'A'),
            float(msa.consensus[pos] == 'C'),
            float(msa.consensus[pos] == 'G'),
            float(msa.consensus[pos] == 'T'),
            msa.weightsA[pos]/weightsACGT,
            msa.weightsC[pos]/weightsACGT,
            msa.weightsG[pos]/weightsACGT,
            msa.weightsT[pos]/weightsACGT,
            msa.countsA[pos]/countsACGT,
            msa.countsC[pos]/countsACGT,
            msa.countsG[pos]/countsACGT,
            msa.countsT[pos]/countsACGT,
            props.avg_support,
            props.min_support,
            float(props.max_coverage)/opt.estimatedCoverage,
            float(props.min_coverage)/opt.estimatedCoverage,
            float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(c_end-c_begin), // absolute shift (compatible with differing read lengths)
            float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(a_end-a_begin),
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(a_end-a_begin), // relative overlap (ratio of a or c length in case of diff. read len)
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(c_end-c_begin),
            float(std::max(a_begin-pos, pos-a_end))/(a_end-a_begin),
            float(std::max(a_begin-pos, pos-a_end))/(c_end-c_begin)
        };
    }
};

struct extract_anchor_transformed_normed_weights {
    using features_t = std::array<float, 37>;

    constexpr operator auto() {
        return u8"37 extract_anchor_transformed_normed_weights";
    }

    features_t operator()(const CpuErrorCorrectorTask& task, int i, const ProgramOptions& opt) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int pos = a_begin + i;
        char orig = task.decodedAnchor[i];
        float countsACGT = msa.countsA[pos] + msa.countsC[pos] + msa.countsG[pos] + msa.countsT[pos];
        float weightsACGT = msa.weightsA[pos] + msa.weightsC[pos] + msa.weightsG[pos] + msa.weightsT[pos];
        return {
            float(orig == 'A'),
            float(orig == 'C'),
            float(orig == 'G'),
            float(orig == 'T'),
            float(msa.consensus[pos] == 'A'),
            float(msa.consensus[pos] == 'C'),
            float(msa.consensus[pos] == 'G'),
            float(msa.consensus[pos] == 'T'),
            orig == 'A'?msa.countsA[pos]/countsACGT:0,
            orig == 'C'?msa.countsC[pos]/countsACGT:0,
            orig == 'G'?msa.countsG[pos]/countsACGT:0,
            orig == 'T'?msa.countsT[pos]/countsACGT:0,
            orig == 'A'?msa.weightsA[pos]/weightsACGT:0,
            orig == 'C'?msa.weightsC[pos]/weightsACGT:0,
            orig == 'G'?msa.weightsG[pos]/weightsACGT:0,
            orig == 'T'?msa.weightsT[pos]/weightsACGT:0,
            msa.consensus[pos] == 'A'?msa.countsA[pos]/countsACGT:0,
            msa.consensus[pos] == 'C'?msa.countsC[pos]/countsACGT:0,
            msa.consensus[pos] == 'G'?msa.countsG[pos]/countsACGT:0,
            msa.consensus[pos] == 'T'?msa.countsT[pos]/countsACGT:0,
            msa.consensus[pos] == 'A'?msa.weightsA[pos]/weightsACGT:0,
            msa.consensus[pos] == 'C'?msa.weightsC[pos]/weightsACGT:0,
            msa.consensus[pos] == 'G'?msa.weightsG[pos]/weightsACGT:0,
            msa.consensus[pos] == 'T'?msa.weightsT[pos]/weightsACGT:0,
            msa.weightsA[pos]/weightsACGT,
            msa.weightsC[pos]/weightsACGT,
            msa.weightsG[pos]/weightsACGT,
            msa.weightsT[pos]/weightsACGT,
            msa.countsA[pos]/countsACGT,
            msa.countsC[pos]/countsACGT,
            msa.countsG[pos]/countsACGT,
            msa.countsT[pos]/countsACGT,
            task.msaProperties.avg_support,
            task.msaProperties.min_support,
            float(task.msaProperties.max_coverage)/opt.estimatedCoverage,
            float(task.msaProperties.min_coverage)/opt.estimatedCoverage,
            float(std::max(a_begin-pos, pos-a_end))/(a_end-a_begin)
        };
    }
};

struct extract_cands_transformed_normed_weights {
    using features_t = std::array<float, 42>;

    constexpr operator auto() {
        return u8"42 extract_cands_transformed_normed_weights";
    }

    features_t operator()(const CpuErrorCorrectorTask& task, size_t i, const ProgramOptions& opt, size_t cand, size_t offset) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int c_begin = a_begin + task.alignmentShifts[cand];
        int c_end = c_begin + task.candidateSequencesLengths[cand];
        int pos = c_begin + i;
        char orig = task.decodedCandidateSequences[offset+i];
        float countsACGT = msa.countsA[pos] + msa.countsC[pos] + msa.countsG[pos] + msa.countsT[pos];
        float weightsACGT = msa.weightsA[pos] + msa.weightsC[pos] + msa.weightsG[pos] + msa.weightsT[pos];
        MSAProperties props = msa.getMSAProperties(c_begin, c_end, opt.estimatedErrorrate, opt.estimatedCoverage, opt.m_coverage);

        return {
            float(orig == 'A'),
            float(orig == 'C'),
            float(orig == 'G'),
            float(orig == 'T'),
            float(msa.consensus[pos] == 'A'),
            float(msa.consensus[pos] == 'C'),
            float(msa.consensus[pos] == 'G'),
            float(msa.consensus[pos] == 'T'),
            orig == 'A'?msa.countsA[pos]/countsACGT:0,
            orig == 'C'?msa.countsC[pos]/countsACGT:0,
            orig == 'G'?msa.countsG[pos]/countsACGT:0,
            orig == 'T'?msa.countsT[pos]/countsACGT:0,
            orig == 'A'?msa.weightsA[pos]/weightsACGT:0,
            orig == 'C'?msa.weightsC[pos]/weightsACGT:0,
            orig == 'G'?msa.weightsG[pos]/weightsACGT:0,
            orig == 'T'?msa.weightsT[pos]/weightsACGT:0,
            msa.consensus[pos] == 'A'?msa.countsA[pos]/countsACGT:0,
            msa.consensus[pos] == 'C'?msa.countsC[pos]/countsACGT:0,
            msa.consensus[pos] == 'G'?msa.countsG[pos]/countsACGT:0,
            msa.consensus[pos] == 'T'?msa.countsT[pos]/countsACGT:0,
            msa.consensus[pos] == 'A'?msa.weightsA[pos]/weightsACGT:0,
            msa.consensus[pos] == 'C'?msa.weightsC[pos]/weightsACGT:0,
            msa.consensus[pos] == 'G'?msa.weightsG[pos]/weightsACGT:0,
            msa.consensus[pos] == 'T'?msa.weightsT[pos]/weightsACGT:0,
            msa.weightsA[pos]/weightsACGT,
            msa.weightsC[pos]/weightsACGT,
            msa.weightsG[pos]/weightsACGT,
            msa.weightsT[pos]/weightsACGT,
            msa.countsA[pos]/countsACGT,
            msa.countsC[pos]/countsACGT,
            msa.countsG[pos]/countsACGT,
            msa.countsT[pos]/countsACGT,
            props.avg_support,
            props.min_support,
            float(props.max_coverage)/opt.estimatedCoverage,
            float(props.min_coverage)/opt.estimatedCoverage,
            float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(c_end-c_begin), // absolute shift (compatible with differing read lengths)
            float(std::max(std::abs(c_begin-a_begin), std::abs(a_end-c_end)))/(a_end-a_begin),
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(a_end-a_begin), // relative overlap (ratio of a or c length in case of diff. read len)
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(c_end-c_begin),
            float(std::max(a_begin-pos, pos-a_end))/(a_end-a_begin),
            float(std::max(a_begin-pos, pos-a_end))/(c_end-c_begin)
        };
    }
};

struct extract_anchor_v2 {

    constexpr operator auto() {
        return u8"11 extract_anchor_v2";
    }

    using features_t = std::array<float, 11>;

    features_t operator()(const CpuErrorCorrectorTask& task, int i, const ProgramOptions& opt) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int pos = a_begin + i;
        float weightsACGT = msa.weightsA[pos] + msa.weightsC[pos] + msa.weightsG[pos] + msa.weightsT[pos];
        float weightCons = 0, countCons = 0;

        switch (msa.consensus[pos]) {
            case 'A':
                weightCons = msa.weightsA[pos];
                countCons = msa.countsA[pos];
                break;
            case 'C':
                weightCons = msa.weightsC[pos];
                countCons = msa.countsC[pos];
                break;
            case 'G':
                weightCons = msa.weightsG[pos];
                countCons = msa.countsG[pos];
                break;
            case 'T':
                weightCons = msa.weightsT[pos];
                countCons = msa.countsT[pos];
                break;
        }
        return {
            float(msa.origCoverages[pos]) / msa.coverage[pos],
            msa.origWeights[pos] / weightsACGT,
            msa.origWeights[pos] / msa.origCoverages[pos],
            countCons / msa.coverage[pos],
            msa.support[pos],
            weightCons / countCons,
            weightsACGT / msa.coverage[pos],
            task.msaProperties.avg_support,
            task.msaProperties.min_support,
            float(task.msaProperties.max_coverage)/opt.estimatedCoverage,
            float(task.msaProperties.min_coverage)/opt.estimatedCoverage
        };
    }
};

struct extract_cands_v2 {

    constexpr operator auto() {
        return u8"12 extract_cands_v2";
    }

    using features_t = std::array<float, 12>;

    features_t operator()(const CpuErrorCorrectorTask& task, size_t i, const ProgramOptions& opt, size_t cand, size_t offset) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int c_begin = a_begin + task.alignmentShifts[cand];
        int c_end = c_begin + task.candidateSequencesLengths[cand];
        int pos = c_begin + i;
        float weightsACGT = msa.weightsA[pos] + msa.weightsC[pos] + msa.weightsG[pos] + msa.weightsT[pos];
        MSAProperties props = msa.getMSAProperties(c_begin, c_end, opt.estimatedErrorrate, opt.estimatedCoverage, opt.m_coverage);
        float weightCons = 0, countCons = 0, weightOrig = 0, countOrig = 0;

        switch (msa.consensus[pos]) {
            case 'A':
                weightCons = msa.weightsA[pos];
                countCons = msa.countsA[pos];
                break;
            case 'C':
                weightCons = msa.weightsC[pos];
                countCons = msa.countsC[pos];
                break;
            case 'G':
                weightCons = msa.weightsG[pos];
                countCons = msa.countsG[pos];
                break;
            case 'T':
                weightCons = msa.weightsT[pos];
                countCons = msa.countsT[pos];
                break;
        }

        switch (task.decodedCandidateSequences[offset+i]) {
            case 'A':
                weightOrig = msa.weightsA[pos];
                countOrig = msa.countsA[pos];
                break;
            case 'C':
                weightOrig = msa.weightsC[pos];
                countOrig = msa.countsC[pos];
                break;
            case 'G':
                weightOrig = msa.weightsG[pos];
                countOrig = msa.countsG[pos];
                break;
            case 'T':
                weightOrig = msa.weightsT[pos];
                countOrig = msa.countsT[pos];
                break;
        }

        return {
            countOrig / msa.coverage[pos],
            weightOrig / weightsACGT,
            weightOrig / countOrig,

            countCons / msa.coverage[pos],
            msa.support[pos],
            weightCons / countCons,

            weightsACGT / msa.coverage[pos],
            props.avg_support,
            props.min_support,
            float(props.max_coverage)/opt.estimatedCoverage,
            float(props.min_coverage)/opt.estimatedCoverage,
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(std::max(a_end, c_end)-std::min(a_begin, c_begin)) // jaccard
        };
    }
};

struct extract_anchor_support {

    constexpr operator auto() {
        return u8"1 extract_anchor_support";
    }

    using features_t = std::array<float, 1>;

    features_t operator()(const CpuErrorCorrectorTask& task, int i, const ProgramOptions&) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int pos = a_begin + i;

        return {
            msa.support[pos],
        };
    }
};

struct extract_cands_support {

    constexpr operator auto() {
        return u8"1 extract_cands_support";
    }

    using features_t = std::array<float, 1>;

    features_t operator()(const CpuErrorCorrectorTask& task, size_t i, const ProgramOptions&, size_t cand, size_t) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int c_begin = a_begin + task.alignmentShifts[cand];
        int pos = c_begin + i;

        return {
            msa.support[pos],
        };
    }
};

struct extract_anchor_v3 {

    constexpr operator auto() {
        return u8"13 extract_anchor_v3";
    }

    using features_t = std::array<float, 13>;

    features_t operator()(const CpuErrorCorrectorTask& task, int i, const ProgramOptions& opt) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int pos = a_begin + i;
        float weightsACGT = msa.weightsA[pos] + msa.weightsC[pos] + msa.weightsG[pos] + msa.weightsT[pos];
        float weightCons = 0, countCons = 0;

        switch (msa.consensus[pos]) {
            case 'A':
                weightCons = msa.weightsA[pos];
                countCons = msa.countsA[pos];
                break;
            case 'C':
                weightCons = msa.weightsC[pos];
                countCons = msa.countsC[pos];
                break;
            case 'G':
                weightCons = msa.weightsG[pos];
                countCons = msa.countsG[pos];
                break;
            case 'T':
                weightCons = msa.weightsT[pos];
                countCons = msa.countsT[pos];
                break;
        }
        return {
            float(msa.origCoverages[pos]) / msa.coverage[pos],
            msa.origWeights[pos] / weightsACGT,
            msa.origWeights[pos] / msa.origCoverages[pos],
            countCons / msa.coverage[pos],
            msa.support[pos],
            weightCons / countCons,
            weightsACGT / msa.coverage[pos],
            msa.coverage[pos] / opt.estimatedCoverage,
            weightsACGT / opt.estimatedCoverage,
            task.msaProperties.avg_support,
            task.msaProperties.min_support,
            float(task.msaProperties.max_coverage)/opt.estimatedCoverage,
            float(task.msaProperties.min_coverage)/opt.estimatedCoverage
        };
    }
};

struct extract_cands_v3 {

    constexpr operator auto() {
        return u8"14 extract_cands_v3";
    }

    using features_t = std::array<float, 14>;

    features_t operator()(const CpuErrorCorrectorTask& task, size_t i, const ProgramOptions& opt, size_t cand, size_t offset) noexcept {   
        auto& msa = task.multipleSequenceAlignment;
        int a_begin = msa.anchorColumnsBegin_incl;
        int a_end = msa.anchorColumnsEnd_excl;
        int c_begin = a_begin + task.alignmentShifts[cand];
        int c_end = c_begin + task.candidateSequencesLengths[cand];
        int pos = c_begin + i;
        float weightsACGT = msa.weightsA[pos] + msa.weightsC[pos] + msa.weightsG[pos] + msa.weightsT[pos];
        MSAProperties props = msa.getMSAProperties(c_begin, c_end, opt.estimatedErrorrate, opt.estimatedCoverage, opt.m_coverage);
        float weightCons = 0, countCons = 0, weightOrig = 0, countOrig = 0;

        switch (msa.consensus[pos]) {
            case 'A':
                weightCons = msa.weightsA[pos];
                countCons = msa.countsA[pos];
                break;
            case 'C':
                weightCons = msa.weightsC[pos];
                countCons = msa.countsC[pos];
                break;
            case 'G':
                weightCons = msa.weightsG[pos];
                countCons = msa.countsG[pos];
                break;
            case 'T':
                weightCons = msa.weightsT[pos];
                countCons = msa.countsT[pos];
                break;
        }

        switch (task.decodedCandidateSequences[offset+i]) {
            case 'A':
                weightOrig = msa.weightsA[pos];
                countOrig = msa.countsA[pos];
                break;
            case 'C':
                weightOrig = msa.weightsC[pos];
                countOrig = msa.countsC[pos];
                break;
            case 'G':
                weightOrig = msa.weightsG[pos];
                countOrig = msa.countsG[pos];
                break;
            case 'T':
                weightOrig = msa.weightsT[pos];
                countOrig = msa.countsT[pos];
                break;
        }

        return {
            countOrig / msa.coverage[pos],
            weightOrig / weightsACGT,
            weightOrig / countOrig,

            countCons / msa.coverage[pos],
            msa.support[pos],
            weightCons / countCons,

            weightsACGT / msa.coverage[pos],

            msa.coverage[pos] / opt.estimatedCoverage,
            weightsACGT / opt.estimatedCoverage,

            props.avg_support,
            props.min_support,
            float(props.max_coverage)/opt.estimatedCoverage,
            float(props.min_coverage)/opt.estimatedCoverage,
            float(std::min(a_end, c_end)-std::max(a_begin, c_begin))/(std::max(a_end, c_end)-std::min(a_begin, c_begin)) // jaccard
        };
    }
};



} //namespace detail


//--------------------------------------------------------------------------------

using anchor_extractor = detail::extract_anchor_v3;
using cands_extractor = detail::extract_cands_v3;

using anchor_clf_t = ForestClf<anchor_extractor>;
using cands_clf_t = ForestClf<cands_extractor>;

using ClfAgent = clf_agent<anchor_clf_t, cands_clf_t, anchor_extractor, cands_extractor>;

} // namespace care

#endif
