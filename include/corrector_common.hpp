#ifndef CARE_CORRECTOR_COMMON_HPP
#define CARE_CORRECTOR_COMMON_HPP

#include <config.hpp>
#include <correctedsequence.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <vector>

namespace care{

    class CorrectionOutput{
    public:
        std::vector<TempCorrectedSequence> anchorCorrections;
        std::vector<TempCorrectedSequence> candidateCorrections;
    };

    class EncodedCorrectionOutput{
    public:
        EncodedCorrectionOutput() = default;
        EncodedCorrectionOutput(const EncodedCorrectionOutput&) = default;
        EncodedCorrectionOutput(EncodedCorrectionOutput&&) = default;

        EncodedCorrectionOutput(const CorrectionOutput& rhs){
            encodedAnchorCorrections.resize(rhs.anchorCorrections.size());
            encodedCandidateCorrections.resize(rhs.candidateCorrections.size());

            for(std::size_t i = 0; i < rhs.anchorCorrections.size(); i++){
                rhs.anchorCorrections[i].encodeInto(encodedAnchorCorrections[i]);
            }

            for(std::size_t i = 0; i < rhs.candidateCorrections.size(); i++){
                rhs.candidateCorrections[i].encodeInto(encodedCandidateCorrections[i]);
            }
        }

        EncodedCorrectionOutput& operator=(EncodedCorrectionOutput rhs){
            std::swap(*this, rhs);
            return *this;
        }

        EncodedCorrectionOutput& operator=(const CorrectionOutput& rhs){
            EncodedCorrectionOutput tmp(rhs);
            std::swap(*this, tmp);
            return *this;
        }

        std::vector<EncodedTempCorrectedSequence> encodedAnchorCorrections;
        std::vector<EncodedTempCorrectedSequence> encodedCandidateCorrections;
    };

    class ReadCorrectionFlags{
    public:
        ReadCorrectionFlags() = default;

        ReadCorrectionFlags(std::size_t numReads)
            : size(numReads), flags(std::make_unique<std::uint8_t[]>(numReads)){
            std::fill(flags.get(), flags.get() + size, 0);
        }

        std::size_t sizeInBytes() const noexcept{
            return size * sizeof(std::uint8_t);
        }

        bool isCorrectedAsHQAnchor(std::int64_t position) const noexcept{
            return (flags[position] & readCorrectedAsHQAnchor()) > 0;
        }

        bool isNotCorrectedAsAnchor(std::int64_t position) const noexcept{
            return (flags[position] & readCouldNotBeCorrectedAsAnchor()) > 0;
        }

        void setCorrectedAsHqAnchor(std::int64_t position) const noexcept{
            flags[position] = readCorrectedAsHQAnchor();
        }

        void setCouldNotBeCorrectedAsAnchor(std::int64_t position) const noexcept{
            flags[position] = readCouldNotBeCorrectedAsAnchor();
        }

    private:
        static constexpr std::uint8_t readCorrectedAsHQAnchor() noexcept{ return 1; };
        static constexpr std::uint8_t readCouldNotBeCorrectedAsAnchor() noexcept{ return 2; };

        std::size_t size;
        std::unique_ptr<std::uint8_t[]> flags{};
    };


}


#endif