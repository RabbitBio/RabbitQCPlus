#ifndef CARE_CORRECTION_RESULT_OUTPUT_HPP
#define CARE_CORRECTION_RESULT_OUTPUT_HPP

#include <config.hpp>

#include <readlibraryio.hpp>
#include <serializedobjectstorage.hpp>
#include <options.hpp>
#include <concurrencyhelpers.hpp>

#include <string>
#include <vector>

namespace care{

    void constructOutputFileFromCorrectionResults(
        const std::vector<std::string>& originalReadFiles,
        SerializedObjectStorage& partialResults, 
        FileFormat outputFormat,
        const std::vector<std::string>& outputfiles,
        bool showProgress,
        const ProgramOptions& programOptions
    );



    void constructOutputFileFromCorrectionResultsOutToQueue(
        const std::vector<std::string>& originalReadFiles,
        SerializedObjectStorage& partialResults, 
        FileFormat outputFormat,
        const std::vector<std::string>& outputfiles,
        bool showProgress,
        const ProgramOptions& programOptions,
        moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q,
        std::atomic_int *producerDone
    );



}


#endif
