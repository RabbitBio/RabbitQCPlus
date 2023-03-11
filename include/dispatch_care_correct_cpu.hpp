#ifndef CARE_DISPATCH_CARE_CORRECT_CPU_HPP
#define CARE_DISPATCH_CARE_CORRECT_CPU_HPP

#include <options.hpp>
#include <concurrencyhelpers.hpp>

namespace care{

    void performCorrection(ProgramOptions programOptions);

    void performCorrectionOutToQueue(ProgramOptions programOptions, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q1, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q2, std::atomic_int *producerDone, std::atomic_int *careStartWrite, int* changNum);

}



#endif
