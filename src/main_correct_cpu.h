#ifndef __MAIN_CORRECT_CPU_H__
#define __MAIN_CORRECT_CPU_H__


int main_correction(int argc, char** argv, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q1, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q2, std::atomic_int *producerDone, std::atomic_int *careStartWrite, int* changNum);

#endif
