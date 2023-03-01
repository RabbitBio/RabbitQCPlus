#ifndef __MAIN_CORRECT_CPU_H__
#define __MAIN_CORRECT_CPU_H__


int main_correction(int argc, char** argv, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q, std::atomic_int *producerDone);

#endif
