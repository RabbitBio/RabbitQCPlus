
#ifndef RERABBITQC_PRAGZIP_H
#define RERABBITQC_PRAGZIP_H


#include "Globals.h"


int main_pragzip(int argc, char *argv[], moodycamel::ReaderWriterQueue<std::pair<char *, int>>*Q, std::atomic_int *producerDone);

#endif//RERABBITQC_PRAGZIP_H
