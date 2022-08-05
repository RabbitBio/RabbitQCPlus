//
// Created by ylf9811 on 2021/11/3.
//

#ifndef RERABBITQC_PUGZ_H
#define RERABBITQC_PUGZ_H


#include "Globals.h"


void main_pugz(std::string in_name, int threads, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q,
               std::atomic_int *producerDone);

#endif//RERABBITQC_PUGZ_H
