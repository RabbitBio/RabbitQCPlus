//
// Created by ylf9811 on 2021/11/3.
//

#ifndef RERABBITQC_PUGZ_H
#define RERABBITQC_PUGZ_H


#include "Globals.h"



void main_pugz(string in_name, int threads, moodycamel::ReaderWriterQueue<pair<char *, int>> *Q,
               atomic_int *producerDone);

#endif //RERABBITQC_PUGZ_H
