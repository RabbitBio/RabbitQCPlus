//
// Created by ylf9811 on 2021/10/14.
//

#ifndef RERABBITQC_PIGZ_H
#define RERABBITQC_PIGZ_H

#include "Globals.h"


int main_pigz(int argc, char **argv, moodycamel::ReaderWriterQueue<pair<char *, int>> *Q,
              std::atomic_int *wDone, pair<char *, int> &L);

#endif //RERABBITQC_PIGZ_H
