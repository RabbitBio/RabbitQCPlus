//
// Created by ylf9811 on 2021/7/12.
//

#ifndef RABBITQCPLUS_REPOTER_H
#define RABBITQCPLUS_REPOTER_H

#include "Reference.h"

class Repoter {
public:
    Repoter();


    static void PrintRef(neoReference &ref);

    static void PrintRef(Reference &ref);

    static Reference GetRevRef(neoReference &ref);
};


#endif //RABBITQCPLUS_REPOTER_H
