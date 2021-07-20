//
// Created by ylf9811 on 2021/7/14.
//

#ifndef RABBITQCPLUS_POLYX_H
#define RABBITQCPLUS_POLYX_H


#include "Reference.h"

class PolyX {
public:
    PolyX();

    ~PolyX();

    static void trimPolyG(neoReference &r1, neoReference &r2, int compareReq);

    static void trimPolyG(neoReference &ref, int compareReq);

    static void trimPolyX(neoReference &r1, neoReference &r2, int compareReq);

    static void trimPolyX(neoReference &ref, int compareReq);


};


#endif //RABBITQCPLUS_POLYX_H
