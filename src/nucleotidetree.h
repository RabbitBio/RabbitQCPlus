//
// Created by yanli on 2021/7/20.
//

#ifndef RERABBITQC_NUCLEOTIDETREE_H
#define RERABBITQC_NUCLEOTIDETREE_H

#include <cstring>
#include <sstream>
#include "Globals.h"

class NucleotideNode {
public:
    NucleotideNode();

    ~NucleotideNode();

    void dfs();

public:
    int count;
    char base;
    NucleotideNode *children[8];
};

class NucleotideTree {
public:
    NucleotideTree();

    ~NucleotideTree();

    void addSeq(std::string seq);

    std::string getDominantPath(bool &reachedLeaf);

    static bool test();

private:
    NucleotideNode *mRoot;
};


#endif //RERABBITQC_NUCLEOTIDETREE_H
