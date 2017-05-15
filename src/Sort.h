/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#ifndef _SORT
#define _SORT
#include "PairsFileRecord.h"
#include "PairsFileHeader.h"
#include <vector>
#include <string>

inline int Spaceship(const int &lhs, const int &rhs) {
    if (lhs < rhs) return -1;
    if (rhs < lhs) return 1;
    return 0;
}

class PairsFileSorter{
public:
    void AddHeader(PairsFileHeader &header);
    bool AddRecord(const std::string &line);
    void Sort(const int compareOrder [], const unsigned int n, const int threads);
    void PrintRecords();
    
private:
    bool CompareFunction (const PairsRecord *a, const PairsRecord *b);
    
    std::vector<PairsRecord*> records;
    std::vector<Chromosome> chromosomes;
    std::map<std::string, int> chromMap;
    std::vector<int> order;
    unsigned int nfields;
    int shape;
};

#endif
