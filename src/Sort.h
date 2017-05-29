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
#include <memory>
inline int Spaceship(const int &lhs, const int &rhs) {
    if (lhs < rhs) return -1;
    if (rhs < lhs) return 1;
    return 0;
}

class PairsFileSorter{
public:
    PairsFileSorter(PairsFileHeader &header);
    bool AddRecord(const std::string &line);
    void Sort(const int compareOrder [], const unsigned int n, const unsigned int threads);
    void PrintRecords(std::shared_ptr<std::ostream> &output);
    ~PairsFileSorter(){
        for (auto it = records.begin(); it != records.end(); ++it)
            delete (*it);
    }
private:
    bool CompareFunction (const PairsRecord *a, const PairsRecord *b);
    void ThreadMergeSort(std::vector<unsigned int> &range, unsigned int l, unsigned int r);
    
    std::string headerRepresentation;
    std::vector<PairsRecord*> records;
    std::vector<Chromosome> chromosomes;
    std::map<std::string, int> chromMap;
    std::vector<int> order;
    unsigned int nfields;
    int shape;
};

#endif
