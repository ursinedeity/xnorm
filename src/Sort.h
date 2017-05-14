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

class PairsFileSorter{
public:
    PairsFileSorter(PairsFileHeader &header);
    bool AddRecord(const std::string &line);
    void Sort(const int compareOrder [], const int threads);
private:
    std::vector<PairsRecord*> records;
    std::vector<Chromosome> chromosomes;
    std::map<std::string, int> chromMap;
    unsigned int nfields;
    int shape;
};

#endif
