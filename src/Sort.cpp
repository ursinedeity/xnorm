/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "Sort.h"
#include <sstream>

PairsFileSorter::PairsFileSorter(PairsFileHeader &header){
    header.sort_chromosome();
    this->chromosomes = header.chromosomes;
    this->chromMap = header.chromMap;
    this->nfields = header.columns.size();
    if (header.shape.substr(8) == "upper triangle")
        this->shape = 1;
    else{
        if (header.shape.substr(8) == "lower triangle")
            this->shape = -1;
        else
            this->shape = 0;
    }
}
bool PairsFileSorter::AddRecord(const std::string &line){
    std::istringstream linestream(line);
    std::string readID;
    std::string *extrafields;
    unsigned int combinedPos[4];
    
    linestream >> readID;
    for (int i=0; i<2; ++i){
        std::string chrom;
        unsigned int pos;
        
        linestream >> chrom >> pos;

        auto it = chromMap.find(chrom);
        if (it == chromMap.end()) return false;
        
        combinedPos[i*2] = chromMap[chrom];
        combinedPos[i*2+1] = pos;
    }
    if (nfields > 5){
        extrafields = new std::string[nfields-5];
        for (unsigned int i=0; i<nfields-5; ++i) linestream >> extrafields[i];
    }else{
        extrafields = NULL;
    }
    
    records.push_back(new PairsRecord(readID,combinedPos,extrafields,shape));
    return true;
}
