/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "Sort.h"
#include <sstream>
#include <iostream>
#include <algorithm>
/*
bool PairsFileSorter::CompareFunction (const PairsRecord *a, const PairsRecord *b){
    static int g_Order [] = {0,2,1,3};
    for (int i=0; i<4; i++){
        int r = Spaceship(a->pos[g_Order[i]],b->pos[g_Order[i]]);
        if (r!=0) return 1-r;
    }
    return false;
}
*/

bool PairsFileSorter::CompareFunction (const PairsRecord *a, const PairsRecord *b){
    for (auto it = this->order.begin(); it != this->order.end(); ++it){
        if ((*it) > 0){
            int r = Spaceship(a->pos[(*it)-1],b->pos[(*it)-1]);
            if (r!=0) return 1-r;
        }else{
            if ((*it) == 0) return a->readID < b->readID;
            else{
                //TODO
            }
        }
    }
    return false;
}


void PairsFileSorter::Sort(const int compareOrder [], const unsigned int n, const int threads){
    for (unsigned int i=0; i<n; ++i){
        order.push_back(compareOrder[i]);
    }
    
    std::cerr << "Sort start\n";
    //std::sort(records.begin(),records.end(), [this](const PairsRecord *a, const PairsRecord *b) {return CompareFunction(a, b); });
    using namespace std::placeholders;
    std::sort(records.begin(),records.end(), std::bind(&PairsFileSorter::CompareFunction, this, _1, _2));
    std::cerr << "Sort end\n";
    
}

void PairsFileSorter::AddHeader(PairsFileHeader &header){
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

void PairsFileSorter::PrintRecords(){
    for (auto it = records.begin(); it != records.end(); ++it){
        std::cout << (*it)->readID;
        for (unsigned int i=0; i<2; ++i){
            std::cout << '\t' << chromosomes[(*it)->pos[i*2]].chrom << '\t' << (*it)->pos[i*2+1];
        }
        for (unsigned int i=0; i<nfields-5; ++i){
            std::cout << '\t' << (*it)->extrafields[i];
        }
        std::cout << '\n';
    }
}