/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#ifndef _CHROMOSOME
#define _CHROMOSOME
#include <string>

inline bool ChromosomeOrder(const std::string &chrom1, const std::string &chrom2){
//     if (chrom1.length() != chrom2.length()){
//         if (chrom1.length() < chrom2.length()){ //chr? chr??
//             if (chrom1[3] > 'A') return false;//chrX chr11
//             else return true;                 //chr2 chr11
//         }else{                                  //chr?? chr?
//             if (chrom2[3] > 'A') return true; //chr11 chrX  
//             else return false;                //chr11 chr2
//         }
//     }else
    return chrom1 < chrom2;
}

struct Chromosome{
    std::string chrom;
    int length,order;
    Chromosome(const std::string &chrom, const unsigned int &length){
        this->chrom = chrom;
        this->length = length;
        this->order = -1;
    }
    Chromosome(const std::string &chrom, const unsigned int &length, const int &order){
        this->chrom = chrom;
        this->length = length;
        this->order = order;
    }
    bool operator < (const Chromosome &a) const{
        if ((order >= 0) && (a.order >= 0)) return order < a.order;
        else return ChromosomeOrder(chrom,a.chrom);
    }
};

#endif
