/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "Builder.h"
#include <sstream>
#include <iostream>
bool CooBuilder::AddPairsLine(const std::string &line){
    std::string chrom, readID;
    unsigned int pos;
    unsigned int bin[2];
    std::istringstream linestream(line);
    
    linestream >> readID;
    for (int i=0; i<2; ++i){
        linestream >> chrom >> pos;
        
        int chromnum = genome.get_chromnum(chrom);
        if (chromnum < 0)
            return false;
        
        bin[i] = pos/resolution + index.offset[chromnum];
    }
    
    if ((bin[0] == ci) && (bin[1] == cj))
        cv += 1;
    else{
        if (cv > 0){
            Ai.push_back(ci);
            Aj.push_back(cj);
            Ax.push_back(cv);
        }
        ci = bin[0]; cj = bin[1];
        cv = 0;
    }
    
    return true;
    
}


alab::Matrix * CooBuilder::BuildCsrMatrix(){
    alab::Matrix *m = new alab::Matrix;
    std::cout << index.size << ' ' << Ax.size() << '\n';
    m->LoadCoo(Ai.data(), Aj.data(), Ax.data(), index.size, Ax.size());
    return m;
}

void CooBuilder::PrintCoo(){
    for (int i = 0; i < Ax.size(); ++i){
        std::cout << Ai[i] << ' ' << Aj[i] << ' ' << Ax[i] << '\n';
    }
}
