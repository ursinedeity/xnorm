/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "PairsFileRecord.h"
#include <sstream>

PairsRecord::PairsRecord(const std::string &readID, const unsigned int pos [], std::string *extrafields, const int shape){
    this->readID = readID;
    for (int i=0; i<4; i++) this->pos[i] = pos[i];
    
    if (this->pos[0] == this->pos[2]){
        if ((this->pos[1] - this->pos[3])* shape > 0) SwapPos();
    }else{
        if ((this->pos[0] - this->pos[2])* shape > 0) SwapPos();
    }
    
    this->extrafields = extrafields;
}

void PairsRecord::SwapPos(){
    unsigned int t;
    t=pos[0]; pos[0]=pos[2]; pos[2]=t;
    t=pos[1]; pos[1]=pos[3]; pos[3]=t;
}
