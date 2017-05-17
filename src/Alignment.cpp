/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */
#include "Alignment.h"
#include <sstream>

//Takes a line from sam file and parse it
void Alignment::ParseAlignment(const std::string &line, const int &mapqCutoff){
    std::istringstream linestream(line);
    std::string tmp;
    
    linestream >> qname;
    int flag;
    linestream >> flag;
    unmapped = flag & UNMAPPED; //T for unmapped
    strand = flag & REVERSE; //T for reverse
    
    linestream >> rname;
    linestream >> pos;
    linestream >> mapq;
    linestream >> cigar;
    linestream >> tmp >> tmp >> tmp; //RNEXT PNEXT TLEN
    linestream >> sequence;
    if (unmapped | (mapq < mapqCutoff))
        isValid = false;
    else
        isValid = true;
}

bool Alignment::HasCuttingSiteSignature(const std::string &ligationJunction){
    unsigned int oneSideLength = ligationJunction.length() / 2;
    if (strand){ //reverse strand
        return sequence.compare(sequence.length()-oneSideLength, oneSideLength, ligationJunction.substr(0, oneSideLength)) == 0;
    }else{
        return sequence.compare(0, oneSideLength, ligationJunction.substr(oneSideLength, oneSideLength)) == 0;
    }
        
}

bool Alignment::operator > (const Alignment &a) const{
    if (rname == a.rname)
        return pos > a.pos;
    else
        return not ChromosomeOrder(rname,a.rname);
}

std::string Alignment::Representation(){
    std::string repr;
    repr += qname + '\t';
    repr += rname + '\t';
    repr += std::to_string(pos) + '\t';
    repr += cigar + '\t' + sequence ;
    return repr;
}
