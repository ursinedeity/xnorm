/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the license.txt file included.
 */
#include "xnorm.h"
#include <sstream>

bool ChromosomeOrder(const std::string &chrom1, const std::string &chrom2){
    if (chrom1.length() != chrom2.length()){
        if (chrom1.length() < chrom2.length()){ //chr? chr??
            if (chrom1[3] > 'A') return false;//chrX chr11
            else return true;                 //chr2 chr11
        }else{                                  //chr?? chr?
            if (chrom2[3] > 'A') return true; //chr11 chrX  
            else return false;                //chr11 chr2
        }
    }else
        return chrom1 <= chrom2;
}

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

bool Alignment::HasCuttingSiteSignature(const std::string &enzymeSite){
    if (sequence.compare(0,enzymeSite.length(),enzymeSite) == 0)
        return true;
    else
        return false;
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
