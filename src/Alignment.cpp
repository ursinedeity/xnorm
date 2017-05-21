/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */
#include "Alignment.h"
#include <sstream>

bool CompareAllowOneMismatch(const std::string &s, unsigned int l, const std::string &t, int linc, int rinc){
    //compare from left side
    int i,j;
    for (i = 0; (i < t.length()) && (s[i+l] == t[i]); ++i);
    if (i == t.length()) 
        return true;
    for (j = t.length()-1; (j >= 0) && (s[j+l] == t[j]); ++j);
    if ((i == j) || (j-i == 1))
        return true;
    
    for (i = 0; (i < t.length()-1) && (s[i+l+linc] == t[i]); ++i);
    for (j = t.length()-1; (j > 0) && (s[j+l+rinc] == t[j]); ++j);

    if (i >= j)
        return true;
    
    return false;
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

bool Alignment::HasCuttingSiteSignature(const std::string &enzymeSite, const std::string &enzymeSiteReverse){
    unsigned int oneSideLength = enzymeSiteReverse.length();
    if (strand){ //reverse strand
        return CompareAllowOneMismatch(sequence, sequence.length()-oneSideLength, enzymeSiteReverse, 1, 0);
        //return sequence.compare(sequence.length()-oneSideLength, oneSideLength, enzymeSiteReverse ) == 0;
    }else{
        return CompareAllowOneMismatch(sequence, 0, enzymeSite, 0, -1);
        //return sequence.compare(0, oneSideLength, enzymeSite) == 0;
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
