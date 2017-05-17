/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#ifndef _ALIGNMENT
#define _ALIGNMENT
#include "Chromosome.h"
#define UNMAPPED 4
#define REVERSE 16

struct Alignment{
    bool unmapped,strand,isValid;
    int mapq,pos;
    std::string qname,rname,cigar;
    std::string sequence;
    
    void ParseAlignment(const std::string &line, const int &mapqCutoff);
    
    bool HasCuttingSiteSignature(const std::string &enzymeSite, const std::string &enzymeSiteReverse);
    bool operator > (const Alignment &a) const;
    //alignment representation
    std::string Representation();
};

#endif
