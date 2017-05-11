/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#ifndef _ALIGNMENT
#define _ALIGNMENT

#include <string>
#define UNMAPPED 4
#define REVERSE 16

bool ChromosomeOrder(const std::string &chrom1, const std::string &chrom2);

struct Alignment{
    bool unmapped,strand,isValid;
    int mapq,pos;
    std::string qname,rname,cigar;
    std::string sequence;
    
    void ParseAlignment(const std::string &line, const int &mapqCutoff);
    
    bool HasCuttingSiteSignature(const std::string &enzymeSite);
    bool operator > (const Alignment &a) const;
    //alignment representation
    std::string Representation();
};

#endif
