/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the license.txt file included.
 */

#ifndef _XNORM_T
#define _XNORM_T

#include <iostream>
#include <fstream>

#include <string>
#include <vector>

#define UNMAPPED 4
#define REVERSE 16
#define MAXINSERT 2147183647

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


/**********************************************************/

struct Chromosome{
    std::string chrom;
    int length;
    Chromosome(const std::string &chrom, const int &length){
        this->chrom = chrom;
        this->length = length;
    }
    bool operator < (const Chromosome &a) const{
        return ChromosomeOrder(chrom,a.chrom);
    }
};

/**********************************************************/
//Sorting pair of alignments and process them
class PairReads{
public:
    bool AddAlignment(const std::string &line);
    void AddHeader(const std::string &line);
    
    PairReads(const std::string &pairFile,
              const std::string &controlFile,
              const std::string &enzymeSite,
              const std::string &ligationJunction,
              const int &insertLength,
              const int &mapqCutoff);
    ~PairReads();
    void PrintStats();
private:
    std::ofstream pairlist,control;
    std::string enzymeSite,ligationJunction;
    int insertLength,mapqCutoff;
    Alignment *last,*current;
    bool isWaiting;
    
    std::vector<Chromosome> chromosomes;
    bool headerChrInProcess=true;
    
    int totalCount=0,hicCount=0,ctlCount=0,rlgCount=0,sglCount=0,jkCount=0;
    void ProcessPair();
    void WriteControl(const Alignment &alignment);
    void WritePairlist();
    void OutputHeader();
};

/*********************************************************/
#endif //_XNORM_T
