/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#ifndef _PAIRREADS
#define _PAIRREADS
#include "Alignment.h"
#include "PairsFileHeader.h"
#include <fstream>
#include <vector>

#define MAXINSERT 2147183647

class PairReads{
public:
    bool AddAlignment(const std::string &line);
    void AddHeader(const std::string &line);
    
    PairReads(const std::string &pairFile,
              const std::string &controlFile,
              const std::string &enzyme,
              const unsigned int &insertLength,
              const unsigned int &mapqCutoff);
    ~PairReads();
    void set_genome_assembly(const std::string &genome);
    void PrintStats();
    void PrintArgs();
private:
    std::ofstream pairlist,control;
    std::string enzymeSite,ligationJunction;
    unsigned int insertLength,mapqCutoff;
    Alignment *last,*current;
    bool isWaiting;
    
    PairsFileHeader header;
    bool headerChrInProcess;
    
    unsigned int totalCount,hicCount,ctlCount,rlgCount,sglCount,jkCount;
    void ProcessPair();
    void WriteControl(const Alignment &alignment);
    void WritePairlist();
    void OutputHeader();
};

#endif
