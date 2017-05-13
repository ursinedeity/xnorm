/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#ifndef _PAIRREADS
#define _PAIRREADS
#include "Alignment.h"
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
              const int &insertLength,
              const int &mapqCutoff);
    ~PairReads();
    void PrintStats();
    void PrintArgs();
private:
    std::ofstream pairlist,control;
    std::string enzymeSite,ligationJunction;
    int insertLength,mapqCutoff;
    Alignment *last,*current;
    bool isWaiting;
    
    std::vector<Chromosome> chromosomes;
    bool headerChrInProcess;
    
    int totalCount,hicCount,ctlCount,rlgCount,sglCount,jkCount;
    void ProcessPair();
    void WriteControl(const Alignment &alignment);
    void WritePairlist();
    void OutputHeader();
};

#endif
