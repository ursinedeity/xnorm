/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the license file included.
 */
#include "PairReads.h"
#include <iostream>
#include <sstream>
#include <algorithm>

int InsertSize(const Alignment &p, const Alignment &q){
    if (p.rname != q.rname)
        return MAXINSERT;
    else
        return abs(q.pos - p.pos);
}

int CigarFirstPos(const std::string &cigar){
    if (cigar[0] == '*')
        return 0;
    else
        return std::stoi(cigar);    
}

bool FindLocalLigationJunction(const std::string &s, const int &p, const std::string &j){
//     std::string local = s.substr(std::max(0,          p-int(j.length())),
//                                         std::min(int(j.length()), int(s.length())-p));
    std::size_t found = s.find(j);
    if (found != std::string::npos) return true;
    else return false;
}

/**********************************************/
//Class PairReads

void PairReads::ProcessPair(){
    
    totalCount++;
    if (totalCount % 10000000 == 0){
        std::cerr << "Read pairs processed: " << totalCount/1000000 << 'M' << std::endl;
    }
    if (last->isValid | current->isValid){
        if ((InsertSize(*last,*current) > insertLength) && (last->isValid & current->isValid)){
            //valid HiC pairs
            hicCount++;
            WritePairlist();
            
        }else{
            if (last->HasCuttingSiteSignature(enzymeSite) || current->HasCuttingSiteSignature(enzymeSite)){
                //potential control reads
                
                //Check if there is a religation site in any of the sequences
                unsigned int plast = CigarFirstPos(last->cigar);
                unsigned int pcurrent = CigarFirstPos(current->cigar);
                //Check last
                if ((plast != 0) && (plast != last->sequence.length())){
                    if (FindLocalLigationJunction(last->sequence, plast, ligationJunction)){
                        rlgCount++;
                        //std::cout << last->Representation() << std::endl;
                        return;
                    }
                }
                //Check current
                if ((pcurrent != 0) && (pcurrent != current->sequence.length())){
                    if (FindLocalLigationJunction(current->sequence, pcurrent, ligationJunction)){
                        rlgCount++;
                        //std::cout << current->Representation() << std::endl;
                        return;
                    }
                }
                
                //Check whether they are control pairs or reads
                if (last->isValid & current->isValid){
                    if (last->pos*(last->strand*2-1) + current->pos*(current->strand*2-1) <= 0){
                        //check if it's FR orientaion
                        ctlCount = ctlCount+2;
                        WriteControl(*last);
                        WriteControl(*current);
                    }
                }else{
                    if (last->isValid){ctlCount++; WriteControl(*last);} 
                    if (current->isValid){ctlCount++; WriteControl(*current);}
                }
            }else{
                if (last->isValid & current->isValid){
                    //religation
                    rlgCount++;
                    
                }else{
                    //single end reads
                    sglCount++;
                    
                }
            }
        }
    }else{
        //junk reads
        jkCount++;
        
    }
}


PairReads::PairReads(const std::string &pairFile,
                     const std::string &controlFile,
                     const std::string &enzymeSite,
                     const std::string &ligationJunction,
                     const int &insertLength,
                     const int &mapqCutoff)
{
    pairlist.open(pairFile.c_str());
    if (!pairlist.is_open()){
        std::cerr << "Unable to write file " << pairFile << std::endl;
        exit(1);
    }
    control.open(controlFile.c_str());
    if (!control.is_open()){
        std::cerr << "Unable to write file " << controlFile << std::endl;
        exit(1);
    }
    this->enzymeSite = enzymeSite;
    this->ligationJunction = ligationJunction;
    this->insertLength = insertLength;
    this->mapqCutoff = mapqCutoff;
    //initialize variables
    isWaiting = false;
    headerChrInProcess=true;
    totalCount=0,hicCount=0,ctlCount=0,rlgCount=0,sglCount=0,jkCount=0;
    last = new Alignment; current = new Alignment;
}

PairReads::~PairReads(){
    pairlist.close();
    control.close();
    delete last;
    delete current;
}

bool Mates(const Alignment &p, const Alignment &q){
    if (p.qname == q.qname)
        return true;
    else
        return false;
}

bool PairReads::AddAlignment(const std::string &line){
    //write header to pairs file
    if (headerChrInProcess){OutputHeader();}
    
    current->ParseAlignment(line,mapqCutoff);
    
    if (isWaiting){
        if (Mates(*last,*current)){
            ProcessPair();
            isWaiting = false;
            return false;
        }else{
            std::cerr << "Mate mismatch found. Alignments are not fully paired.\n";
        }
    }
    static Alignment *tmp;
    tmp = last; last = current; current = tmp;
    //SwapAlignments(last,current);
    isWaiting = true;
    return true;
}

void PairReads::AddHeader(const std::string &line){
    std::istringstream linestream(line);
    std::string tag,field[2],chrom="";
    int length=0;
    
    linestream >> tag >> field[0] >> field[1];
    if (tag == "@SQ"){
        headerChrInProcess = true;
        for (int i=0; i<2; ++i){
            if (field[i].substr(0,2) == "SN") chrom = field[i].substr(3);
            if (field[i].substr(0,2) == "LN") length = std::stoi(field[i].substr(3));
        }
        chromosomes.push_back(Chromosome(chrom, length));
    }
}

void PairReads::OutputHeader(){
    std::sort (chromosomes.begin(), chromosomes.end());
    pairlist << "## pairs format v1.0\n";
    pairlist << "#sorted: none\n";
    pairlist << "#shape: upper triangle\n";
    
    if (chromosomes.size() == 0) std::cerr << "chromosome sizes not found, header invalid.\n";
    
    for (std::vector<Chromosome>::iterator it=chromosomes.begin(); it!=chromosomes.end(); ++it){
        pairlist << "#chromsize: " << it->chrom << " " << it->length << '\n';
    }
    pairlist << "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n";
    headerChrInProcess = false;
}

void PairReads::WritePairlist(){
    
    if (*last > *current){
        static Alignment *tmp;
        tmp = last; last = current; current = tmp;
    }
    
    pairlist << last->qname << '\t' 
             << last->rname << '\t' << last->pos << '\t'
             << current->rname << '\t' << current->pos << '\t'
             << "+-"[last->strand] << '\t' << "+-"[current->strand] << '\n';
    
}
void PairReads::WriteControl(const Alignment &p){
    control << p.qname << '\t' << p.rname << '\t' << p.pos << '\t' << "+-"[p.strand] << '\n';
}

void PairReads::PrintStats(){
    std::cerr << "total number of read pairs processed: " << totalCount << std::endl;
    std::cerr << "number of Hi-C contacts: " << hicCount << std::endl;
    std::cerr << "number of single end read pairs: " << sglCount << std::endl;
    std::cerr << "number of control reads: " << ctlCount << std::endl;
    std::cerr << "number of re-ligation read pairs: " << rlgCount << std::endl;
    std::cerr << "number of junk read pairs: " << jkCount << std::endl;
    
}
