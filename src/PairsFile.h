/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#ifndef _PAIRSFILE
#define _PAIRSFILE
#include "Chromosome.h"
#include <vector>
#include <iostream>
#include <map>

class PairsFileHeader{
public:
    int ParseHeader(const std::string &line);
    
    void set_chrom_order(const std::string order [], const unsigned int &size);
    std::vector<int> get_field_order(const std::string order [], const unsigned int &size);
    
    void set_version(const std::string &s);
    void set_sorted(const std::string &s);
    void set_shape(const std::string &s);
    void set_genome_assembly(const std::string &s);
    void set_command(const std::string &s);
    void set_columns(const std::string columns [], const unsigned int &size);
    void set_chromosomes(const std::string chrom [], const unsigned int length [], const unsigned int &size);
    
    void add_column(const std::string &c);
    void add_chromosome(const std::string &c, const unsigned int length);
    
    void sort_chromosome();
    
    PairsFileHeader();
    
    int get_chrom_id(const std::string &chrom);
    std::string get_chrom_by_id(const unsigned int &id);
    
    std::string Representation();
    std::string Representation(bool sort);
private:
    std::string version,sorted,shape,genome_assembly,command,misc;
    std::vector<std::string> columns;
    std::vector<Chromosome> chromosomes;
    std::map<std::string, int> chromMap;
    std::map<std::string, int> chromOrder;
};

#endif
