/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "PairsFileHeader.h"
#include <algorithm>
#include <sstream>

int PairsFileHeader::ParseHeader(const std::string &line){
    std::istringstream linestream(line);
    std::string f,key,tmp;
    linestream >> key;
    if (key == "##"){
        set_version(line.substr(16));
        return 0;
    }
    
    if (key == "#sorted:"){
        set_sorted(line.substr(9));
        return 0;
    }
    
    if (key == "#shape:"){
        set_shape(line.substr(8));
        return 0;
    }
    
    if (key == "#genome_assembly:"){
        set_genome_assembly(line.substr(18));
        return 0;
    }
    
    if (key == "#command:"){
        set_command(line.substr(10));
        return 0;
    }
    
    if (key == "#chromsize:"){
        linestream >> f;
        int length;
        linestream >> length;
        add_chromosome(f,length);
        return 0;
    }
    
    if (key == "#columns:"){
        while (linestream >> f) add_column(f);
        return 0;
    }
    misc += line + "\n";
    return 1;
}





PairsFileHeader::PairsFileHeader(){
    set_version("v1.0");
}

void PairsFileHeader::set_chrom_order(const std::string order[], const unsigned int &size){
    chromOrder.clear();
    for (unsigned int i=0; i<size; ++i) chromOrder[order[i]] = i;    
}

std::vector<int> PairsFileHeader::get_field_order(const std::string forder [], const unsigned int &size){
    std::vector<int> order;
    for (unsigned int i=0; i<size; ++i){
        auto p = std::find(columns.begin(),columns.end(),forder[i]);
        if (p != columns.end())
            order.push_back(p - columns.begin());
        else
            order.push_back(-1);
    }
    return order;
}


std::string PairsFileHeader::Representation(bool sort){
    std::string repr = version;
    repr += sorted + shape + genome_assembly;
    
    if (chromosomes.size() == 0) std::cerr << "chromosome sizes not defined, header might be invalid.\n";
    
    if (sort) sort_chromosome();
    
    for (auto it = chromosomes.begin(); it != chromosomes.end(); ++it)
        repr += "#chromsize: " + it->chrom + " " + std::to_string(it->length) + "\n";
    
    repr += "#columns:";
    for (auto it = columns.begin(); it != columns.end(); ++it)
        repr += " " + *it;
    
    repr += "\n" + misc;
    return repr;
}

std::string PairsFileHeader::Representation(){
    return Representation(false);
}






int PairsFileHeader::get_chrom_id(const std::string &chrom){
    auto it = chromMap.find(chrom);
    if (it != chromMap.end())
        return it->second;
    else
        return -1;
}

std::string PairsFileHeader::get_chrom_by_id(const unsigned int &id){
    if (id < chromosomes.size())
        return chromosomes[id].chrom;
    else
        return "";
}








void PairsFileHeader::set_version(const std::string &s){
    version = "## pairs format " + s + "\n";    
}

void PairsFileHeader::set_sorted(const std::string &s){
    sorted = "#sorted: " + s + "\n";
}

void PairsFileHeader::set_shape(const std::string &s){
    shape = "#shape: " + s + "\n";
}

void PairsFileHeader::set_genome_assembly(const std::string &s){
    genome_assembly = "#genome_assembly: " + s + "\n";
}

void PairsFileHeader::set_command(const std::string &s){
    command = "#command:" + s + "\n";
}
void PairsFileHeader::set_columns(const std::string columns [], const unsigned int &size){
    for (unsigned int i=0; i<size; ++i) add_column(columns[i]);
}

void PairsFileHeader::set_chromosomes(const std::string chrom [], const unsigned int length [], const unsigned int &size){
    for (unsigned int i=0; i<size; ++i) add_chromosome(chrom[i], length[i]);
}

void PairsFileHeader::add_column(const std::string &c){
    columns.push_back(c);
}

void PairsFileHeader::add_chromosome(const std::string &c, const unsigned int length){
    if (chromOrder.empty()){
        chromosomes.push_back(Chromosome(c,length));
        chromMap[c] = chromosomes.size()-1;
    }else{
        auto it = chromOrder.find(c);
        if (it != chromOrder.end()){
            chromosomes.push_back(Chromosome(c, length, it->second));
            chromMap[c] = chromosomes.size()-1;
        }
    }
}





void PairsFileHeader::sort_chromosome(){
    std::sort (chromosomes.begin(), chromosomes.end());
    for (unsigned int i = 0; i < chromosomes.size(); ++i) chromMap[chromosomes[i].chrom] = i;
}

alab::Genome PairsFileHeader::MakeGenome(){
    sort_chromosome();
    std::vector<std::string> chroms;
    std::vector<unsigned int> lengths;
    for (auto &chr : chromosomes){
        chroms.push_back(chr.chrom);
        lengths.push_back(chr.length);
    }
    
    return alab::Genome(genome_assembly, chroms, lengths);
}

/********************************************************************/
// int main(){
//     PairsFileHeader p;
//     p.set_genome_assembly("hg38");
//     std::string  chrom []= {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
//                             "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
//                             "chrM","chrX","chrY"};
//     int length [] = {195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,16299,171031299,91744698};
//     std::string order [] = {"chrX","chr2","chr3","chr4","chr5"};
//     //p.set_chrom_order(order,5);
//     p.set_chromosomes(chrom,length,22);
//     std::cout << p.Representation(true) << std::endl;
//     std::cout << p.get_chrom_id("chr5") << " " << p.get_chrom_by_id(5) << std::endl;
// }
