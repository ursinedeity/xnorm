#include "Alab.h"
#include <cmath>
#include <iostream>
namespace alab{
Genome::Genome(const std::string &assembly, const std::vector<std::string> &chroms, const std::vector<unsigned int> &origins, const std::vector<unsigned int> &lengths){
    if (chroms.size() != lengths.size()){
        std::cerr << "Chroms and lengths sizes don't match!\n";
        exit(1);        
    }
    
    this->assembly = assembly;
    this->chroms = chroms;
    this->origins = origins;
    this->lengths = lengths;
    for (unsigned int i = 0; i < chroms.size(); ++i)
        chromMap[chroms[i]] = i;
}
Genome::Genome(const std::string &assembly, const std::vector<std::string> &chroms, const std::vector<unsigned int> &lengths){
    this->assembly = assembly;
    this->chroms = chroms;
    this->lengths = lengths;
    
    for (unsigned int i = 0; i < chroms.size(); ++i){
        this->origins.push_back(0);
        chromMap[chroms[i]] = i;
    }
}

std::string Genome::get_chrom(unsigned int c){
    return chroms[c];
}

int Genome::get_chromnum(const std::string & chrom){
    auto it = chromMap.find(chrom);
    if (it != chromMap.end())
        return it->second;
    else
        return -1;
}

Index Genome::BinInfo(unsigned int resolution){
    std::vector<unsigned int> binSize, chromList, binLabel, startList, endList;
    
    for (auto &l : lengths)
        binSize.push_back( l/resolution + ( l%resolution != 0) );
    
    for (unsigned int i = 0; i < chroms.size(); ++i){
        for (unsigned int j = 0; j < binSize[i]; ++j){
            chromList.push_back(i);
            binLabel.push_back(j + origins[i]/resolution);
        }
    }
    
    
    for (unsigned int &label : binLabel){
        startList.push_back (label*resolution );
        endList.push_back( (label+1)*resolution );
    }
    
    return Index(chromList, startList, endList, binSize);
}

Index::Index(const std::vector<unsigned int> &chromList, 
             const std::vector<unsigned int> &startList, 
             const std::vector<unsigned int> &endList, 
             const std::vector<unsigned int> &chrom_sizes){
    this->chrom = chromList;
    this->start = startList;
    this->end = endList;
    num_chrom = chromList.size();
    
    for (unsigned int i = 0; i < chromList.size(); ++i){
        copy.push_back(0);
        label.push_back("");
    }
    
    this->chrom_sizes = chrom_sizes;
    offset.push_back(0);
    for (unsigned int i = 0, sum = 0; i < chrom_sizes.size(); ++i){
        sum += chrom_sizes[i];
        offset.push_back(sum);
        size = sum;
    }    
}

// IndexItem Index::operator[](int i){
//   return IndexItem(chrom[i], copy[i], start[i], end[i], label[i]);
// }


void Matrix::LoadCoo(const unsigned int *Ai, const unsigned int *Aj, const float *Ax, unsigned int size, unsigned int nnz){
    indptr = new unsigned int [size+1];
    indices = new unsigned int [nnz];
    data = new float [nnz];
    this->size = size;
    this->nnz = nnz;
    
    std::fill(indptr, indptr + size, 0);
    
    for (unsigned int i = 0; i < nnz; ++i)
        indptr[Ai[i]]++;
    
    for (unsigned int i = 0, cumsum = 0; i < size; ++i){
        unsigned int tmp = indptr[i];
        indptr[i] = cumsum;
        cumsum += tmp;
    }
    indptr[size] = nnz;
    
    for (unsigned int i = 0; i < nnz; ++i){
        unsigned int row = Ai[i];
        unsigned int dest = indptr[row];
        
        indices[dest] = Aj[i];
        data[dest] = Ax[i];
        
        indptr[row]++;
    }
    
    for (unsigned int i = 0, last = 0; i <= size; ++i){
        unsigned int tmp = indptr[i];
        indptr[i] = last;
        last = tmp;
    }
    
}

};
