 
#ifndef __ALAB
#define __ALAB

#include <vector>
#include <string>
#include <map>
namespace alab{

class Index;
// Keep the index. All the variables are public since this is supposed 
// to be a mere container, and to keep syntax simple.
class Genome{
public:
    Genome(){}
    std::string assembly;
    std::vector<std::string> chroms;
    std::vector<unsigned int> origins;
    std::vector<unsigned int> lengths; 
    Genome(const std::string &assembly, 
           const std::vector<std::string> &chroms, 
           const std::vector<unsigned int> &origins, 
           const std::vector<unsigned int> &lengths);
    Genome(const std::string &assembly, 
           const std::vector<std::string> &chroms, 
           const std::vector<unsigned int> &lengths);
    std::string get_chrom(unsigned int c);
    int get_chromnum(const std::string & chrom);
    Index BinInfo(unsigned int resolution);
private:
    std::map<std::string, int> chromMap;
    
}; // struct genome

// A container of references, just in case one wants
// to use operator [] on index.
struct IndexItem{
    int& chrom;  // chromosome id
    int& copy;   // chromosome copy
    int& start;  // starting base pair
    int& end;    // ending base pair
    std::string& label;  // label for gap, domains, etc.
    IndexItem(int& ch, int& cp, int& st, int& en, std::string& lb):
        chrom(ch), copy(cp), start(st), end(en), label(lb) {};
};


// Keep the index. All the variables are public since this is supposed 
// to be a mere container, and to keep syntax simple.
// I'm not considering weird cases where records are inserted in a 
// non-consecutive way.
class Index{
public:
    Index(){}
    Index(const std::vector<unsigned int> &chromList, 
          const std::vector<unsigned int> &startList, 
          const std::vector<unsigned int> &endList,
          const std::vector<unsigned int> &chrom_sizes
         );

    std::vector<unsigned int> chrom;  // chromosome ids
    std::vector<unsigned int> copy;   // chromosome copies
    std::vector<unsigned int> start;  // starting basepairs
    std::vector<unsigned int> end;    // ending basepairs
    std::vector<std::string> label; // labels

    unsigned int num_chrom,size;  // number of chromosomes
    std::vector<unsigned int> chrom_sizes;  // sizes of chromosomes
    std::vector<unsigned int> offset;  // start position for each chromosome
    
    // returns an IndexItem, a container of references
    // to the values. Keep in mind this could break
    // badly if IndexItem is used after Index goes out of 
    // scope, for example.
    //IndexItem operator[](int i);
    
    // pushes back a record
    void push_back(unsigned int chrom, unsigned int copy, unsigned int start, unsigned int end, 
                   const std::string& label = "");

    // update num_chrom and the chrom_* vectors. Call if the data
    // vectors have been directly modified. 
    void rebuild_internal_index();

}; // class Index


class Matrix{
public:
    unsigned int *indptr;
    unsigned int *indices;
    float *data;
    unsigned int size, nnz;
    
    void LoadCoo(const unsigned int *Ai, const unsigned int *Aj, const float *Ax, unsigned int size, unsigned int nnz);
    
    ~Matrix(){
        delete [] indptr;
        delete [] indices;
        delete [] data;
    }
};
}; // namespace alab
#endif
