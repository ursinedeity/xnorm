/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "Main.h"
#include "PairReads.h"
#include "Sort.h"
#include "BuildMatrix.h"
#include <unistd.h>
#include <memory>

void ParseSam(int argc, char* argv[]);
void SortOut(int argc, char* argv[]);
void Merge(int argc, char* argv[]);

int main(int argc, char* argv[]){
    if (argc < 2){
        std::cerr << "Unrecognized command.\n";
        return 1;
    }
    
    if (std::string(argv[1]) == "ps"){ // parse sam file to extract Hi-C pairs and control
        ParseSam(argc,argv);
        return 0;
    }
    
    if (std::string(argv[1]) == "sort"){ // sort Hi-C pairs file and control file 
        SortOut(argc,argv);
        return 0;
    }
    
    if (std::string(argv[1]) == "bm"){ // Merge already sorted Hi-C pairs file and control file 
        BuildMatrix(argc,argv);
        return 0;
    }
    
    if (std::string(argv[1]) == "merge"){ // Merge already sorted Hi-C pairs file and control file 
        Merge(argc,argv);
        return 0;
    }
    std::cerr << "Unrecognized command.\n";
    return 1;
}

void BuildMatrix(int argc, char* argv[]){
    if (argc < 2){
        std::cerr << "Unrecognized command.\n";
        exit(1);
    }else
        BuildMatrixFromPairs(argv[2]);
}
/******************************************************************/

std::vector<std::string> split(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  return tokens;
}

std::shared_ptr<std::istream> GetInput(const char* inputFile){
    std::shared_ptr<std::istream> input;
    if (inputFile[0] == '-'){
        if (isatty(fileno(stdin))){
             std::cerr << "No support for terminal input.\n";
             exit(1);
        }
        std::ios_base::sync_with_stdio(false);
        input.reset(&std::cin, [](...){});
        
    }else
        input.reset(new std::ifstream(inputFile));
    
    return input;    
}

std::shared_ptr<std::ostream> GetOutput(const char* outputFile){
    std::shared_ptr<std::ostream> output;
    if (outputFile[0] == '-'){
        output.reset(&std::cout, [](...){});
    }else
        output.reset(new std::ofstream(outputFile));
    
    return output;    
}

/*****************************************************************/
void SortOut(int argc, char* argv[]){
    int optc;
    unsigned int threads=1;
    std::string outfile("-");
    std::string corder,genome_assembly;
    std::string forder("chr1-chr2-pos1-pos2");
    std::string shape("upper triangle");
    std::vector<std::string> chromOrder,fieldOrder;
    
    while((optc = getopt(argc,argv,"o:c:f:g:s:t:")) != -1){
        switch (optc) {
            case 'o':
                outfile = optarg;
                break;
            case 'c':
                corder = optarg;
                chromOrder = split(corder,'-');
                break;
            case 'f':
                forder = optarg;
                break;
            case 'g':
                genome_assembly = optarg;
                break;
            case 's':
                if (optarg[0] == 'l') shape = "lower triangle";
                break;
            case 't':
                threads = std::stoi(optarg);
                break;
            case '?':
                exit(1);
            default:
            break;
        }
    }
    fieldOrder = split(forder,'-');
    
    // set default input from STDIN
    char* inputFile = (char*)"-";
    
    // if there is an operand other than sort
    optind++;
    if (optind < argc){
        inputFile = argv[optind++];
    }
    
    std::shared_ptr<std::istream> input = GetInput(inputFile);
    
    PairsFileHeader header;
    
    
    //using arguments to get chromosome order e.g. -c chr1-chr2-chr3-chr4-...
    //unspecified chromosome in order will be discared.
    header.set_chrom_order(chromOrder.data(),chromOrder.size());
    header.set_shape(shape);
    
    //process pairs file
    std::string line;
    for (;std::getline(*input,line) && line[0] == '#';){
        header.ParseHeader(line);
    }
    header.set_sorted(forder);
    if (!genome_assembly.empty()) header.set_genome_assembly(genome_assembly);
    
    //sorter must have a header to initialize.
    PairsFileSorter sorter(header);
    sorter.AddRecord(line);
    for (;std::getline(*input,line);){
        sorter.AddRecord(line);
    }
    
    std::vector<int>  order = header.get_field_order(fieldOrder.data(),fieldOrder.size());
//     for (auto it = order.begin(); it != order.end(); ++it) std::cerr << *it << ' ';
//     std::cerr << std::endl;
    
    sorter.Sort(order.data(),order.size(),threads);

    std::shared_ptr<std::ostream> output = GetOutput(outfile.c_str());
    sorter.PrintRecords(output);
}

void Merge(int argc, char* argv[]){

    
}
/******************************************************************/
void ParseSam(int argc, char* argv[]){
    int optc;
    unsigned int insertLength = 1000;
    unsigned int mapqCutoff = 30;
    std::string pairlistFile("output.pairs");
    std::string controlFile("output.ctl");
    std::string enzyme("^GATC");
    std::string genome_assembly;
    // parse arguments
    // l:inseart length default 1000bp
    // q:map quality cutoff default 30
    // c:control filename
    // p:pairlist filename
    // e:enzyme cutting site
    while((optc = getopt(argc,argv,"c:e:g:l:p:q:")) != -1){
        switch (optc) {
            case 'c':
                controlFile = optarg;
                break;
            case 'e':
                enzyme = optarg;
                break;
            case 'g':
                genome_assembly = optarg;
                break;
            case 'l':
                insertLength = std::stoi(optarg);
                break;
            case 'p':
                pairlistFile = optarg;
                break;
            case 'q':
                mapqCutoff = std::stoi(optarg);
                break;
            case '?':
                exit(1);
            default:
            break;
        }
    }
    /*******************Handling Input**********************/
    // set default input from STDIN
    char* inputFile = (char*)"-";
    
    // if there is an operand other than ps
    optind++;
    if (optind < argc){
        inputFile = argv[optind++];
    }
    
    //process sam file
    
    std::shared_ptr<std::istream> input = GetInput(inputFile);
    
    PairReads worker(pairlistFile,controlFile,enzyme,insertLength,mapqCutoff);
    worker.PrintArgs();
    
    if (!genome_assembly.empty()) worker.set_genome_assembly(genome_assembly);
    for (std::string line; std::getline(*input,line);){
        if (line[0] == '@') worker.AddHeader(line);
        else worker.AddAlignment(line);
    }
    worker.PrintStats();
    
}
