/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "Main.h"
#include "PairReads.h"
#include <unistd.h>

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
    
    if (std::string(argv[1]) == "merge"){ // Merge already sorted Hi-C pairs file and control file 
        Merge(argc,argv);
        return 0;
    }
    std::cerr << "Unrecognized command.\n";
    return 1;
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

void SortOut(int argc, char* argv[]){
    int optc;
    std::string outfile,corder,forder;
    std::vector<std::string> chromOrder,fieldOrder;
    while((optc = getopt(argc,argv,"o:c:f:")) != -1){
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
                fieldOrder = split(forder,'-');
                break;
            case '?':
                exit(1);
            default:
            break;
        }
    }
    
    // set default input from STDIN
    char* input = (char*)"-";
    
    // if there is an operand other than sort
    optind++;
    if (optind < argc){
        input = argv[optind++];
    }
    
    PairsFileHeader header;
    header.set_chrom_order(chromOrder.data(),chromOrder.size());
    //process pairs file
    if (input[0] != '-'){
        std::ifstream fin (input);
        if (!fin.is_open()){
            std::cerr << "Unable to open file " << input << '\n';
            exit(1);
        }
        
        for (std::string line; std::getline(fin,line);){
            if (line[0] == '#') header.ParseHeader(line);
            else ;
        }

    }else{
        // Detect if it's terminal
        // Throw an error and exit for terminal input
        if (isatty(fileno(stdin))){
            std::cerr << "No support for terminal input.\n";
            exit(1);
        }
        
        //important performance improvement
        //unable stdio sync
        std::ios_base::sync_with_stdio(false);
        
        for (std::string line; std::getline(std::cin,line);){
            if (line[0] == '#') header.ParseHeader(line);
            else ;
        }
    }
    //std::vector<int>  order = header.get_field_order(fieldOrder.data(),fieldOrder.size());
    
    std::cout << header.Representation(true);
    //for (auto it = order.begin(); it != order.end(); ++it) std::cout << *it;
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
    
    // parse arguments
    // l:inseart length default 1000bp
    // q:map quality cutoff default 30
    // c:control filename
    // p:pairlist filename
    // e:enzyme cutting site
    while((optc = getopt(argc,argv,"c:e:l:p:q:")) != -1){
        switch (optc) {
            case 'c':
                controlFile = optarg;
                break;
            case 'e':
                enzyme = optarg;
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
    
    // set default input from STDIN
    char* input = (char*)"-";
    
    // if there is an operand other than ps
    optind++;
    if (optind < argc){
        input = argv[optind++];
    }
    
    //process sam file
    if (input[0] != '-'){
        std::ifstream fin (input);
        if (!fin.is_open()){
            std::cerr << "Unable to open file " << input << '\n';
            exit(1);
        }
        
        PairReads worker(pairlistFile,controlFile,enzyme,insertLength,mapqCutoff);
        worker.PrintArgs();
        
        for (std::string line; std::getline(fin,line);){
            if (line[0] == '@') worker.AddHeader(line);
            else worker.AddAlignment(line);
        }
        worker.PrintStats();
    }else{
        // Detect if it's terminal
        // Throw an error and exit for terminal input
        if (isatty(fileno(stdin))){
            std::cerr << "No support for terminal input.\n";
            exit(1);
        }
        
        PairReads worker(pairlistFile,controlFile,enzyme,insertLength,mapqCutoff);
        worker.PrintArgs();
        //important performance improvement
        //unable stdio sync
        std::ios_base::sync_with_stdio(false);
        
        for (std::string line; std::getline(std::cin,line);){
            if (line[0] == '@') worker.AddHeader(line);
            else worker.AddAlignment(line);
        }
        worker.PrintStats();
    }
    
}
