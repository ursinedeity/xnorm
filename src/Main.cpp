/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "Main.h"
#include "PairReads.h"
#include <unistd.h>

void ParseSam(int argc, char* argv[]){
    int optc;
    int insertLength = 1000;
    int mapqCutoff = 30;
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

int main(int argc, char* argv[]){
    if (argc < 2){
        std::cerr << "Unrecognized command.\n";
        return 1;
    }
    
    if (std::string(argv[1]) == "ps"){
        ParseSam(argc,argv);
        return 0;
    }
    
    std::cerr << "Unrecognized command.\n";
    return 1;
}

/******************************************************************/
