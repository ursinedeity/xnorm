/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "BuildMatrix.h"
#include "PairsFileHeader.h"
#include <string>

void BuildMatrixFromPairs(const char* pairsfile_gz){
    pairix_t *tb;
    
    tb = load_from_file(pairsfile_gz);
    
    const ti_conf_t *idxconf;
    idxconf = ti_get_conf(tb->idx);
    
    int len;
    const char * s;
    sequential_iter_t *siter = ti_query_general(tb, 0, 0, 0);
    while ((s = sequential_ti_read(siter, &len)) != 0) {
    if ((int)(*s) != idxconf->meta_char) break;
        std::cout << s << std::endl;
    }
    
    
    
    
    
    
    
    destroy_sequential_iter(siter);
    bgzf_seek(tb->fp, 0, SEEK_SET);
    sit = ti_query_2d_general(tb, "chr1",12000000,13000000,"chr1",12000000,13000000);
    
    int n = 0;
    while((s = sequential_ti_read(sit, &len)) != 0){
        n++;
    }
    destroy_sequential_iter(sit);
    std::cout << n << std::endl;
    n=0;
    sit = ti_query_2d_general(tb, "chr1",11000000,12000000,"chr1",11000000,12000000);
    while((s = sequential_ti_read(sit, &len)) != 0){
        n++;
    }
    destroy_sequential_iter(sit);
    std::cout << n << std::endl;
    ti_close(tb);
}
