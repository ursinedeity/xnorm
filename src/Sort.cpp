/**
 * Copyright 2017 University of Southern California
 * Licensed under GPLv2 or any later version
 * Refer to the LICENSE file included.
 */

#include "Sort.h"

PairsFileSorter::PairsFileSorter(PairsFileHeader &header){
    header.sort_chromosome();
    this->chromosomes = header.chromosomes;
    this->chromMap = header.chromMap;
}
