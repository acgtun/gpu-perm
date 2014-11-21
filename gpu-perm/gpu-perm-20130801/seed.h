#pragma once
#ifndef SEED_H_
#define SEED_H_
#include "stdafx.h"
#include "refin.h"

__host__ __device__  SIZE_T GetHashValue(InBits r);
__host__ __device__  SIZE_T GetKmer(const InBits * refGenome, SIZE_T nRefSize, SIZE_T nRefStart, SIZE_T kmerLen, InBits * r);
__host__ __device__  int GetF2SeedForBits(InBits r, SIZE_T len, InBits * ret);

#endif /* SEED_H_ */
