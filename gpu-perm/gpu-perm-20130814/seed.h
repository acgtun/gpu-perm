#pragma once
#ifndef SEED_H_
#define SEED_H_
#include "stdafx.h"
#include "refin.h"

__device__ SIZE_T GetHashValue(const char * strVal);
__host__ __device__ SIZE_T GetHashValue(InBits r);
__host__ __device__ SIZE_T GetKmer(const CReference * refGenome, SIZE_T nRefStart, SIZE_T kmerLen, InBits * r);
__host__ __device__ int GetF2SeedForBits(InBits r, SIZE_T len, InBits * ret);

#endif /* SEED_H_ */
