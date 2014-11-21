#pragma once
#ifndef SEED_H_
#define SEED_H_
#include "stdafx.h"
#include "refin.h"

SIZE_T GetHashValue(InBits r);
SIZE_T GetKmer(const RefGenome * refGenome, SIZE_T nRefStart, SIZE_T kmerLen,
		InBits * r);
int GetF2SeedForBits(InBits r, SIZE_T len, InBits * ret);

#endif /* SEED_H_ */
