#ifndef MATCHING_H_
#define MATCHING_H_

#include "hash.h"
#include "refin.h"
#include "option.h"
#include "stdafx.h"
#include "seed.h"

typedef struct {
	InBits readInBits;
	SIZE_T readLen;
} CRead;

typedef struct {
	CRead * reads;
	SIZE_T nReadsNum;
} CReadArray;

__device__ void MappingOneRead(const MapOpt mapOpt, const CReference * refGenome, const CHashTable * hashTable, CRead read, CResult * oneResult);

#endif /* MATCHING_H_ */
