#ifndef MATCHING_H_
#define MATCHING_H_

#include "hash.h"
#include "refin.h"
#include "option.h"
#include "stdafx.h"
#include "seed.h"

typedef struct {
	char * strReads;
	InBits * readsInBits;
	SIZE_T * readsLen;
	SIZE_T nReads;
} ReadsMatch;

__device__ void MappingOneRead(int nMaxMismatch, const InBits * refGenome, SIZE_T RefSize, const SIZE_T * hashCounter, const SIZE_T * hashIndex,
		InBits readInBits, int readLen, ResultMatchedReads * oneResult);

#endif /* MATCHING_H_ */
