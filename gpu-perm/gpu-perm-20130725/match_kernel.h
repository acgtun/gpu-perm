#ifndef MATCHING_H_
#define MATCHING_H_

#include "hash.h"
#include "refin.h"
#include "option.h"
#include "stdafx.h"
#include "seed.h"

typedef struct {
	SIZE_T nMismatch[100];
	SIZE_T nStartPos[100];
} ResultMatchedReads;

typedef struct {
	char * strReads;
	InBits * readsInBits;
	SIZE_T * readsLenAretN;
	int nReads;
} ReadsMatch;

__global__ void CUDAMapping(const MatchOpt * d_matchOpt,
		const RefGenome * d_refGenome, const HashTable * d_hashTable,
		const ReadsMatch * d_readsMatch, ResultMatchedReads * d_result);

#endif /* MATCHING_H_ */
