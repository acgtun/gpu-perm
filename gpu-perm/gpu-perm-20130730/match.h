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

__host__ __device__ void MappingOneRead(MatchOpt matchOpt, const RefGenome * refGenome,
		const HashTable * hashTable, InBits readInBits, int readLen,
		ResultMatchedReads * oneResult);

#endif /* MATCHING_H_ */
