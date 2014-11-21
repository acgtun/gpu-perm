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
	int nRet;
} ResultMatchedReads;

typedef struct {
	char * strReads;
	InBits * readsInBits;
	SIZE_T * readsLen;
	int nReads;
} ReadsMatch;

void MappingOneRead(const MatchOpt * matchOpt, const RefGenome * refGenome,
		const HashTable * hashTable, InBits readInBits, int readLen,
		ResultMatchedReads * oneResult);
void CPU_Matching(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable);
#endif /* MATCHING_H_ */
