#pragma once
#ifndef MATCHING_H_
#define MATCHING_H_

#include "hash.h"
#include "refin.h"
#include "option.h"
#include "stdafx.h"
#include "seed.h"

typedef struct {
	char readInStr[320];
	SIZE_T readLen;
} CReadInStr;

typedef struct {
	InBits readInBits[5];
	SIZE_T readLen;
} CReadInBits;

typedef struct {
	SIZE_T lower;
	SIZE_T upper;
} CRegion;

__device__ void Match(const CReference * refGenome, const CHashTable * hashTable, const CReadInBits * oneRead, CRegion * oneResult);
__device__ void reverseCompliment(char * strRead, char * strReadO, int len);

#endif /* MATCHING_H_ */
