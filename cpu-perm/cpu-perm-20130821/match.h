#pragma once
#ifndef MATCHING_H_
#define MATCHING_H_

#include "hash.h"
#include "refin.h"
#include "option.h"
#include "stdafx.h"
#include "seed.h"

typedef struct {
	char readInStr[MAX_READ_LEN];
	SIZE_T readLen;
} CReadInStr;

//typedef struct {
//	InBits readInBits[5];
//	SIZE_T readLen;
//} CReadInBits;

typedef struct {
	SIZE_T lower;
	SIZE_T upper;
} CRegion;

typedef struct {
	SIZE_T pos;
	SIZE_T readid;
	char org_rev;
} CMatch;

void Match(const CReference * refGenome, const CHashTable * hashTable, const CReadInStr * oneRead, const SIZE_T & nStartPos, CRegion * oneResult);

#endif /* MATCHING_H_ */
