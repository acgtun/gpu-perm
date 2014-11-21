#pragma once
#ifndef HASH_H_
#define HASH_H_

#include "stdafx.h"
#include "refin.h"
#include "bitscode.h"
#include "seed.h"

typedef struct {
	SIZE_T * counter;
	SIZE_T * index;
	SIZE_T nSizeCounter;
	SIZE_T nSizeIndex;
} CHashTable;

void MakeHashTable(const CReference * refGenome, CHashTable * hashTable, Option & opt);
void ReadIndexAndRef(CReference * refGenome, CHashTable * hashTable, const Option & opt);

#endif /* HASH_H_ */
