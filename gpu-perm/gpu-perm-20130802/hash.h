#pragma once
#ifndef HASH_H_
#define HASH_H_

#include "stdafx.h"
#include "refin.h"
#include "bitscode.h"
#include "seed.h"
#include "common.h"

typedef struct {
	SIZE_T * counter;
	SIZE_T * index;
} CHashTable;

void MakeHashTable(const CReference * refGenome, CReference * d_refGenome, CHashTable * d_hashTable);

#endif /* HASH_H_ */
