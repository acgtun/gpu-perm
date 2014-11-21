#pragma once
#ifndef HASH_H_
#define HASH_H_

#include "stdafx.h"
#include "refin.h"
#include "bitscode.h"
#include "seed.h"

//typedef struct {
//	SIZE_T * counter;
//	SIZE_T * index;
//} HashTable;

void MakeHashTable(const InBits * refGenome, SIZE_T nRefSize,
		SIZE_T ** hashCounter, SIZE_T ** hashIndex);

#endif /* HASH_H_ */
