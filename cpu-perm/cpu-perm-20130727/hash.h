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
	SIZE_T NO_OF_BUCKET;
} HashTable;

void MakeHashTable(const RefGenome * refGenome, HashTable * hashTable);
void CountBucketSize(const RefGenome * refGenome, HashTable * hashTable);
void HashToBucket(const RefGenome * refGenome, HashTable * hashTable);
void SortEachBucket(const RefGenome * refGenome, HashTable * hashTable);
int Compare(const void * a, const void * b);

#endif /* HASH_H_ */
