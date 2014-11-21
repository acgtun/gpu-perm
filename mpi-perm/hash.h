#pragma once
#ifndef HASH_H_
#define HASH_H_

#include "sdk.h"
#include "option.h"
#include "iofile.h"

typedef struct {
  uint32_t * counter;
  uint32_t * index;
  uint32_t nSizeCounter;
  uint32_t nSizeIndex;
} CHashTable;

typedef struct {
  char * refSeq;
  uint32_t nRefSize;
} CReference;

void BuildIndex(const Option & opt, CReference * refGenome,
                CHashTable * hashTable);
void ReadIndexAndRef(const Option & opt, CReference * refGenome,
                     CHashTable * hashTable);

#endif /* HASH_H_ */
