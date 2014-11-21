#pragma once
#ifndef BITSCODE_H_
#define BITSCODE_H_

#include "stdafx.h"

typedef struct {
	WORD_SIZE ub;
	WORD_SIZE lb;
} InBits;

//1000000000000000000000000000000000000000000000000000000000000000000000000000
//static const unsigned short getonebit[66];

__host__ __device__ void EncodeRead(const char * strReads, InBits * readsInBits, int len);
__host__ __device__ void DecodeRead(char * strReads, int readLen, const InBits * readsInBits);
__host__ __device__ void DecodeReadReverse(char * strReads, int readLen, const InBits * readsInBits);

__device__ void printWORD(WORD_SIZE word, SIZE_T len);
__device__ void printWORD2File(FILE * fout, WORD_SIZE word, SIZE_T len);

__host__ __device__ void Swap(char * strVal, int len);
__device__ char complimentBase(char nt);
__device__ SIZE_T bitsStrNCompare(InBits r1, InBits r2, SIZE_T N);

#endif /* BITSCODE_H_ */
