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

void EncodeRead(const char * strReads, InBits * readsInBits,
		int len);
void DecodeRead(char * strReads, int readLen,
		const InBits * readsInBits);
void DecodeReadReverse(char * strReads, int readLen,
		const InBits * readsInBits);

void printWORD(WORD_SIZE word, SIZE_T len);
void printWORD2File(FILE * fout, WORD_SIZE word, SIZE_T len);

#endif /* BITSCODE_H_ */
