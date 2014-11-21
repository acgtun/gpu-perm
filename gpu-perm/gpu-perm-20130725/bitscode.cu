#include "bitscode.h"

/*
 * uiReadLength must < WORD_SIZE
 * WORD_SIZE is 32 bp in 32 big machine and 64 bp in 64 bit machine
 * Each base is encoded into 2 bits: A -> 00, C->01, G->10 and T->11.
 * These two digits are located in two word, for bits operation.
 * The first nucleotide is encoded as the last digit.
 */

__device__ __host__ void EncodeRead(const char * strRead, InBits * readsInBits,
		int len) {

	readsInBits->ub = 0;
	readsInBits->lb = 0;

	// A 00
	// C 01
	// G 10
	// T 11

	for (int i = len - 1; i >= 0; i--) {
		if (strRead[i] == 'A' || strRead[i] == 'a') {

		} else if (strRead[i] == 'C' || strRead[i] == 'c') {
			readsInBits->lb++;
		} else if (strRead[i] == 'G' || strRead[i] == 'g') {
			readsInBits->ub++;
		} else if (strRead[i] == 'T' || strRead[i] == 't') {
			readsInBits->ub++;
			readsInBits->lb++;
		} else {
			//printf("Not A, C, G, T\n");
			//anything else as A
		}
		if (i != 0) {
			readsInBits->ub <<= 1; //left shift 1
			readsInBits->lb <<= 1;
		}
	}
}

void printWORD(WORD_SIZE word, SIZE_T len) {
	printf("\n");
	for (SIZE_T i = 0; i < len; i++) {
		if ((word & 0x01) == 1) {
			printf("1");
		} else {
			printf("0");
		}
		word >>= 1;
	}
	printf("\n");
}

void printWORD2File(FILE * fout, WORD_SIZE word, SIZE_T len) {
	for (SIZE_T i = 0; i < len; i++) {
		if ((word & 0x01) == 1) {
			fprintf(fout, "1");
		} else {
			fprintf(fout, "0");
		}
		word >>= 1;
	}
	fprintf(fout, "\n");
}

__device__ __host__ void DecodeRead(char * strReads, int readLen,
		const InBits * readsInBits) {
	WORD_SIZE UpperBits = readsInBits->ub;
	WORD_SIZE LowerBits = readsInBits->lb;
	int strReadsl = 0;
	for (int i = 0; i < readLen; i++) {
		WORD_SIZE c = (UpperBits & 0x01) << 1 | (LowerBits & 0x01);
		switch (c) {
		case 0x00:
			strReads[strReadsl++] = 'A';
			break;
		case 0x01:
			strReads[strReadsl++] = 'C';
			break;
		case 0x02:
			strReads[strReadsl++] = 'G';
			break;
		case 0x03:
			strReads[strReadsl++] = 'T';
			break;
		default:
			strReads[strReadsl++] = 'N';
			break;
		}
		LowerBits >>= 1;
		UpperBits >>= 1;
	}
	strReads[strReadsl] = 0;
}

__device__ __host__ void Swap(char * strVal, int len) {
	char chr;
	for (int i = 0; i < len / 2; i++) {
		chr = strVal[i];
		strVal[i] = strVal[len - i - 1];
		strVal[len - i - 1] = chr;
	}
}

__device__ __host__ void DecodeReadReverse(char * strRead, int readLen,
		const InBits * readsInBits) {
	DecodeRead(strRead, readLen, readsInBits);
	Swap(strRead, readLen);
}
