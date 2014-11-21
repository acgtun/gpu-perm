#include "match.h"
#include "seed.h"

__device__ int CMP(const CReadInBits * oneRead, SIZE_T index, const CReference * refGenome) {

	InBits ref; //+22
	SIZE_T indexInWords = index / wordSize;
	SIZE_T bitsShift = index % wordSize;
	ref.ub = refGenome->refInBits[indexInWords].ub >> bitsShift;
	ref.lb = refGenome->refInBits[indexInWords].lb >> bitsShift;

	SIZE_T indexInWordsread = 0;
	InBits read = oneRead->readInBits[indexInWordsread];

	SIZE_T cmpbits = 0;

	WORD_SIZE code = 0x17; //0010111

	int bits1;
	if (oneRead->readLen >= wordSize)
		bits1 = wordSize;
	else
		bits1 = oneRead->readLen;

	int bits2 = wordSize - bitsShift;

	while (cmpbits < 1000) {
		if (cmpbits % 7 == 0)
			code = 0x17;
		if (code & 0x01) {
			SIZE_T c1 = (read.ub & 0x01) << 1 | (read.lb & 0x01);
			SIZE_T c2 = (ref.ub & 0x01) << 1 | (ref.lb & 0x01);

			if (c1 > c2) {
				return 1;
			} else if (c1 < c2) {
				return -1;
			}
		}
		code >>= 1;
		read.ub >>= 1;
		read.lb >>= 1;
		ref.ub >>= 1;
		ref.lb >>= 1;
		cmpbits++;
		bits1--;
		bits2--;

		if (index + cmpbits >= refGenome->nRefSize)
			return 1;
		if (cmpbits >= oneRead->readLen)
			return 0; //if len2 > len1, the string equals to certain part of genome

		if (bits1 == 0) {
			indexInWords++;
			ref = refGenome->refInBits[indexInWords];
			if (index + cmpbits + wordSize < refGenome->nRefSize)
				bits1 = wordSize;
			else
				bits1 = refGenome->nRefSize - cmpbits - index;
		}

		if (bits2 == 0) {
			indexInWordsread++;
			read = oneRead->readInBits[indexInWordsread];
			if (cmpbits + wordSize < oneRead->readLen)
				bits2 = wordSize;
			else
				bits2 = oneRead->readLen - cmpbits;
		}
	}
	return 0;
}

__device__ SIZE_T LowerBound(SIZE_T low, SIZE_T high, const CReadInBits * oneRead, const CReference * refGenome,
		const CHashTable * hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high) / 2;
		if (CMP(oneRead, hashTable->index[mid], refGenome) <= 0)
			high = mid;
		else
			low = mid + 1;
	}
	return low;
}

__device__ SIZE_T UpperBound(SIZE_T low, SIZE_T high, const CReadInBits * oneRead, const CReference * refGenome,
		const CHashTable * hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high + 1) / 2;
		if (CMP(oneRead, hashTable->index[mid], refGenome) >= 0)
			low = mid;
		else
			high = mid - 1;
	}
	return low;
}

__device__ void Match(const CReference * refGenome, const CHashTable * hashTable, const CReadInBits * oneRead, CRegion * oneResult) {

	SIZE_T hashValue = GetHashValue(oneRead->readInBits[0]);

	SIZE_T l = hashTable->counter[hashValue];
	SIZE_T u = hashTable->counter[hashValue + 1] - 1;

	SIZE_T lower = LowerBound(l, u, oneRead, refGenome, hashTable);
	SIZE_T upper = UpperBound(l, u, oneRead, refGenome, hashTable);
}
