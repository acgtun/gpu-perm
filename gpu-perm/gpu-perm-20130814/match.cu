#include "match.h"
#include "seed.h"

__device__ int CMP(const CReadInStr * oneRead, SIZE_T index, const CReference * refGenome) {

	InBits ref; //+22
	SIZE_T indexInWords = index / wordSize;
	SIZE_T bitsShift = index % wordSize;
	ref.ub = refGenome->refInBits[indexInWords].ub >> bitsShift;
	ref.lb = refGenome->refInBits[indexInWords].lb >> bitsShift;

	SIZE_T cmpbits = 0;
	WORD_SIZE code = 0x17; //0010111
	int bits2 = wordSize - bitsShift;

	while (cmpbits < oneRead->readLen) {
		if (cmpbits % 7 == 0)
			code = 0x17;
		if (code & 0x01) {
			char c1 = oneRead->readInStr[cmpbits];
			//SIZE_T c2 = (ref.ub & 0x01) << 1 | (ref.lb & 0x01);
			char c2 = getNT((ref.ub & 0x01) << 1 | (ref.lb & 0x01));
			//printf("c2 = %c\n");

			if (c1 > c2) {
				return 1;
			} else if (c1 < c2) {
				return -1;
			}
		}
		code >>= 1;
		ref.ub >>= 1;
		ref.lb >>= 1;
		cmpbits++;
		bits2--;

		if (index + cmpbits >= refGenome->nRefSize)
			return 1;
		if (cmpbits >= oneRead->readLen)
			return 0; //if len2 > len1, the string equals to certain part of genome

		if (bits2 == 0) {
			indexInWords++;
			ref = refGenome->refInBits[indexInWords];
			if (index + cmpbits + wordSize < refGenome->nRefSize)
				bits2 = wordSize;
			else
				bits2 = refGenome->nRefSize - cmpbits - index;
		}
	}
	return 0;
}

__device__ SIZE_T LowerBound(SIZE_T low, SIZE_T high, const CReadInStr * oneRead, const CReference * refGenome,
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

__device__ SIZE_T UpperBound(SIZE_T low, SIZE_T high, const CReadInStr * oneRead, const CReference * refGenome,
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

__device__ void Match(const CReference * refGenome, const CHashTable * hashTable, const CReadInStr * oneRead, CRegion * oneResult) {

	SIZE_T hashValue = GetHashValue(oneRead->readInStr);

	SIZE_T l = hashTable->counter[hashValue];
	SIZE_T u = hashTable->counter[hashValue + 1] - 1;

	oneResult->lower = LowerBound(l, u, oneRead, refGenome, hashTable);
	oneResult->upper = UpperBound(l, u, oneRead, refGenome, hashTable);
}
