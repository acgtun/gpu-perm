#include "match_kernel.h"
#include "seed.h"

__device__ InBits GetShiftRead(SIZE_T shift, InBits read) {
	InBits r;
	r.ub = read.ub >> shift;
	r.lb = read.lb >> shift;
	return r;
}

__device__ int cuSTRCMP(char * str1, int len1, char * str2, int len2) {
	int i = 0;
	while (i < len1 && i < len2) {
		if (str1[i] < str2[i])
			return -1;
		else if (str1[i] > str2[i])
			return 1;
		i++;
	}
	if (i == len1 && i == len2)
		return 0;
	if (i == len2)
		return 1;
	else
		return -1;
}

__device__ int CMP(char * strRead, int len, SIZE_T index,
		const RefGenome * d_refGenome) {
	InBits r, v;
	SIZE_T s = GetKmer(d_refGenome, index, wordSize, &r);
	int ss = GetF2SeedForBits(r, s, &v);

	char strRef64[MAX_READ_LEN];
	DecodeReadReverse(strRef64, ss, &v);

	return cuSTRCMP(strRead, len, strRef64, ss);
}

__device__ SIZE_T bitsSetCount(WORD_SIZE bits) {
	/* magic function to calculate how many ones are there */
	SIZE_T c; // c accumulates the total bits set in v
	for (c = 0; bits; c++) {
		bits &= bits - 1; // clear the least significant bit set
	}
	return c;
}

__device__ SIZE_T bitsStrNCompare(InBits r1, InBits r2, SIZE_T N) {
	/* compare only the last N bases (bits) */
	WORD_SIZE bits = (r1.ub ^ r2.ub) | (r1.lb ^ r2.lb);
	bits <<= (wordSize - N);
	return ((SIZE_T) bitsSetCount(bits));
}

__device__ SIZE_T LowerBound(SIZE_T low, SIZE_T high, char * strRead, int s,
		const RefGenome * d_refGenome, const HashTable * d_hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high) / 2;
		if (CMP(strRead, s, d_hashTable->index[mid], d_refGenome) <= 0)
			high = mid;
		else
			low = mid + 1;
	}
	return low;
}

__device__ SIZE_T UpperBound(SIZE_T low, SIZE_T high, char * strRead, int s,
		const RefGenome * d_refGenome, const HashTable * d_hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high + 1) / 2;
		if (CMP(strRead, s, d_hashTable->index[mid], d_refGenome) >= 0)
			low = mid;
		else
			high = mid - 1;
	}
	return low;
}

__global__ void CUDAMapping(const MatchOpt * d_matchOpt,
		const RefGenome * d_refGenome, const HashTable * d_hashTable,
		const ReadsMatch *d_readsMatch, ResultMatchedReads * d_result) {
	SIZE_T readID = threadIdx.x + blockDim.x * blockIdx.x;
	if (readID >= d_readsMatch->nReads)
		return;

	int readLen = d_readsMatch->readsLenAretN[readID];
	int nRet = 0;
	for (int i = 0; i <= NUMBER_OF_SHIFT; i++) {
		InBits read = GetShiftRead(i, d_readsMatch->readsInBits[readID]);
		SIZE_T hashValue = GetHashValue(read);
		InBits ret;
		int len = readLen - i;
		int s = GetF2SeedForBits(read, len, &ret);
		char strRead[MAX_READ_LEN];
		DecodeReadReverse(strRead, s, &ret);

		SIZE_T l = d_hashTable->counter[hashValue];
		SIZE_T u = d_hashTable->counter[hashValue + 1] - 1;

		SIZE_T lower = LowerBound(l, u, strRead, s, d_refGenome, d_hashTable);
		SIZE_T upper = UpperBound(l, u, strRead, s, d_refGenome, d_hashTable);
		printf("lower = %d upper = %d\n", lower, upper);
		for (SIZE_T j = lower; j <= upper; j++) {
			SIZE_T s = GetKmer(d_refGenome, d_hashTable->index[j] - i, wordSize,
					&ret);
			SIZE_T nDiff = bitsStrNCompare(ret,
					d_readsMatch->readsInBits[readID], readLen);
			if (nDiff <= d_matchOpt->nMaxMismatch) {
				d_result[readID].nMismatch[nRet] = nDiff;
				d_result[readID].nStartPos[nRet] = d_hashTable->index[j] - i;
				printf("readID = %d %d %d\n", d_hashTable->index[j], i);
				nRet++;
			}
		}
	}
	d_readsMatch->readsLenAretN[readID] = nRet;
}
