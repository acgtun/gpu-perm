#include "match.h"
#include "seed.h"

__host__ __device__ InBits GetShiftRead(SIZE_T shift, InBits read) {
	InBits r;
	r.ub = read.ub >> shift;
	r.lb = read.lb >> shift;
	return r;
}

__host__ __device__ int STRCMP(const char * str1, int len1, const char * str2,
		int len2) {
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
		return 0; //if len2 > len1, the string equals to certain part of genome
}

__host__ __device__ int CMP(char * strRead, int len, SIZE_T index,
		const RefGenome * refGenome) {
	InBits r, v;
	SIZE_T s = GetKmer(refGenome, index, wordSize, &r);
	int ss = GetF2SeedForBits(r, s, &v);

	char strRef64[MAX_READ_LEN];
	DecodeReadReverse(strRef64, ss, &v);

	return STRCMP(strRead, len, strRef64, ss);
}

__host__ __device__ SIZE_T LowerBound(SIZE_T low, SIZE_T high, char * strRead,
		int s, const RefGenome * refGenome, const HashTable * hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high) / 2;
		if (CMP(strRead, s, hashTable->index[mid], refGenome) <= 0)
			high = mid;
		else
			low = mid + 1;
	}
	return low;
}

__host__ __device__ SIZE_T UpperBound(SIZE_T low, SIZE_T high, char * strRead,
		int s, const RefGenome * refGenome, const HashTable * hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high + 1) / 2;
		if (CMP(strRead, s, hashTable->index[mid], refGenome) >= 0)
			low = mid;
		else
			high = mid - 1;
	}
	return low;
}

__host__ __device__ void reverseCompliment(char * strRead, InBits readInBits,
		int len) {
	DecodeReadReverse(strRead, len, &readInBits);
	for (int i = 0; i < len; i++) {
		strRead[i] = complimentBase(strRead[i]);
	}
}

__host__ __device__ int ResultExist(const ResultMatchedReads * oneResult,
		SIZE_T index) {
	for (SIZE_T i = 0; i < oneResult->nRet; i++) {
		if (oneResult->nStartPos[i] == index)
			return 1;
	}
	return 0;
}

__host__ __device__ void Mapping(const MatchOpt * matchOpt,
		const RefGenome * refGenome, const HashTable * hashTable,
		InBits readInBits, int readLen, ResultMatchedReads * oneResult) {

	for (int i = 0; i <= NUMBER_OF_SHIFT; i++) {
		InBits read = GetShiftRead(i, readInBits);
		SIZE_T hashValue = GetHashValue(read);
		InBits ret;
		int len = readLen - i;
		int s = GetF2SeedForBits(read, len, &ret);
		char strRead[MAX_READ_LEN];
		DecodeReadReverse(strRead, s, &ret);

		SIZE_T l = hashTable->counter[hashValue];
		SIZE_T u = hashTable->counter[hashValue + 1] - 1;

		SIZE_T lower = LowerBound(l, u, strRead, s, refGenome, hashTable);
		SIZE_T upper = UpperBound(l, u, strRead, s, refGenome, hashTable);

		for (SIZE_T j = lower; j <= upper; j++) {
			if (ResultExist(oneResult, hashTable->index[j] - i))
				continue;

			int s = GetKmer(refGenome, hashTable->index[j] - i, readLen, &ret);
			if (s != readLen)
				continue;
			if (oneResult->nRet >= 200) {
				printf("Array touch the Boundary!\n");
				break;
			}

			SIZE_T nDiff = bitsStrNCompare(ret, readInBits, readLen);
			if (nDiff <= matchOpt->nMaxMismatch) {
				oneResult->nMismatch[oneResult->nRet] = nDiff;
				oneResult->nStartPos[oneResult->nRet] = hashTable->index[j] - i;
				oneResult->nRet++;
			}
		}
	}
}

__host__ __device__ void MappingOneRead(const MatchOpt * matchOpt,
		const RefGenome * refGenome, const HashTable * hashTable,
		InBits readInBits, int readLen, ResultMatchedReads * oneResult) {
	oneResult->nRet = 0;
	Mapping(matchOpt, refGenome, hashTable, readInBits, readLen, oneResult);
	char strRead[MAX_READ_LEN];
	reverseCompliment(strRead, readInBits, readLen);
	InBits readInBits_rev;
	EncodeRead(strRead, &readInBits_rev, readLen);
	Mapping(matchOpt, refGenome, hashTable, readInBits_rev, readLen, oneResult);
}
