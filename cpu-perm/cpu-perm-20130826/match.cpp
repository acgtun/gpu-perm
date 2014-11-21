#include "match.h"

void reverseCompliment(char * strRead_rev, const char * strRead, const SIZE_T & len) {
	for (SIZE_T i = 0; i < len; i++) {
		strRead_rev[i] = complimentBase(strRead[len - i - 1]);
	}
	strRead_rev[len] = 0;
}

void Reverse_Kernel(const CRead * reads, CRead * reads_rev, const SIZE_T & nReadsNum) {
	for (SIZE_T i = 0; i < nReadsNum; i++) {
		reverseCompliment(reads_rev[i].readInStr, reads[i].readInStr, reads[i].readLen);
		reads_rev[i].readLen = reads[i].readLen;
	}
}

int GetNTFromGenome(const SIZE_T & index, const CReference * refGenome) {
	if (index >= refGenome->nRefSize)
		return 0;

	InBits ref;
	SIZE_T indexInWords = index / wordSize;
	SIZE_T bitsShift = index % wordSize;
	ref.ub = refGenome->refInBits[indexInWords].ub >> bitsShift;
	ref.lb = refGenome->refInBits[indexInWords].lb >> bitsShift;

	return getNT((ref.ub & 0x01) << 1 | (ref.lb & 0x01));
}

SIZE_T LowerBound(SIZE_T low, SIZE_T high, const char & chr, const SIZE_T & cmpbits, const CReference * refGenome, const CHashTable * hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high) / 2;
		if (GetNTFromGenome(hashTable->index[mid] + cmpbits, refGenome) >= chr)
			high = mid;
		else
			low = mid + 1;
	}
	return low;
}

SIZE_T UpperBound(SIZE_T low, SIZE_T high, const char & chr, const SIZE_T & cmpbits, const CReference * refGenome, const CHashTable * hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high + 1) / 2;
		if (GetNTFromGenome(hashTable->index[mid] + cmpbits, refGenome) <= chr)
			low = mid;
		else
			high = mid - 1;
	}
	return low;
}

void NoRegion(CRegion * ret) {
	ret->lower = 1;
	ret->upper = 0;
}

void Match(const CReference * refGenome, const CHashTable * hashTable, const CRead * read, const SIZE_T & nStartPos, CRegion * ret) {

	SIZE_T hashValue = GetHashValue(read->readInStr, nStartPos);
	SIZE_T l = hashTable->counter[hashValue];
	if (hashTable->counter[hashValue + 1] == 0)
		return NoRegion(ret);
	SIZE_T u = hashTable->counter[hashValue + 1] - 1;
	if (l > u)
		return NoRegion(ret);

	WORD_SIZE code = 0x17; //0010111
	code >>= 1;
	for (SIZE_T i = nStartPos + 22, cmpbits = 22; i < read->readLen; i++, cmpbits++) {
		if (cmpbits % 7 == 0)
			code = 0x17;
		if (!(code & 0x01)) {
			code >>= 1;
			continue;
		}
		l = LowerBound(l, u, read->readInStr[i], cmpbits, refGenome, hashTable);
		u = UpperBound(l, u, read->readInStr[i], cmpbits, refGenome, hashTable);
		code >>= 1;
		if (l > u)
			return NoRegion(ret);
	}

	ret->lower = l;
	ret->upper = u;
}
