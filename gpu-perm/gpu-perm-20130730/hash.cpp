#include "hash.h"

const RefGenome * globalRefGenome;

void MakeHashTable(const RefGenome * refGenome, HashTable * hashTable) {
	LOG_INFO
	;
	/* 64M buckets, 256 MB  4^13, 67108864, 0100000000000000000000000000 */
	//hashTable->NO_OF_BUCKET = 0x4000000;
	hashTable->counter = (SIZE_T *) malloc(sizeof(SIZE_T) * (NO_OF_BUCKET + 1));
	memset(hashTable->counter, 0x00, sizeof(SIZE_T) * (NO_OF_BUCKET + 1));

	/* count each bucket size */
	TIME_INFO(CountBucketSize(refGenome, hashTable), "count bucket size");
	/* put each element into a bucket */
	TIME_INFO(HashToBucket(refGenome, hashTable), "hash to bucket");
	/* sort each bucket */
	TIME_INFO(SortEachBucket(refGenome, hashTable), "sort each bucket");
}

void CountBucketSize(const RefGenome * refGenome, HashTable * hashTable) {
	LOG_INFO
	;
	InBits r;
	for (SIZE_T i = 0; i < *(refGenome->nRefSize); i++) {
		if (GetKmer(refGenome, i, 22, &r) < 22)
			break;
		SIZE_T hashValue = GetHashValue(r);
		hashTable->counter[hashValue]++;
	}
	for (SIZE_T i = 1; i <= NO_OF_BUCKET; i++) {
		hashTable->counter[i] += hashTable->counter[i - 1];
	}

	for (SIZE_T i = NO_OF_BUCKET - 1; i >= 1; i--) {
		hashTable->counter[i] = hashTable->counter[i - 1];
	}
	hashTable->counter[0] = 0;
}

void HashToBucket(const RefGenome * refGenome, HashTable * hashTable) {
	LOG_INFO
	;
	hashTable->index = (SIZE_T *) malloc(
			sizeof(SIZE_T) * (*(refGenome->nRefSize) + 1));
	InBits r;
	for (SIZE_T i = 0; i < *(refGenome->nRefSize); i++) {
		if (GetKmer(refGenome, i, 22, &r) < 22)
			break;
		SIZE_T hashValue = GetHashValue(r);
		hashTable->index[hashTable->counter[hashValue]++] = i; /* make sure -- is not less than zero */
	}

	for (SIZE_T i = NO_OF_BUCKET - 1; i >= 1; i--) {
		hashTable->counter[i] = hashTable->counter[i - 1];
	}
	hashTable->counter[0] = 0;
}

void testHashTalbe(const RefGenome * refGenome, const HashTable * hashTable) {
	LOG_INFO
	;
	InBits r;
	FILE * fsee = fopen("test_out_sort.txt", "wb");
	LOG_INFO
	for (SIZE_T i = 0; i < NO_OF_BUCKET; i++) {
		int start = hashTable->counter[i];

		for (SIZE_T j = start; j < hashTable->counter[i + 1]; j++) {
			SIZE_T s = GetKmer(refGenome, hashTable->index[j], 64, &r);
			InBits value;
			int ss = GetF2SeedForBits(r, s, &value);
			char strRead[MAX_READ_LEN];
			DecodeReadReverse(strRead, ss, &value);
			fprintf(fsee, " j = %08d, i=%d %08d ", j, i, hashTable->index[j]);
			fprintf(fsee, "strRead = %s ", strRead);
			DecodeRead(strRead, s, &r);
			fprintf(fsee, "strRead = %s\n", strRead);
		}
		fprintf(fsee, "\n");
	}
	fclose(fsee);
	LOG_INFO
	;
}

void SortEachBucket(const RefGenome * refGenome, HashTable * hashTable) {
	LOG_INFO
	;
	globalRefGenome = refGenome;
	for (SIZE_T i = 0; i < NO_OF_BUCKET; i++) {
		int start = hashTable->counter[i];
		int n = hashTable->counter[i + 1] - start;
		if (n > 1) {
			qsort(&(hashTable->index[start]), n, sizeof(SIZE_T), Compare);
		}
	}
	globalRefGenome = NULL;
	//testHashTalbe(refGenome, hashTable);
}

int Compare(const void * a, const void * b) {

	SIZE_T index1 = *(SIZE_T *) a;
	SIZE_T index2 = *(SIZE_T *) b;

	InBits r1, r2;
	SIZE_T s1 = GetKmer(globalRefGenome, index1, wordSize, &r1);
	SIZE_T s2 = GetKmer(globalRefGenome, index2, wordSize, &r2);
	InBits v1, v2;
	int ss1 = GetF2SeedForBits(r1, s1, &v1);
	int ss2 = GetF2SeedForBits(r2, s2, &v2);

	char strRead1[MAX_READ_LEN], strRead2[MAX_READ_LEN];
	DecodeReadReverse(strRead1, ss1, &v1);
	DecodeReadReverse(strRead2, ss2, &v2);

	return strcmp(strRead1, strRead2);
}
