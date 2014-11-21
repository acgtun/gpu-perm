#include "hash.h"

const InBits * globalRefGenome;
SIZE_T GLOBALnRefSize;

void CountBucketSize(const InBits * refGenome, SIZE_T nRefSize,
		SIZE_T * hashCounter) {
	LOG_INFO;
	InBits r;
	for (SIZE_T i = 0; i < nRefSize; i++) {
		if (GetKmer(refGenome, nRefSize, i, 22, &r) < 22)
			break;
		SIZE_T hashValue = GetHashValue(r);
		hashCounter[hashValue]++;
	}
	for (SIZE_T i = 1; i <= NO_OF_BUCKET; i++) {
		hashCounter[i] += hashCounter[i - 1];
	}

	for (SIZE_T i = NO_OF_BUCKET - 1; i >= 1; i--) {
		hashCounter[i] = hashCounter[i - 1];
	}
	hashCounter[0] = 0;
}

void HashToBucket(const InBits * refGenome, SIZE_T nRefSize,
		SIZE_T * hashCounter, SIZE_T ** hashIndex) {
	LOG_INFO;
	*hashIndex = (SIZE_T *) malloc(sizeof(SIZE_T) * (nRefSize + 1));
	if (*hashIndex == NULL)
		MEMORY_ALLOCATE_ERROR
//	FILE * fsee = fopen("seeeee.txt", "wb");
//	for (int i = 0; i < NO_OF_BUCKET + 1; i++) {
//		fprintf(fsee, "hashCounter[i%d] = %d\n", i, hashCounter[i]);
//	}
//	fclose(fsee);

	InBits r;
	//cout << nRefSize << endl;
	for (SIZE_T i = 0; i < nRefSize; i++) {
		//cout << "i = " << i << " " << nRefSize << endl;
		if (GetKmer(refGenome, nRefSize, i, 22, &r) < 22)
			break;
		//cout << "good" << endl;
		SIZE_T hashValue = GetHashValue(r);
		//cout << hashValue << endl;
		//cout << "good too" << endl;
		//cout << hashCounter[hashValue] << endl;
		(*hashIndex)[hashCounter[hashValue]++] = i; /* make sure -- is not less than zero */
		//cout << "good too too" << endl;
	}
	LOG_INFO
	for (SIZE_T i = NO_OF_BUCKET - 1; i >= 1; i--) {
		hashCounter[i] = hashCounter[i - 1];
	}
	hashCounter[0] = 0;
}

void testHashTalbe(const InBits * refGenome, SIZE_T nRefSize,
		SIZE_T * hashCounter, SIZE_T * hashIndex) {
	LOG_INFO;
	InBits r;
	FILE * fsee = fopen("test_out_sort.txt", "wb");
	LOG_INFO
	for (SIZE_T i = 0; i < NO_OF_BUCKET; i++) {
		int start = hashCounter[i];
		for (SIZE_T j = start; j < hashCounter[i + 1]; j++) {
			SIZE_T s = GetKmer(refGenome, nRefSize, hashIndex[j], 64, &r);
			InBits value;
			int ss = GetF2SeedForBits(r, s, &value);
			char strRead[MAX_READ_LEN];
			DecodeReadReverse(strRead, ss, &value);
			fprintf(fsee, " j = %08d, i=%d %08d ", j, i, hashIndex[j]);
			fprintf(fsee, "strRead = %s ", strRead);
			DecodeRead(strRead, s, &r);
			fprintf(fsee, "strRead = %s\n", strRead);
		}
		fprintf(fsee, "\n");
	}
	fclose(fsee);
	LOG_INFO;
}

int Compare(const void * a, const void * b) {

	SIZE_T index1 = *(SIZE_T *) a;
	SIZE_T index2 = *(SIZE_T *) b;

	InBits r1, r2;
	SIZE_T s1 = GetKmer(globalRefGenome, GLOBALnRefSize, index1, wordSize, &r1);
	SIZE_T s2 = GetKmer(globalRefGenome, GLOBALnRefSize, index2, wordSize, &r2);
	InBits v1, v2;
	int ss1 = GetF2SeedForBits(r1, s1, &v1);
	int ss2 = GetF2SeedForBits(r2, s2, &v2);

	char strRead1[MAX_READ_LEN], strRead2[MAX_READ_LEN];
	DecodeReadReverse(strRead1, ss1, &v1);
	DecodeReadReverse(strRead2, ss2, &v2);

	return strcmp(strRead1, strRead2);
}

void SortEachBucket(const InBits * refGenome, SIZE_T nRefSize,
		SIZE_T * hashCounter, SIZE_T * hashIndex) {
	LOG_INFO;
	globalRefGenome = refGenome;
	GLOBALnRefSize = nRefSize;
	for (SIZE_T i = 0; i < NO_OF_BUCKET; i++) {
		int start = hashCounter[i];
		int n = hashCounter[i + 1] - start;
		if (n > 1) {
			qsort(&(hashIndex[start]), n, sizeof(SIZE_T), Compare);
		}
	}
	globalRefGenome = NULL;
	//testHashTalbe(refGenome, hashTable);
}
void MakeHashTable(const InBits * refGenome, SIZE_T nRefSize,
		SIZE_T ** hashCounter, SIZE_T ** hashIndex) {
	LOG_INFO
	/* 64M buckets, 256 MB  4^13, 67108864, 0100000000000000000000000000 */
	//hashTable->NO_OF_BUCKET = 0x4000000;
	*hashCounter = (SIZE_T *) malloc(sizeof(SIZE_T) * (NO_OF_BUCKET + 1));
	if (*hashCounter == NULL)
		MEMORY_ALLOCATE_ERROR
	cout << sizeof(*hashCounter) << endl;
	cout << NO_OF_BUCKET << endl;
	LOG_INFO
	memset(*hashCounter, 0x00, sizeof(SIZE_T) * (NO_OF_BUCKET + 1));
	LOG_INFO
	/* count each bucket size */
	TIME_INFO(CountBucketSize(refGenome, nRefSize, *hashCounter),
			"count bucket size");
	/* put each element into a bucket */
	TIME_INFO(HashToBucket(refGenome, nRefSize, *hashCounter,hashIndex),
			"hash to bucket");
	/* sort each bucket */
	TIME_INFO(SortEachBucket(refGenome, nRefSize, *hashCounter, *hashIndex),
			"sort each bucket");
}
