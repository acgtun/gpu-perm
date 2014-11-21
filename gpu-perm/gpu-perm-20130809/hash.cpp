#include "hash.h"

const CReference * globalRefGenome;

void CountBucketSize(const CReference * refGenome, CHashTable * hashTable, const set<SIZE_T> & setConsectiveN) {
	LOG_INFO;
	InBits r;
	for (SIZE_T i = 0; i < refGenome->nRefSize; i++) {
		if (i % NUMBER_OF_SPACE != 0)
			continue;
		if (setConsectiveN.find(i) != setConsectiveN.end())
			continue;
		if (GetKmer(refGenome, i, 22, &r) < 22)
			break;
		SIZE_T hashValue = GetHashValue(r);
		hashTable->counter[hashValue]++;
	}

	for (SIZE_T i = 1; i <= NO_OF_BUCKET; i++) {
		hashTable->counter[i] += hashTable->counter[i - 1];
	}
	hashTable->nSizeIndex = hashTable->counter[NO_OF_BUCKET];
	cout << "hashTable->nSizeIndex = " << hashTable->nSizeIndex << endl;
	for (SIZE_T i = NO_OF_BUCKET - 1; i >= 1; i--) {
		hashTable->counter[i] = hashTable->counter[i - 1];
	}
	hashTable->counter[0] = 0;
}

void HashToBucket(const CReference * refGenome, CHashTable * hashTable, const set<SIZE_T> & setConsectiveN) {
	LOG_INFO;
	MEMORY_ALLOCATE_CHECK(hashTable->index = (SIZE_T * ) malloc(sizeof(SIZE_T) * hashTable->nSizeIndex));
	InBits r;
	for (SIZE_T i = 0; i < refGenome->nRefSize; i++) {
		if (i % NUMBER_OF_SPACE != 0)
			continue;
		if (setConsectiveN.find(i) != setConsectiveN.end())
			continue;
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

void testHashTalbe(const CReference * refGenome, const CHashTable * hashTable) {
	LOG_INFO;
	InBits r;
	FILE * fsee = fopen("test_out_sort.txt", "wb");
	LOG_INFO;
	for (SIZE_T i = 0; i < NO_OF_BUCKET; i++) {
		int start = hashTable->counter[i];
		if (hashTable->counter[i + 1] <= start)
			continue;
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
	LOG_INFO;
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

void SortEachBucket(const CReference * refGenome, CHashTable * hashTable) {
	LOG_INFO;
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

void WriteIndexAndRef(const CReference * refGenome, const CHashTable * hashTable, const Option & opt) {
	LOG_INFO;
	FILE * fout = fopen(opt.indexFile, "wb");
	cout << "write reference and index to " << opt.indexFile << endl;
	fwrite(&(refGenome->nRefSize), sizeof(SIZE_T), 1, fout);
	fwrite(&(refGenome->nRefSizeInWordSize), sizeof(SIZE_T), 1, fout);
	fwrite(refGenome->refInBits, sizeof(InBits), refGenome->nRefSizeInWordSize, fout);
	fwrite(&(hashTable->nSizeCounter), sizeof(SIZE_T), 1, fout);
	fwrite(hashTable->counter, sizeof(SIZE_T), hashTable->nSizeCounter, fout);
	fwrite(&(hashTable->nSizeIndex), sizeof(SIZE_T), 1, fout);
	fwrite(hashTable->index, sizeof(SIZE_T), hashTable->nSizeIndex, fout);
	fclose(fout);

	LOG_INFO;
	FILE * foutnotb = fopen("texxxst.index", "w");
	fwrite(&(refGenome->nRefSize), sizeof(SIZE_T), 1, foutnotb);
	fwrite(&(refGenome->nRefSizeInWordSize), sizeof(SIZE_T), 1, foutnotb);
	fwrite(refGenome->refInBits, sizeof(InBits), refGenome->nRefSizeInWordSize, foutnotb);
	fwrite(&(hashTable->nSizeCounter), sizeof(SIZE_T), 1, foutnotb);
	fwrite(hashTable->counter, sizeof(SIZE_T), hashTable->nSizeCounter, foutnotb);
	fwrite(&(hashTable->nSizeIndex), sizeof(SIZE_T), 1, foutnotb);
	fwrite(hashTable->index, sizeof(SIZE_T), hashTable->nSizeIndex, foutnotb);
	fclose(foutnotb);
}

void ReadIndexAndRef(CReference * refGenome, CHashTable * hashTable, const Option & opt) {
	LOG_INFO;
	FILE * fin = fopen(opt.indexFile, "rb");
	FILE_OPEN_CHECK(fin);
	SIZE_T size;
	FREAD_CHECK(fread(&(refGenome->nRefSize), sizeof(SIZE_T), 1, fin), 1);
	FREAD_CHECK(fread(&(refGenome->nRefSizeInWordSize), sizeof(SIZE_T), 1, fin), 1);
	MEMORY_ALLOCATE_CHECK(refGenome->refInBits = (InBits * ) malloc(sizeof(InBits) * refGenome->nRefSizeInWordSize));
	FREAD_CHECK(fread(refGenome->refInBits, sizeof(InBits), refGenome->nRefSizeInWordSize, fin), refGenome->nRefSizeInWordSize);

	FREAD_CHECK(fread(&(hashTable->nSizeCounter), sizeof(SIZE_T), 1, fin), 1);
	MEMORY_ALLOCATE_CHECK(hashTable->counter = (SIZE_T * ) malloc(sizeof(SIZE_T) * hashTable->nSizeCounter));
	FREAD_CHECK(fread(hashTable->counter, sizeof(SIZE_T), hashTable->nSizeCounter, fin), hashTable->nSizeCounter);

	FREAD_CHECK(fread(&(hashTable->nSizeIndex), sizeof(SIZE_T), 1, fin), 1);
	MEMORY_ALLOCATE_CHECK(hashTable->index = (SIZE_T * ) malloc(sizeof(SIZE_T) * hashTable->nSizeIndex));
	FREAD_CHECK(fread(hashTable->index, sizeof(SIZE_T), hashTable->nSizeIndex, fin), hashTable->nSizeIndex);
	fclose(fin);
}

void MakeHashTable(const CReference * refGenome, CHashTable * hashTable, Option & opt) {
	LOG_INFO;
	/* 64M buckets, 256 MB  4^13, 67108864, 0100000000000000000000000000 */
	//hashTable->NO_OF_BUCKET = 0x4000000;
	cout << "hf4" << endl;
	hashTable->nSizeCounter = NO_OF_BUCKET;
	SIZE_T sizeCounter = sizeof(SIZE_T) * (hashTable->nSizeCounter + 1);
	MEMORY_ALLOCATE_CHECK(hashTable->counter = (SIZE_T * ) malloc(sizeCounter));
	memset(hashTable->counter, 0x00, sizeCounter);
	cout << "hf3" << endl;
	/* count each bucket size */
	TIME_INFO(CountBucketSize(refGenome, hashTable, opt.setConsectiveN), "count bucket size");
	cout << "hf6" << endl;
	/* put each element into a bucket */
	TIME_INFO(HashToBucket(refGenome, hashTable, opt.setConsectiveN), "hash to bucket");
	cout << "sizecon = " << opt.setConsectiveN.size() << endl;
	opt.setConsectiveN.clear();
	cout << "hf5" << endl;
	/* sort each bucket */
	TIME_INFO(SortEachBucket(refGenome, hashTable), "sort each bucket");
	cout << "hf7" << endl;
	TIME_INFO(WriteIndexAndRef(refGenome, hashTable, opt), "write reference and index");
	cout << "hf8" << endl;
}
