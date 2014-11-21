#include "hash.h"

const CReference * globalRefGenome;

void CountBucketSize(const CReference * refGenome, CHashTable * hashTable,
		const set<SIZE_T> & setConsectiveN) {
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

	SIZE_T maxBucket = hashTable->counter[0];
	SIZE_T k = 0;
	for (SIZE_T i = 1; i <= NO_OF_BUCKET; i++) {
		if (hashTable->counter[i] > maxBucket) {
			maxBucket = hashTable->counter[i];
			k = i;
		}
		hashTable->counter[i] += hashTable->counter[i - 1];
	}
	printf("k = %d\n", k);
	printf("maxBucket = %d\n", maxBucket);
	for (SIZE_T i = 0; i < 40; i++) {
		cout << "i = " << i << " " << hashTable->counter[i] << endl;
	}

	hashTable->nSizeIndex = hashTable->counter[NO_OF_BUCKET];
	cout << "hashTable->nSizeIndex = " << hashTable->nSizeIndex << endl;

	for (SIZE_T i = NO_OF_BUCKET - 1; i >= 1; i--) {
		hashTable->counter[i] = hashTable->counter[i - 1];
	}
	hashTable->counter[0] = 0;
}

void HashToBucket(const CReference * refGenome, CHashTable * hashTable,
		const set<SIZE_T> & setConsectiveN) {
	LOG_INFO;
	MEMORY_ALLOCATE_CHECK(
			hashTable->index = (SIZE_T * ) malloc(sizeof(SIZE_T) * hashTable->nSizeIndex));
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
		SIZE_T start = hashTable->counter[i];
		if (hashTable->counter[i + 1] <= start)
			continue;
		for (SIZE_T j = start; j < hashTable->counter[i + 1]; j++) {
			SIZE_T s = GetKmer(refGenome, hashTable->index[j], 64, &r);
			InBits value;
			int ss = GetF2SeedForBits(r, s, &value);
			char strRead[MAX_READ_LEN];
			DecodeReadReverse(strRead, ss, &value);
			fprintf(fsee, " j = %09d, i=%d %09d ", j, i, hashTable->index[j]);
			fprintf(fsee, "%s ", strRead);
			DecodeRead(strRead, s, &r);
			fprintf(fsee, "%s\n", strRead);
		}
		fprintf(fsee, "\n");
	}
	fclose(fsee);
	LOG_INFO;
}

int SortHashCMP(const void * a, const void * b) {

	SIZE_T nRefStart1 = *(SIZE_T *) a;
	SIZE_T nRefStart2 = *(SIZE_T *) b;
	InBits r1; //+22
	SIZE_T indexInWords1 = nRefStart1 / wordSize;
	SIZE_T bitsShift1 = nRefStart1 % wordSize;
	r1.ub = globalRefGenome->refInBits[indexInWords1].ub >> bitsShift1;
	r1.lb = globalRefGenome->refInBits[indexInWords1].lb >> bitsShift1;

	InBits r2;
	SIZE_T indexInWords2 = nRefStart2 / wordSize;
	SIZE_T bitsShift2 = nRefStart2 % wordSize;
	r2.ub = globalRefGenome->refInBits[indexInWords2].ub >> bitsShift2;
	r2.lb = globalRefGenome->refInBits[indexInWords2].lb >> bitsShift2;

	///////////////////////////////////////////////////////////////////
	InBits r;
	SIZE_T s = GetKmer(globalRefGenome, nRefStart1, 64, &r);
	InBits value;

	int ss = GetF2SeedForBits(r, s, &value);
	char strRead[MAX_READ_LEN];
	DecodeReadReverse(strRead, ss, &value);
	//printf("strRead = %s ", strRead);
	DecodeRead(strRead, s, &r);
	//printf("strRead = %s\n", strRead);

	s = GetKmer(globalRefGenome, nRefStart2, 64, &r);
	ss = GetF2SeedForBits(r, s, &value);
	DecodeReadReverse(strRead, ss, &value);
	//printf("strRead = %s ", strRead);
	DecodeRead(strRead, s, &r);
	//printf("strRead = %s\n", strRead);
	/////////////////////////////////////////////////////////////////////

	int cmpbits = 0;

	WORD_SIZE code = 0x17; //0010111

	int bits1 = wordSize - bitsShift1;
	int bits2 = wordSize - bitsShift2;
	//printf("(s = %d %d)", nRefStart1, nRefStart2);
	while (cmpbits < 1000) {
		if (cmpbits % 7 == 0)
			code = 0x17;
		if (code & 0x01) {
			SIZE_T c1 = (r1.ub & 0x01) << 1 | (r1.lb & 0x01);
			SIZE_T c2 = (r2.ub & 0x01) << 1 | (r2.lb & 0x01);
			//cout << "(" << c1 << "," << c2 << ")";
			if (c1 > c2) {
				//printf("x\n");
				return 1;
			} else if (c1 < c2) {
				//printf("g\n");
				return -1;
			}
		}
		code >>= 1;
		r1.ub >>= 1;
		r1.lb >>= 1;
		r2.ub >>= 1;
		r2.lb >>= 1;
		cmpbits++;
		bits1--;
		bits2--;

		if (nRefStart1 + cmpbits >= globalRefGenome->nRefSize)
			return -1;
		if (nRefStart2 + cmpbits >= globalRefGenome->nRefSize)
			return 1;

		if (bits1 == 0) {
			indexInWords1++;
			r1 = globalRefGenome->refInBits[indexInWords1];
			bits1 = wordSize;
		}

		if (bits2 == 0) {
			indexInWords2++;
			r2 = globalRefGenome->refInBits[indexInWords2];
			bits2 = wordSize;
		}
	}
	//printf("\n");
	return 0;
}

void SortEachBucket(const CReference * refGenome, CHashTable * hashTable) {
	LOG_INFO;
	globalRefGenome = refGenome;
	for (SIZE_T i = 0; i < NO_OF_BUCKET; i++) {
		SIZE_T start = hashTable->counter[i];
		SIZE_T n = hashTable->counter[i + 1] - start;
		if (n > 1) {
			qsort(&(hashTable->index[start]), n, sizeof(SIZE_T), SortHashCMP);
			//cout << endl;
		}
	}
	globalRefGenome = NULL;
}

void WriteIndexAndRef(const CReference * refGenome,
		const CHashTable * hashTable, const Option & opt) {
	LOG_INFO;
	FILE * fout = fopen(opt.indexFile, "wb");
	cout << "write reference and index to " << opt.indexFile << endl;
	fwrite(&(refGenome->nRefSize), sizeof(SIZE_T), 1, fout);
	fwrite(&(refGenome->nRefSizeInWordSize), sizeof(SIZE_T), 1, fout);
	fwrite(refGenome->refInBits, sizeof(InBits), refGenome->nRefSizeInWordSize,
			fout);
	fwrite(&(hashTable->nSizeCounter), sizeof(SIZE_T), 1, fout);
	fwrite(hashTable->counter, sizeof(SIZE_T), hashTable->nSizeCounter, fout);
	fwrite(&(hashTable->nSizeIndex), sizeof(SIZE_T), 1, fout);
	fwrite(hashTable->index, sizeof(SIZE_T), hashTable->nSizeIndex, fout);
	fclose(fout);
}

void ReadIndexAndRef(CReference * refGenome, CHashTable * hashTable,
		const Option & opt) {
	LOG_INFO;
	FILE * fin = fopen(opt.indexFile, "rb");
	FILE_OPEN_CHECK(fin);
	FREAD_CHECK(fread(&(refGenome->nRefSize), sizeof(SIZE_T), 1, fin), 1);
	FREAD_CHECK(fread(&(refGenome->nRefSizeInWordSize), sizeof(SIZE_T), 1, fin),
			1);
	MEMORY_ALLOCATE_CHECK(
			refGenome->refInBits = (InBits * ) malloc(sizeof(InBits) * refGenome->nRefSizeInWordSize));
	FREAD_CHECK(
			fread(refGenome->refInBits, sizeof(InBits), refGenome->nRefSizeInWordSize, fin),
			refGenome->nRefSizeInWordSize);

	FREAD_CHECK(fread(&(hashTable->nSizeCounter), sizeof(SIZE_T), 1, fin), 1);
	MEMORY_ALLOCATE_CHECK(
			hashTable->counter = (SIZE_T * ) malloc(sizeof(SIZE_T) * hashTable->nSizeCounter));
	FREAD_CHECK(
			fread(hashTable->counter, sizeof(SIZE_T), hashTable->nSizeCounter, fin),
			hashTable->nSizeCounter);

	FREAD_CHECK(fread(&(hashTable->nSizeIndex), sizeof(SIZE_T), 1, fin), 1);
	MEMORY_ALLOCATE_CHECK(
			hashTable->index = (SIZE_T * ) malloc(sizeof(SIZE_T) * hashTable->nSizeIndex));
	FREAD_CHECK(
			fread(hashTable->index, sizeof(SIZE_T), hashTable->nSizeIndex, fin),
			hashTable->nSizeIndex);
	fclose(fin);
}

void MakeHashTable(const CReference * refGenome, CHashTable * hashTable,
		Option & opt) {
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
	TIME_INFO(CountBucketSize(refGenome, hashTable, opt.setConsectiveN),
			"count bucket size");
	cout << "hf6" << endl;
	/* put each element into a bucket */
	TIME_INFO(HashToBucket(refGenome, hashTable, opt.setConsectiveN),
			"hash to bucket");
	cout << "sizecon = " << opt.setConsectiveN.size() << endl;
	opt.setConsectiveN.clear();
	cout << "hf5" << endl;
	/* sort each bucket */
	TIME_INFO(SortEachBucket(refGenome, hashTable), "sort each bucket");
	testHashTalbe(refGenome, hashTable);
	cout << "hf7" << endl;
	if (opt.bSaveIndex) {
		TIME_INFO(WriteIndexAndRef(refGenome, hashTable, opt),
				"write reference and index");
	}
	cout << "hf8" << endl;
}
