#include "hash.h"

void CountBucketSize(const CReference * refGenome, CHashTable * hashTable) {
	uint32_t hashValue;
	for (uint32_t i = 0; i < refGenome->nRefSize - HASHLEN; i++) {
		hashValue = GetHashValue(&(refGenome->refSeq[i]));
		hashTable->counter[hashValue]++;
	}

	for (uint32_t i = 1; i <= NO_OF_BUCKET; i++) {
		hashTable->counter[i] += hashTable->counter[i - 1];
	}

	hashTable->nSizeIndex = hashTable->counter[NO_OF_BUCKET];
	INFO("The size of Hash Table Index Array is", hashTable->nSizeIndex);

	for (uint32_t i = NO_OF_BUCKET - 1; i >= 1; i--) {
		hashTable->counter[i] = hashTable->counter[i - 1];
	}
	hashTable->counter[0] = 0;
}

void HashToBucket(const CReference * refGenome, CHashTable * hashTable) {
	INFO("The memory of Hash Table Index Array is", sizeof(uint32_t) * (hashTable->nSizeIndex / GB), "GB");
	MEMORY_ALLOCATE_CHECK(hashTable->index = (uint32_t * ) malloc(sizeof(uint32_t) * hashTable->nSizeIndex));
	uint32_t hashValue;
	for (uint32_t i = 0; i < refGenome->nRefSize - HASHLEN; i++) {
		hashValue = GetHashValue(&(refGenome->refSeq[i]));
		hashTable->index[hashTable->counter[hashValue]] = i;
		hashTable->counter[hashValue]++;
	}

	for (uint32_t i = NO_OF_BUCKET - 1; i >= 1; i--) {
		hashTable->counter[i] = hashTable->counter[i - 1];
	}
	hashTable->counter[0] = 0;
}

void hashTableTest(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {
#ifdef debug1
	ofstream fout("hashtableTest.txt");
	string refSeq = refGenome->refSeq;
	for (uint32_t i = 1; i <= NO_OF_BUCKET; i++) {
		for (uint32_t j = hashTable->counter[i]; j < hashTable->counter[i + 1]; j++) {
			fout << "i = " << i << " j = " << j << " " << hashTable->index[j] << " "
				<< refSeq.substr(hashTable->index[j], 50) << endl;
		}
	}
#endif
}

void WriteIndexAndRef(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {
	FILE * fout = fopen(opt.indexFile.c_str(), "wb");
	INFO("write reference and index to", opt.indexFile);
	fwrite(&(refGenome->nRefSize), sizeof(uint32_t), 1, fout);
	fwrite(refGenome->refSeq, sizeof(char), refGenome->nRefSize, fout);
	fwrite(&(hashTable->nSizeCounter), sizeof(uint32_t), 1, fout);
	fwrite(hashTable->counter, sizeof(uint32_t), hashTable->nSizeCounter + 1, fout);
	fwrite(&(hashTable->nSizeIndex), sizeof(uint32_t), 1, fout);
	fwrite(hashTable->index, sizeof(uint32_t), hashTable->nSizeIndex, fout);
	fclose(fout);
}

void ReadIndexAndRef(const Option & opt, CReference * refGenome, CHashTable * hashTable) {
	printf("\n");
	INFO("Read reference index from", opt.refFile);
	FILE * fin = fopen(opt.refFile.c_str(), "rb");
	FILE_OPEN_CHECK(fin);
	FREAD_CHECK(fread(&(refGenome->nRefSize), sizeof(uint32_t), 1, fin), 1);
	MEMORY_ALLOCATE_CHECK(refGenome->refSeq = (char * ) malloc(sizeof(char) * (refGenome->nRefSize + 1)));
	FREAD_CHECK(fread(refGenome->refSeq, sizeof(char), refGenome->nRefSize, fin), refGenome->nRefSize);
	refGenome->refSeq[refGenome->nRefSize] = 0;
	INFO("The size of the reference genome is", refGenome->nRefSize);

	FREAD_CHECK(fread(&(hashTable->nSizeCounter), sizeof(uint32_t), 1, fin), 1);
	INFO("The size of the hash table counter array is", hashTable->nSizeCounter);
	MEMORY_ALLOCATE_CHECK(hashTable->counter = (uint32_t * ) malloc(sizeof(uint32_t) * (hashTable->nSizeCounter + 1)));
	FREAD_CHECK(fread(hashTable->counter, sizeof(uint32_t), hashTable->nSizeCounter + 1, fin), hashTable->nSizeCounter + 1);
	FREAD_CHECK(fread(&(hashTable->nSizeIndex), sizeof(uint32_t), 1, fin), 1);
	INFO("The size of the hash table index array is", hashTable->nSizeIndex);
	MEMORY_ALLOCATE_CHECK(hashTable->index = (uint32_t * ) malloc(sizeof(uint32_t) * hashTable->nSizeIndex));
	FREAD_CHECK(fread(hashTable->index, sizeof(uint32_t), hashTable->nSizeIndex, fin), hashTable->nSizeIndex);
	fclose(fin);
}

const CReference * globalRefGenome;
bool SortHashCMP(const uint32_t & a, const uint32_t & b) {
	uint32_t i = a, j = b, cur = 0;
	while (i < globalRefGenome->nRefSize && j < globalRefGenome->nRefSize) {
		if (globalRefGenome->refSeq[i] < globalRefGenome->refSeq[j])
			return true;
		else if (globalRefGenome->refSeq[i] > globalRefGenome->refSeq[j])
			return false;
		i++;
		j++;
		cur++;
		if (cur > 10010)
			break; //otherwise sort will take a lot time
	}
	if (i == globalRefGenome->nRefSize)
		return true;
	return false;
}

void SortEachBucket(const CReference * refGenome, CHashTable * hashTable) {
	INFO("Sort each bucket...");
	globalRefGenome = refGenome;
	for (uint32_t i = 0; i < NO_OF_BUCKET; i++) {
		if (hashTable->counter[i] == hashTable->counter[i + 1])
			continue;
		sort(hashTable->index + hashTable->counter[i], hashTable->index + hashTable->counter[i + 1], SortHashCMP);
	}
	globalRefGenome = NULL;
}

void BuildHashTable(const Option & opt, const CReference * refGenome, CHashTable * hashTable) {
	printf("\n");
	INFO("Build hash table...");
	/* 64M buckets, 256 MB  4^13, 67108864, 0100000000000000000000000000 */
	hashTable->nSizeCounter = NO_OF_BUCKET;
	INFO("The size of Hash Table Counter Array is", sizeof(uint32_t) * (hashTable->nSizeCounter + 1) / GB, "GB");
	uint32_t sizeCounter = sizeof(uint32_t) * (hashTable->nSizeCounter + 1);
	MEMORY_ALLOCATE_CHECK(hashTable->counter = (uint32_t * ) malloc(sizeCounter));
	memset(hashTable->counter, 0x00, sizeCounter);

	TIME_INFO(CountBucketSize(refGenome, hashTable), "Count bucket size"); /* count each bucket size */
	TIME_INFO(HashToBucket(refGenome, hashTable), "Hash to bucket"); /* put each element into a bucket */
	TIME_INFO(SortEachBucket(refGenome, hashTable), "Sort each bucket");

	if (opt.bSaveIndex == 1) {
		TIME_INFO(WriteIndexAndRef(opt, refGenome, hashTable), "Write reference and index");
	}
}

uint32_t AnalyzeReference(char * strRef, uint32_t refLen) {
	/* This function removes all non-ACGTN characters */
	char strRet[MAX_LINE_LEN];
	uint32_t j = 0;
	for (uint32_t i = 0; i < refLen; i++) {
		if (strRef[i] == '>') {
			i += GetLineFromString(&strRef[i], strRet);
		} else if (isACGT(strRef[i]) || strRef[i] == 'N' || strRef[i] == 'n') {
			strRef[j++] = toupper(strRef[i]);
		}
	}
	strRef[j] = 0;

	/* change all 'N' to A,C,G,T with the same probability*/
	srand (time(NULL));int
		r = 0;
	for (uint32_t i = 0; i < j; i++) {
		if (strRef[i] == 'N') {
			r = rand() % 4;
			strRef[i] = getNT(r);
		}
	}

	return j;
}

void GenerateTestData(CReference * refGenome) {
	int len[] = {2000, 3000, 5000, 6000, 7000, 10000};
	char filename[100];
	string strVal = refGenome->refSeq;
	for(int i = 0;i < 6;i++) {
		sprintf(filename, "csci503_test_seqlen_%d.fa", len[i]);
		cout << filename << endl;
		ofstream fout(filename);
		for(int j = 0;j < 10000000;j++) {
			fout << ">read " << j << endl;
			uint32_t start = lrand48() % refGenome->nRefSize;
			if(start + len[i] > refGenome->nRefSize - 10) continue;
			fout << strVal.substr(start, len[i]) << endl;
		}
		fout.close();
	}
	cout << "finish" << endl;
	exit(EXIT_FAILURE);
}

void GetReference(CReference * refGenome, const Option & opt) {
	/* read reference genome from file */
	char * strRef;
	INFO("Read reference genome from", opt.refFile);
	uint32_t refLen = ReadWholeFile(opt.refFile, &strRef);
	INFO("The length of reference file is", refLen);
	refGenome->nRefSize = AnalyzeReference(strRef, refLen);
	MEMORY_ALLOCATE_CHECK(refGenome->refSeq = (char * ) malloc(sizeof(char) * (refGenome->nRefSize + 1)));
	memcpy(refGenome->refSeq, strRef, refGenome->nRefSize);
	refGenome->refSeq[refGenome->nRefSize] = 0;
	free(strRef);
#ifdef debug1
	GenerateTestData(refGenome);
#endif
}

void BuildIndex(const Option & opt, CReference * refGenome, CHashTable * hashTable) {
	GetReference(refGenome, opt);
	BuildHashTable(opt, refGenome, hashTable);
	hashTableTest(opt, refGenome, hashTable);
}
