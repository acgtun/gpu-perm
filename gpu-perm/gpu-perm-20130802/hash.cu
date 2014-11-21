#include "hash.h"

const CReference * globalRefGenome;

void CountBucketSize(const CReference * refGenome, CHashTable * hashTable) {
	LOG_INFO;
	InBits r;
	for (SIZE_T i = 0; i < refGenome->nRefSize; i++) {
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

void HashToBucket(const CReference * refGenome, CHashTable * hashTable) {
	LOG_INFO;
	MEMORY_ALLOCATE_CHECK(hashTable->index = (SIZE_T * ) malloc(sizeof(SIZE_T) * (refGenome->nRefSize + 1)));
	InBits r;
	for (SIZE_T i = 0; i < refGenome->nRefSize; i++) {
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

void testHashTalbe(const CReference * refGenome, CHashTable * hashTable) {
	LOG_INFO;
	InBits r;
	FILE * fsee = fopen("test_out_sort.txt", "wb");
	LOG_INFO;
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
	LOG_INFO;
}

__device__ void QuickSort(const CReference * refGenome, CHashTable * hashTable, SIZE_T low, SIZE_T high) {
	/* The QuickSort function was modified from this website:http://alienryderflex.com/quicksort */
	int beg[QuickSort_MAX_LEVELS], end[QuickSort_MAX_LEVELS], i = 0, L, R, swap;

	beg[0] = low;
	end[0] = high + 1;

	while (i >= 0) {
		L = beg[i];
		R = end[i] - 1;
		if (L < R) {
			//piv=arr[L];
			SIZE_T tmpL = L;
			InBits r, v;
			SIZE_T s = GetKmer(refGenome, L, wordSize, &r);
			int ss = GetF2SeedForBits(r, s, &v);
			char piv[MAX_READ_LEN];
			DecodeReadReverse(piv, ss, &v);
			while (L < R) {
				while (CMP(piv, ss, R, refGenome) <= 0 && L < R)
					R--;
				if (L < R) {
					hashTable->index[L] = hashTable->index[R];
					L++;
				}
				while (CMP(piv, ss, L, refGenome) >= 0 && L < R)
					L++;
				if (L < R) {
					hashTable->index[R] = hashTable->index[L];
					R--;
				}
			}
			hashTable->index[L] = tmpL;
			beg[i + 1] = L + 1;
			end[i + 1] = end[i];
			end[i++] = L;
			if (end[i] - beg[i] > end[i - 1] - beg[i - 1]) {
				swap = beg[i];
				beg[i] = beg[i - 1];
				beg[i - 1] = swap;
				swap = end[i];
				end[i] = end[i - 1];
				end[i - 1] = swap;
			}
		} else {
			i--;
		}
	}
}

__global__ void SortEachBucket_kernel(const CReference refGenome, CHashTable hashTable) {
	SIZE_T stride = blockDim.x * blockIdx.x;
	SIZE_T i = threadIdx.x + stride;
	while (i < NO_OF_BUCKET) {
		int start = hashTable.counter[i];
		int n = hashTable.counter[i + 1] - start;
		if (n > 1) {
			QuickSort(&refGenome, &hashTable, start, hashTable.counter[i + 1] - 1);
		}
		i += stride;
	}
}

void MakeHashTable(const CReference * refGenome, CReference * d_refGenome, CHashTable * d_hashTable) {
	LOG_INFO;
	/* 64M buckets, 256 MB  4^13, 67108864, 0100000000000000000000000000 */
	//hashTable->NO_OF_BUCKET = 0x4000000;
	CHashTable hashTable;
	SIZE_T sizeCounter = sizeof(SIZE_T) * (NO_OF_BUCKET + 1);
	MEMORY_ALLOCATE_CHECK(hashTable.counter = (SIZE_T * ) malloc(sizeCounter));
	memset(hashTable.counter, 0x00, sizeCounter);

	/* count each bucket size */
	TIME_INFO(CountBucketSize(refGenome, &hashTable), "count bucket size");
	/* put each element into a bucket */
	TIME_INFO(HashToBucket(refGenome, &hashTable), "hash to bucket");
	testHashTalbe(refGenome, &hashTable);
	/* sort each bucket */
	HANDLE_ERROR(cudaMalloc((void** )&(d_hashTable->counter), sizeCounter));
	HANDLE_ERROR(cudaMemcpy(d_hashTable->counter, hashTable.counter, sizeCounter, cudaMemcpyHostToDevice));
	SIZE_T sizeIndex = sizeof(SIZE_T) * (refGenome->nRefSize + 1);
	HANDLE_ERROR(cudaMalloc((void** )&(d_hashTable->index), sizeIndex));
	HANDLE_ERROR(cudaMemcpy(d_hashTable->index, hashTable.index, sizeIndex, cudaMemcpyHostToDevice));

	SIZE_T sizeRef = sizeof(InBits) * (refGenome->nRefSizeInWordSize + 1);
	HANDLE_ERROR(cudaMalloc((void** )&(d_refGenome->refInBits), sizeRef));
	HANDLE_ERROR(cudaMemcpy(d_refGenome->refInBits, refGenome->refInBits, sizeRef, cudaMemcpyHostToDevice));
	d_refGenome->nRefSize = refGenome->nRefSize;
	d_refGenome->nRefSizeInWordSize = refGenome->nRefSizeInWordSize;

	free(refGenome->refInBits);
	free(hashTable.counter);
	free(hashTable.index);
	SortEachBucket_kernel<<<BLOCKS, THREADS>>>(*d_refGenome, *d_hashTable);
}
