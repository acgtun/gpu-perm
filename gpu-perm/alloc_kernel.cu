#include "alloc_kernel.h"

__global__ void Match_Kernel(const CReference * refGenome, const CHashTable * hashTable, const CRead * reads,
		const CRead * reads_rev, CRegion * d_region, const SIZE_T NUM) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	SIZE_T stride = blockDim.x * gridDim.x;

	while (i < NUM) {
		if (i % TOTAL_SHIFT_rev < TOTAL_SHIFT) {
			Match(refGenome, hashTable, &reads[i / TOTAL_SHIFT_rev], i % TOTAL_SHIFT, &(d_region[i]));
		} else {
			Match(refGenome, hashTable, &reads_rev[i / TOTAL_SHIFT_rev], i % TOTAL_SHIFT, &(d_region[i]));
		}
		i += stride;
	}
}

//
//SIZE_T NumOfDiff(const MapOpt & mapOpt, const CRead * read, SIZE_T index, const CReference * refGenome) {
//
//	InBits ref; //+22
//	SIZE_T indexInWords = index / wordSize;
//	SIZE_T bitsShift = index % wordSize;
//	ref.ub = refGenome->refInBits[indexInWords].ub >> bitsShift;
//	ref.lb = refGenome->refInBits[indexInWords].lb >> bitsShift;
//
//	SIZE_T cmpbits = 0;
//	int bits2 = wordSize - bitsShift;
//
//	SIZE_T nDiff = 0;
//	while (cmpbits < read->readLen) {
//		char c1 = read->readInStr[cmpbits];
//		char c2 = getNT((ref.ub & 0x01) << 1 | (ref.lb & 0x01));
//
//		if (c1 != c2)
//			nDiff++;
//
//		if (nDiff > mapOpt.nMaxMismatch)
//			return nDiff;
//
//		ref.ub >>= 1;
//		ref.lb >>= 1;
//		cmpbits++;
//		bits2--;
//
//		if (bits2 == 0) {
//			indexInWords++;
//			ref = refGenome->refInBits[indexInWords];
//			if (index + cmpbits + wordSize < refGenome->nRefSize)
//				bits2 = wordSize;
//			else
//				bits2 = refGenome->nRefSize - cmpbits - index;
//		}
//	}
//
//	return nDiff;
//}

//void Evaluation_Kernel(const MapOpt mapOpt, const CReference refGenome, const CReadInStr * d_reads, const CReadInStr * d_reads_rev,
//		CMatch * d_reads_match, SIZE_T * d_result, SIZE_T NUM) {
//	//SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
//	//SIZE_T stride = blockDim.x * gridDim.x;
//	while (i < NUM) {
//		if (d_reads_match[i].org_rev == 'o') {
//			d_result[i] = NumOfDiff(mapOpt, &d_reads[d_reads_match[i].readid], d_reads_match[i].pos, &refGenome);
//		} else {
//			d_result[i] = NumOfDiff(mapOpt, &d_reads_rev[d_reads_match[i].readid], d_reads_match[i].pos, &refGenome);
//		}
//		i += stride;
//	}
//}

//bool vector_cmp(pair<SIZE_T, SIZE_T> i, pair<SIZE_T, SIZE_T> j) {
//	return i.second > j.second;
//}

void ReadReads(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {
	FILE * fout = fopen(opt.outputFile.c_str(), "wb");
	CRead * reads, *reads_rev, *d_reads, *d_reads_rev;
	CRegion * region, *d_region;
	char * strReads;

//	MEMORY_ALLOCATE_CHECK( reads = (CRead *) malloc( MAX_MAPPING_READS * sizeof(CRead)));
//	MEMORY_ALLOCATE_CHECK( reads_rev = (CRead *) malloc( MAX_MAPPING_READS * sizeof(CRead)));
//	MEMORY_ALLOCATE_CHECK( region = (CRegion *) malloc(MAX_MAPPING_READS * TOTAL_SHIFT_rev * sizeof(CRegion)));

	HANDLE_ERROR(cudaHostAlloc((void **)&reads, sizeof(CRead) * MAX_MAPPING_READS, cudaHostAllocDefault));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads, MAX_MAPPING_READS * sizeof(CRead)));
	HANDLE_ERROR(cudaHostAlloc((void **)&reads_rev, sizeof(CRead) * MAX_MAPPING_READS, cudaHostAllocDefault));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads_rev, MAX_MAPPING_READS * sizeof(CRead)));
	HANDLE_ERROR(cudaHostAlloc((void **)&region, MAX_MAPPING_READS * TOTAL_SHIFT_rev * sizeof(CRegion), cudaHostAllocDefault));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_region, MAX_MAPPING_READS * TOTAL_SHIFT_rev * sizeof(CRegion)));
	cudaStream_t stream;
	HANDLE_ERROR(cudaStreamCreate( &stream ));

	/* read reads from the file*/
	INFO("read reads from", opt.readsFile);
	SIZE_T readsLen = ReadWholeFile(opt.readsFile, &strReads);
	cout << "log info4" << endl;
	char strRead[MAX_LINE_LEN];
	SIZE_T nReadsNum = 0;
	SIZE_T readID = 0;
	map<SIZE_T, SIZE_T> mapPosCount;

	for (SIZE_T i = 0; i < readsLen; i++) {
		SIZE_T len = GetLineFromString(&strReads[i], strRead);
		i += len;
		//cout << "log info5" << endl;
		if (strRead[0] != '>' && len != 0) {
			CHECK_READ_LEN(len, nReadsNum);
			strcpy(reads[nReadsNum].readInStr, strRead);
			reads[nReadsNum].readLen = len;
			nReadsNum++;
		}
		//cout << "log info3" << endl;
		if (nReadsNum == MAX_MAPPING_READS || (nReadsNum > 0 && i >= readsLen - 1)) {
			cout << "log info1" << endl;
			Reverse_Kernel(reads, reads_rev, nReadsNum);
			cout << "log info2" << endl;
			SIZE_T NUM = TOTAL_SHIFT_rev * nReadsNum;
			HANDLE_ERROR(cudaMemcpyAsync(&d_reads, &reads, nReadsNum * sizeof(CRead), cudaMemcpyHostToDevice, stream));
			cout << "log info5" << endl;
			HANDLE_ERROR(cudaMemcpyAsync(&d_reads_rev, &reads_rev, nReadsNum * sizeof(CRead), cudaMemcpyHostToDevice, stream));
			cout << "log info6" << endl;
			Match_Kernel<<<BLOCKS, THREADS, 0, stream>>>(refGenome, hashTable, reads, reads_rev, region, NUM);
			cout << "log info4" << endl;
			HANDLE_ERROR( cudaMemcpyAsync(region, d_region, nReadsNum * TOTAL_SHIFT_rev * sizeof(CRegion), cudaMemcpyDeviceToHost, stream));
			cout << "log info3" << endl;
			for (SIZE_T i = 0; i < NUM; i++) {
				for (SIZE_T j = region[i].lower; j <= region[i].upper; j++) {
					mapPosCount[hashTable->index[j] - i % TOTAL_SHIFT]++;
				}
				if ((i + 1) % TOTAL_SHIFT == 0) {
					for (map<SIZE_T, SIZE_T>::iterator it = mapPosCount.begin(); it != mapPosCount.end(); it++) {
						if (it->second > 1) {
							fprintf(fout, "read %d: %d %d\n", i / TOTAL_SHIFT_rev, it->first, it->second);
						}
					}
					mapPosCount.clear();
				}
			}
			readID += nReadsNum;
			nReadsNum = 0;
		}
	}
	HANDLE_ERROR(cudaStreamSynchronize(stream));
	HANDLE_ERROR(cudaStreamDestroy(stream));

	cudaFreeHost(reads);
	cudaFreeHost(reads_rev);
	cudaFreeHost(region);

	cudaFree(d_reads);
	cudaFree(d_reads_rev);
	cudaFree(d_region);

	fclose(fout);
	free(strReads);
	free(reads);
	free(region);
}

void Matching(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {

	CReference d_refGenome;
	CHashTable d_hashTable;

	SIZE_T sizeCounter = sizeof(SIZE_T) * hashTable->nSizeCounter;
	HANDLE_ERROR(cudaMalloc((void ** )&(d_hashTable.counter), sizeCounter));
	HANDLE_ERROR(cudaMemcpy(d_hashTable.counter, hashTable->counter, sizeCounter, cudaMemcpyHostToDevice));
	d_hashTable.nSizeCounter = hashTable->nSizeCounter;

	SIZE_T sizeIndex = sizeof(SIZE_T) * hashTable->nSizeIndex;
	HANDLE_ERROR(cudaMalloc((void ** )&(d_hashTable.index), sizeIndex));
	HANDLE_ERROR(cudaMemcpy(d_hashTable.index, hashTable->index, sizeIndex, cudaMemcpyHostToDevice));
	d_hashTable.nSizeIndex = hashTable->nSizeIndex;

	SIZE_T sizeRef = sizeof(InBits) * refGenome->nRefSizeInWordSize;
	HANDLE_ERROR(cudaMalloc((void ** )&(d_refGenome.refInBits), sizeRef));
	HANDLE_ERROR(cudaMemcpy(d_refGenome.refInBits, refGenome->refInBits, sizeRef, cudaMemcpyHostToDevice));
	d_refGenome.nRefSize = refGenome->nRefSize;
	d_refGenome.nRefSizeInWordSize = refGenome->nRefSizeInWordSize;

	ReadReads(opt, &d_refGenome, &d_hashTable);
	/* free memory*/
	free(refGenome->refInBits);
	free(hashTable->counter);
	free(hashTable->index);
}
