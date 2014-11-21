#include "alloc_kernel.h"

__device__ void reverseCompliment(char * strRead, const char * strReadO, int len) {
	for (int i = 0; i < len; i++) {
		strRead[i] = complimentBase(strReadO[len - i - 1]);
	}
	strRead[len] = 0;
}

__global__ void Reverse_Kernel(const CReadInStr * reads, CReadInStr * reads_rev, SIZE_T nReadsNum) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	SIZE_T stride = blockDim.x * gridDim.x;
	while (i < nReadsNum) {
		reverseCompliment(reads_rev[i].readInStr, reads[i].readInStr, reads[i].readLen);
		reads_rev[i].readLen = reads[i].readLen;
		i += stride;
	}
}

__device__ void GetShiftRead(SIZE_T shift, char * strRead, const SIZE_T len) {
	SIZE_T j = 0;
	for (SIZE_T i = shift; i < len; i++) {
		strRead[j++] = strRead[i];
	}
	strRead[j] = 0;
}

__global__ void Match_Kernel(const CReference refGenome, const CHashTable hashTable, const CReadInStr * d_readInBits,
		const CReadInStr * d_readInBits_rev, CRegion * d_reads_region, SIZE_T nReadsNum) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	SIZE_T stride = blockDim.x * gridDim.x;
	SIZE_T NUM = nReadsNum * TOTAL_SHIFT_rev;
	//printf("nReadsNum = %d\n", nReadsNum);
	//printf(" i = %u nReadsNum = %u NUM = %u\n", i, nReadsNum, NUM);
	while (i < NUM) {
		//printf("nReadsNum = %u\n", nReadsNum);
		SIZE_T readID = i / TOTAL_SHIFT_rev;
		SIZE_T readRem = i % TOTAL_SHIFT_rev;

		CReadInStr read;
		if (readRem < TOTAL_SHIFT) {
			read = d_readInBits[readID];
		} else {
			read = d_readInBits_rev[readID];
		}
		GetShiftRead(i % TOTAL_SHIFT, read.readInStr, read.readLen);
		read.readLen = read.readLen - i % TOTAL_SHIFT;
		if (read.readLen >= 22) {
			Match(&refGenome, &hashTable, &read, &(d_reads_region[i]));
		}
		i += stride;
	}
}

__device__ SIZE_T NumOfDiff(const MapOpt & mapOpt, const CReadInStr * oneRead, SIZE_T index, const CReference * refGenome) {

	InBits ref; //+22
	SIZE_T indexInWords = index / wordSize;
	SIZE_T bitsShift = index % wordSize;
	ref.ub = refGenome->refInBits[indexInWords].ub >> bitsShift;
	ref.lb = refGenome->refInBits[indexInWords].lb >> bitsShift;

	SIZE_T cmpbits = 0;
	int bits2 = wordSize - bitsShift;

	SIZE_T nDiff = 0;
	while (cmpbits < oneRead->readLen) {
		char c1 = oneRead->readInStr[cmpbits];
		char c2 = getNT((ref.ub & 0x01) << 1 | (ref.lb & 0x01));
		//printf("c2 = %c\n");

		if (c1 != c2)
			nDiff++;

		if (nDiff > mapOpt.nMaxMismatch)
			return nDiff;

		ref.ub >>= 1;
		ref.lb >>= 1;
		cmpbits++;
		bits2--;

		if (bits2 == 0) {
			indexInWords++;
			ref = refGenome->refInBits[indexInWords];
			if (index + cmpbits + wordSize < refGenome->nRefSize)
				bits2 = wordSize;
			else
				bits2 = refGenome->nRefSize - cmpbits - index;
		}
	}

	return nDiff;
}

__global__ void Evaluation_Kernel(const MapOpt mapOpt, const CReference refGenome, const CReadInStr * d_reads, const CReadInStr * d_reads_rev,
		CMatch * d_reads_match, SIZE_T * d_result, SIZE_T NUM) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	SIZE_T stride = blockDim.x * gridDim.x;
	while (i < NUM) {
		if (d_reads_match[i].org_rev == 'o') {
			d_result[i] = NumOfDiff(mapOpt, &d_reads[d_reads_match[i].readid], d_reads_match[i].pos, &refGenome);
		} else {
			d_result[i] = NumOfDiff(mapOpt, &d_reads_rev[d_reads_match[i].readid], d_reads_match[i].pos, &refGenome);
		}
		i += stride;
	}
}

void ReadReads(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {
	FILE * fout = fopen(opt.outputFile.c_str(), "wb");
	CReadInStr * reads, *d_reads;
	CReadInStr * d_reads_rev;
	CRegion * reads_region, *d_reads_region;
	CMatch * reads_match, *d_reads_match;
	SIZE_T * result, *d_result;

	HANDLE_ERROR(cudaHostAlloc((void **)&reads, sizeof(CReadInStr) * MAX_MAPPING_READS, cudaHostAllocDefault));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads, MAX_MAPPING_READS * sizeof(CReadInStr)));

	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads_rev, MAX_MAPPING_READS * sizeof(CReadInStr)));

	HANDLE_ERROR(cudaHostAlloc((void **)&reads_region, MAX_MAPPING_READS * TOTAL_SHIFT_rev * sizeof(CRegion), cudaHostAllocDefault));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads_region, MAX_MAPPING_READS * TOTAL_SHIFT_rev * sizeof(CRegion)));

	HANDLE_ERROR(cudaHostAlloc((void **)&reads_match,MAX_MAPPING_READS * 200 * sizeof(CMatch), cudaHostAllocDefault));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads_match, MAX_MAPPING_READS * 200 * sizeof(CMatch)));

	HANDLE_ERROR(cudaHostAlloc((void **)&result, MAX_MAPPING_READS * 200 * sizeof(SIZE_T), cudaHostAllocDefault));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_result, MAX_MAPPING_READS * 200 * sizeof(SIZE_T)));

	/* read reads from the file*/
	char * strReads;
	SIZE_T readsLen = ReadWholeFile(opt.readsFile.c_str(), &strReads);

	char strRead[MAX_LINE_LEN];
	SIZE_T nReadsNum = 0;
	SIZE_T readID = 0;
	cudaStream_t stream;
	HANDLE_ERROR(cudaStreamCreate( &stream ));

	for (SIZE_T i = 0; i < readsLen; i++) {
		SIZE_T len = GetLineFromString(&strReads[i], strRead);
		i += len;
		if (strRead[0] == '>')
			continue;
		if (len == 0)
			continue;
		strcpy(reads[nReadsNum].readInStr, strRead);
		reads[nReadsNum].readLen = len;
		nReadsNum++;

		if (nReadsNum == MAX_MAPPING_READS || (nReadsNum > 0 && i == readsLen - 1)) {
			HANDLE_ERROR(cudaMemcpyAsync(d_reads, reads, nReadsNum * sizeof(CReadInStr), cudaMemcpyHostToDevice, stream));
			Reverse_Kernel<<<BLOCKS, THREADS, 0, stream>>>(d_reads, d_reads_rev, nReadsNum);
			Match_Kernel<<<BLOCKS, THREADS, 0, stream>>>(*refGenome, *hashTable, d_reads, d_reads_rev, d_reads_region, nReadsNum);
			HANDLE_ERROR(
					cudaMemcpyAsync(reads_region, d_reads_region, nReadsNum * TOTAL_SHIFT_rev * sizeof(CRegion), cudaMemcpyDeviceToHost, stream));

			set<SIZE_T> pos;
			SIZE_T matchid = 0;
			for (SIZE_T i = 0; i < TOTAL_SHIFT_rev * nReadsNum; i++) {
				if (i != 0 && i % TOTAL_SHIFT == 0) {
					for (set<SIZE_T>::iterator it = pos.begin(); it != pos.end(); it++) {
						reads_match[matchid].pos = *it;
						reads_match[matchid].readid = i % TOTAL_SHIFT_rev;
						if (i % TOTAL_SHIFT_rev < TOTAL_SHIFT) {
							reads_match[matchid].org_rev = 'o';
						} else {
							reads_match[matchid].org_rev = 'r';
						}
						matchid++;
					}
					pos.clear();
				}
				for (SIZE_T j = reads_region[i].lower; j < reads_region[i].upper; j++) {
					if (pos.size() >= 100)
						break;
					pos.insert(hashTable->index[j] - i % TOTAL_SHIFT);
				}
			}

			//HANDLE_ERROR(cudaMemcpyAsync(d_reads_match, reads_match, matchid * sizeof(CMatch), cudaMemcpyHostToDevice, stream));
			//Evaluation_Kernel<<<BLOCKS, THREADS, 0, stream>>>(opt.mapOpt, *refGenome, d_reads, d_reads_rev, d_reads_match, d_result, matchid);
			//HANDLE_ERROR(cudaMemcpyAsync(reads_match, d_reads_match, matchid * sizeof(CMatch), cudaMemcpyDeviceToHost, stream));
			readID += nReadsNum;
			nReadsNum = 0;
		}
	}

	HANDLE_ERROR(cudaStreamSynchronize(stream));
	HANDLE_ERROR(cudaStreamDestroy(stream));

	fclose(fout);
	free(strReads);
	cudaFreeHost(reads);
	cudaFree(d_reads);

	cudaFree(d_reads_rev);

	cudaFreeHost(reads_region);
	cudaFree(d_reads_region);

	cudaFreeHost(reads_match);
	cudaFree(d_reads_match);

	cudaFreeHost(result);
	cudaFree(d_result);
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
	cudaFree(d_refGenome.refInBits);
	cudaFree(d_hashTable.counter);
	cudaFree(d_hashTable.index);
}
