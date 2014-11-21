#include "alloc_kernel.h"

__device__ void reverseCompliment(char * strRead, char * strReadO, int len) {
	Swap(strReadO, len);
	for (int i = 0; i < len; i++) {
		strRead[i] = complimentBase(strReadO[i]);
	}
}

__global__ void Reverse_Kernel(CReadInStr * reads, CReadInStr * reads_rev, SIZE_T nReadsNum) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	SIZE_T stride = blockDim.x * gridDim.x;
	while (i < nReadsNum) {
		reverseCompliment(reads_rev[i].readInStr, reads[i].readInStr, reads[i].readLen);
		reads_rev[i].readLen = reads[i].readLen;
		i += stride;
	}
}

__device__ void StrEncodeToBits(CReadInBits * readInBits, const char * strRef, int len) {
	readInBits->readLen = len;
	SIZE_T nReadSizeInWordSize = (len - 1) / wordSize + 1;

	char strReads[wordSize + 1];
	for (SIZE_T i = 0; i < nReadSizeInWordSize - 1; i++) {
		memcpy(&strReads, &(strRef[i * wordSize]), wordSize);
		strReads[wordSize] = 0;
		EncodeRead(strReads, &(readInBits->readInBits[i]), wordSize);
	}

	SIZE_T codesize = (nReadSizeInWordSize - 1) * wordSize;
	SIZE_T remSize = len - codesize;
	memcpy(strReads, &(strRef[codesize]), (SIZE_T) remSize);
	strReads[remSize] = 0;
	EncodeRead(strReads, &(readInBits->readInBits[nReadSizeInWordSize - 1]), remSize);
}

__global__ void ReadToBits_Kernel(const CReadInStr * reads, CReadInBits * readInBits, SIZE_T nReadsNum) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	SIZE_T stride = blockDim.x * gridDim.x;
	while (i < nReadsNum) {
		StrEncodeToBits(&(readInBits[i]), reads[i].readInStr, reads[i].readLen);
		i += stride;
	}
}

__device__ void GetShiftRead(SIZE_T shift, const char * strRead, const SIZE_T len, char * strReadShift) {
	SIZE_T j = 0;
	for (SIZE_T i = shift; i < len; i++) {
		strReadShift[j] = strRead[i];
	}
	strReadShift[j] = 0;
}

__global__ void Shfit_Kernel(const CReadInStr * reads, CReadInStr * d_reads_shift, SIZE_T nReadsNum) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	SIZE_T stride = blockDim.x * gridDim.x;
	while (i < nReadsNum) {
		for (SIZE_T j = 0; j < TOTAL_NEED_SHIFT; j++) {
			GetShiftRead(j, reads[i].readInStr, reads[i].readLen, d_reads_shift[TOTAL_NEED_SHIFT * i + j].readInStr);
			d_reads_shift[TOTAL_NEED_SHIFT * i + j].readLen = reads[i].readLen - j;
		}
		i += stride;
	}
}

__global__ void Match_Kernel(const CReference refGenome, const CHashTable hashTable, CReadInBits * d_readInBits, CRegion * d_reads_region,
		SIZE_T nReadsNum) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	SIZE_T stride = blockDim.x * gridDim.x;
	while (i < nReadsNum) {
		Match(&refGenome, &hashTable, &(d_readInBits[i]), &(d_reads_region[i]));
		i += stride;
	}
}

void run_kernel(FILE * fout, const MapOpt & mapOopt, const CReference * refGenome, const CHashTable * hashTable, CReadInStr * d_reads,
		const SIZE_T & nReadsNum, const SIZE_T & readID) {
	CReadInStr * d_reads_rev;
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads_rev, nReadsNum * sizeof(CReadInStr)));
	Reverse_Kernel<<<BLOCKS, THREADS>>>(d_reads, d_reads_rev, nReadsNum);
	LOG_INFO_CPP;
	CReadInStr * d_reads_shift, *d_read_rev_shfit;
	cout << sizeof(CReadInStr) << endl;
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads_shift, nReadsNum * TOTAL_NEED_SHIFT * sizeof(CReadInStr)));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_read_rev_shfit, nReadsNum * TOTAL_NEED_SHIFT * sizeof(CReadInStr)));
	Shfit_Kernel<<<BLOCKS, THREADS>>>(d_reads, d_reads_shift, nReadsNum);
	Shfit_Kernel<<<BLOCKS, THREADS>>>(d_reads_rev, d_read_rev_shfit, nReadsNum);
	cudaFree(d_reads);
	cudaFree(d_reads_rev);
	LOG_INFO_CPP;

	cout << nReadsNum * TOTAL_NEED_SHIFT << endl;
	cout << nReadsNum << endl;
	cout << TOTAL_NEED_SHIFT << endl;

	CReadInBits * d_readInBits, *d_readInBits_rev;
	cout << sizeof(CReadInBits) << endl;
	HANDLE_ERROR(cudaMalloc((void ** ) &d_readInBits, nReadsNum * TOTAL_NEED_SHIFT * sizeof(CReadInBits)));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_readInBits_rev, nReadsNum * TOTAL_NEED_SHIFT * sizeof(CReadInBits)));
	ReadToBits_Kernel<<<BLOCKS, THREADS>>>(d_reads_shift, d_readInBits, nReadsNum * TOTAL_NEED_SHIFT);
	ReadToBits_Kernel<<<BLOCKS, THREADS>>>(d_read_rev_shfit, d_readInBits_rev, nReadsNum * TOTAL_NEED_SHIFT);
	cudaFree(d_reads_shift);
	cudaFree(d_read_rev_shfit);
	LOG_INFO_CPP;

	CRegion * d_reads_region, *d_reads_rev_region;
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads_region, nReadsNum * TOTAL_NEED_SHIFT * sizeof(CRegion)));
	HANDLE_ERROR(cudaMalloc((void ** ) &d_reads_rev_region, nReadsNum * TOTAL_NEED_SHIFT * sizeof(CRegion)));
	Match_Kernel<<<BLOCKS, THREADS>>>(*refGenome, *hashTable, d_readInBits, d_reads_region, nReadsNum * TOTAL_NEED_SHIFT);
	Match_Kernel<<<BLOCKS, THREADS>>>(*refGenome, *hashTable, d_readInBits_rev, d_reads_rev_region, nReadsNum * TOTAL_NEED_SHIFT);
	cudaFree(d_reads_shift);
	cudaFree(d_read_rev_shfit);
	LOG_INFO_CPP;

//HANDLE_ERROR(cudaMemcpy(result, d_result, reads->nReadsNum * sizeof(CResult), cudaMemcpyDeviceToHost));
//OutPutResult(fout, result, reads->nReadsNum, readID);
}

void ReadReads(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {
	FILE * fout = fopen(opt.outputFile.c_str(), "wb");
	CReadInStr * reads, *d_reads;
	MEMORY_ALLOCATE_CHECK(reads = (CReadInStr *) malloc(sizeof(CReadInStr) * MAX_MAPPING_READS));
	/************************************************************************************/
//	CResult * result, *d_result;
//	HANDLE_ERROR(cudaMalloc((void **) &d_result, MAX_MAPPING_READS * sizeof(CResult)));
//	MEMORY_ALLOCATE_CHECK(result = (CResult *) malloc(MAX_MAPPING_READS * sizeof(CResult)));
//	memset(result, 0x00, MAX_MAPPING_READS * sizeof(CResult));
	/************************************************************************************/
	/* read reads from the file*/
	char * strReads;
	SIZE_T readsLen = ReadWholeFile(opt.readsFile.c_str(), &strReads);

	char strRead[MAX_LINE_LEN];
	SIZE_T nReadsNum = 0;
	SIZE_T readID = 0;
	for (SIZE_T i = 0; i < readsLen; i++) {
		SIZE_T len = GetLineFromString(&strReads[i], strRead);
		i += len;
		if (strRead[0] == '>')
			continue;

		//EncodeRead(strRead, &reads.reads[reads.nReadsNum].readInBits, len);
		//StrEncodeToBits(&(reads.reads[reads.nReadsNum]), strRead, len);
		//reads.reads[reads.nReadsNum].readLen = len;
		strcpy(reads[nReadsNum].readInStr, strRead);
		reads[nReadsNum].readLen = len;
		nReadsNum++;

		if (nReadsNum >= MAX_MAPPING_READS || (nReadsNum > 0 && i == readsLen - 1)) {
			HANDLE_ERROR(cudaMalloc((void **) &d_reads, MAX_MAPPING_READS * sizeof(CReadInStr)));
			HANDLE_ERROR(cudaMemcpy(d_reads, reads, nReadsNum * sizeof(CReadInStr), cudaMemcpyHostToDevice));
			run_kernel(fout, opt.mapOpt, refGenome, hashTable, d_reads, nReadsNum, readID);
			readID += nReadsNum;
			nReadsNum = 0;
		}
	}

	fclose(fout);
	free(reads);
	cudaFree(d_reads);
	free(strReads);
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
