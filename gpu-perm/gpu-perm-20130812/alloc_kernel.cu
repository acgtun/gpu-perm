#include "alloc_kernel.h"

__global__ void gpu_kernel(const MapOpt mapOpt, const CReference refGenome, const CHashTable hashTable, CReadArray reads,
		CResult * result) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	if (i >= reads.nReadsNum)
		return;
	MappingOneRead(mapOpt, &refGenome, &hashTable, reads.reads[i], &(result[i]));
}

void run_kernel(FILE * fout, const MapOpt & mapOopt, const CReference * refGenome, const CHashTable * hashTable, const CReadArray * reads,
		CReadArray * d_reads, CResult * result, CResult * d_result, const SIZE_T & readID) {
	d_reads->nReadsNum = reads->nReadsNum;
	HANDLE_ERROR(cudaMemcpy(d_reads->reads, reads->reads, reads->nReadsNum * sizeof(CRead), cudaMemcpyHostToDevice));
	gpu_kernel<<<BLOCKS, THREADS>>>(mapOopt, *refGenome, *hashTable, *d_reads, d_result);
	HANDLE_ERROR(cudaMemcpy(result, d_result, reads->nReadsNum * sizeof(CResult), cudaMemcpyDeviceToHost));
	OutPutResult(fout, result, reads->nReadsNum, readID);
}

void ReadReads(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {
	FILE * fout = fopen(opt.outputFile, "wb");
	CReadArray reads, d_reads;
	MEMORY_ALLOCATE_CHECK(reads.reads = (CRead *) malloc(sizeof(CRead) * MAX_MAPPING_READS));
	HANDLE_ERROR(cudaMalloc((void **) &d_reads.reads, MAX_MAPPING_READS * sizeof(CRead)));
	/************************************************************************************/
	CResult * result, *d_result;
	HANDLE_ERROR(cudaMalloc((void **) &d_result, MAX_MAPPING_READS * sizeof(CResult)));
	MEMORY_ALLOCATE_CHECK(result = (CResult *) malloc(MAX_MAPPING_READS * sizeof(CResult)));
	memset(result, 0x00, MAX_MAPPING_READS * sizeof(CResult));
	/************************************************************************************/
	/* read reads from the file*/
	char * strReads;
	SIZE_T readsLen = ReadWholeFile(opt.readsFile, &strReads);

	char strRead[MAX_LINE_LEN];
	reads.nReadsNum = 0;
	SIZE_T readID = 0;
	for (SIZE_T i = 0; i < readsLen; i++) {
		SIZE_T len = GetLineFromString(&strReads[i], strRead);
		i += len;
		if (strRead[0] == '>')
			continue;

		EncodeRead(strRead, &reads.reads[reads.nReadsNum].readInBits, len);
		reads.reads[reads.nReadsNum].readLen = len;
		reads.nReadsNum++;

		if (reads.nReadsNum >= MAX_MAPPING_READS || (reads.nReadsNum > 0 && i == readsLen - 1)) {
			run_kernel(fout, opt.mapOpt, refGenome, hashTable, &reads, &d_reads, result, d_result, readID);
			readID += reads.nReadsNum;
			reads.nReadsNum = 0;
		}
	}

	fclose(fout);
	free(reads.reads);
	cudaFree(d_reads.reads);
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
