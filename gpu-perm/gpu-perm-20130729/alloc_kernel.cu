#include "alloc_kernel.h"

__global__ void gpu_kernel(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable, ReadsMatch * readsMatch,
		ResultMatchedReads * result) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	printf("i = %d\n", i);
//	printf("readsMatch->nReads = %d\n", readsMatch->nReads);
//	if (i >= readsMatch->nReads)
//		return;

	ResultMatchedReads oneResult;
	MappingOneRead(&opt->matchOpt, refGenome, hashTable,
			readsMatch->readsInBits[i], readsMatch->readsLen[i], &oneResult);
	result[i] = oneResult;
}

void ReadReads(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable) {
	FILE * fout = fopen("test_result_out.txt", "wb");

	/**************************************************************************/
	ResultMatchedReads * result;
	SIZE_T size = (MAX_MAPPING_READS + 1) * sizeof(ResultMatchedReads);
	result = (ResultMatchedReads *) malloc(size);
	if (result == NULL)
		MEMORY_ALLOCATE_ERROR
	;
	memset(result, 0x00, size);

	ReadsMatch readsMatch;
	SIZE_T readsLen = ReadWholeFile(opt->readsFile, &readsMatch.strReads);
	readsMatch.readsInBits = (InBits *) malloc(
			sizeof(InBits) * (MAX_MAPPING_READS + 1));
	if (readsMatch.readsInBits == NULL)
		MEMORY_ALLOCATE_ERROR
	;
	readsMatch.readsLen = (SIZE_T *) malloc(
			sizeof(SIZE_T) * (MAX_MAPPING_READS + 1));
	if (readsMatch.readsLen == NULL)
		MEMORY_ALLOCATE_ERROR
	;

	/**************************************************************************/
	ReadsMatch d_readsMatch;
	ResultMatchedReads * d_result;
	size = (MAX_MAPPING_READS + 1) * sizeof(InBits);
	HANDLE_ERROR(cudaMalloc((void **) &d_readsMatch.readsInBits, size));
	size = (MAX_MAPPING_READS + 1) * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMalloc((void **) &d_readsMatch.readsLen, size));

	size = (MAX_MAPPING_READS + 1) * sizeof(ResultMatchedReads);
	HANDLE_ERROR(cudaMalloc((void **) &d_result, size));
	HANDLE_ERROR(cudaMemset(d_result, 0x00, size));
	/**************************************************************************/

	int readID = 0;
	/* read reads from the file*/
	readsMatch.nReads = 0;
	char strRead[MAX_LINE_LEN];
	for (SIZE_T i = 0; i < readsLen; i++) {
		int len = GetLineFromString(&readsMatch.strReads[i], strRead);
		i += len;
		if (strRead[0] == '>')
			continue;

		EncodeRead(strRead, &readsMatch.readsInBits[readsMatch.nReads], len);
		readsMatch.readsLen[readsMatch.nReads] = len;
		readsMatch.nReads++;

		if (readsMatch.nReads >= MAX_MAPPING_READS) {
			HANDLE_ERROR(
					cudaMemcpy(d_readsMatch.readsInBits, readsMatch.readsInBits, (MAX_MAPPING_READS + 1) * sizeof(InBits), cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpy(d_readsMatch.readsLen, readsMatch.readsLen, (MAX_MAPPING_READS + 1) * sizeof(SIZE_T), cudaMemcpyHostToDevice));
			gpu_kernel<<<BLOCKS, THREADS>>>(opt, refGenome, hashTable,
					&d_readsMatch, d_result);
			//HANDLE_ERROR(
			//		cudaMemcpy(result, d_result, (MAX_MAPPING_READS + 1) * sizeof(ResultMatchedReads), cudaMemcpyDeviceToHost));
			//OutPutResult(fout, result, readsMatch.nReads, readID);
			readID += readsMatch.nReads;
			readsMatch.nReads = 0;
		}
	}
	if (readsMatch.nReads > 0) {
		HANDLE_ERROR(
				cudaMemcpy(d_readsMatch.readsInBits, readsMatch.readsInBits, (MAX_MAPPING_READS + 1) * sizeof(InBits), cudaMemcpyHostToDevice));
		HANDLE_ERROR(
				cudaMemcpy(d_readsMatch.readsLen, readsMatch.readsLen, (MAX_MAPPING_READS + 1) * sizeof(SIZE_T), cudaMemcpyHostToDevice));
		gpu_kernel<<<BLOCKS, THREADS>>>(opt, refGenome, hashTable,
				&d_readsMatch, d_result);
		//HANDLE_ERROR(
		//		cudaMemcpy(result, d_result, (MAX_MAPPING_READS + 1) * sizeof(ResultMatchedReads), cudaMemcpyDeviceToHost));
		//OutPutResult(fout, result, readsMatch.nReads, readID);
		readID += readsMatch.nReads;
		readsMatch.nReads = 0;
	}

	fclose(fout);
	cudaFree(d_readsMatch.readsInBits);
	cudaFree(d_readsMatch.readsLen);
	cudaFree(d_result);
	free(result);
	free(readsMatch.strReads);
	free(readsMatch.readsInBits);
}

void Matching(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable) {
	/* device variables memory allocation start*/

	RefGenome d_refGenome;
	HashTable d_hashTable;

	/* CUDA memory allocation start*/
	d_refGenome.nRefSize = refGenome->nRefSize;
	d_refGenome.nRefSizeInWordSize = refGenome->nRefSizeInWordSize;
	SIZE_T size = refGenome->nRefSizeInWordSize * sizeof(InBits);
	HANDLE_ERROR(cudaMalloc((void **) &d_refGenome.refInBits, size));
	HANDLE_ERROR(
			cudaMemcpy(d_refGenome.refInBits, refGenome->refInBits, size, cudaMemcpyHostToDevice));

	/**************************************************************************/
	d_hashTable.NO_OF_BUCKET = hashTable->NO_OF_BUCKET;
	size = (hashTable->NO_OF_BUCKET + 1) * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMalloc((void **) &d_hashTable.counter, size));
	HANDLE_ERROR(
			cudaMemcpy(d_hashTable.counter, hashTable->counter, size, cudaMemcpyHostToDevice));

	size = (refGenome->nRefSize + 1) * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMalloc((void **) &d_hashTable.index, size));
	HANDLE_ERROR(
			cudaMemcpy(d_hashTable.index, hashTable->index, size, cudaMemcpyHostToDevice));

	/**************************************************************************/
	ReadReads(opt, &d_refGenome, &d_hashTable);
	cudaFree(d_hashTable.counter);
	cudaFree(d_hashTable.index);
	cudaFree(d_refGenome.refInBits);
}
