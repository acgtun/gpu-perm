#include "alloc_kernel.h"

__global__ void gpu_kernel(int nMaxMismatch, const InBits * refGenome, SIZE_T nRefSize, const SIZE_T * hashCounter, const SIZE_T * hashIndex,
		InBits * d_reads, SIZE_T * d_readsLen, SIZE_T nReads, ResultMatchedReads * d_result) {
	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
	if (i >= nReads)
		return;

	MappingOneRead(nMaxMismatch, refGenome, nRefSize, hashCounter, hashIndex, d_reads[i], d_readsLen[i], &(d_result[i]));
}

void ReadReads(const Option & opt, const InBits * refGenome, SIZE_T nRefSize, const SIZE_T * hashCounter, const SIZE_T * hashIndex) {
	FILE * fout = fopen("test_result_out.txt", "wb");
	ResultMatchedReads * result;
	SIZE_T size = MAX_MAPPING_READS * sizeof(ResultMatchedReads);
	result = (ResultMatchedReads *) malloc(size);
	if (result == NULL)
		MEMORY_ALLOCATE_ERROR
			;
	memset(result, 0x00, size);

	ReadsMatch readsMatch;
	SIZE_T readsLen = ReadWholeFile(opt.readsFile, &readsMatch.strReads);
	readsMatch.readsInBits = (InBits *) malloc(sizeof(InBits) * MAX_MAPPING_READS);
	if (readsMatch.readsInBits == NULL)
		MEMORY_ALLOCATE_ERROR
			;
	readsMatch.readsLen = (SIZE_T *) malloc(sizeof(SIZE_T) * MAX_MAPPING_READS);
	if (readsMatch.readsLen == NULL)
		MEMORY_ALLOCATE_ERROR
			;
	/////////////////////////////////////////////////////////////////////////////
	InBits * d_reads;
	SIZE_T * d_readsLen;
	HANDLE_ERROR(cudaMalloc((void **) &d_reads, MAX_MAPPING_READS * sizeof(InBits)));
	HANDLE_ERROR(cudaMalloc((void **) &d_readsLen, MAX_MAPPING_READS * sizeof(SIZE_T)));

	ResultMatchedReads * d_result;
	HANDLE_ERROR(cudaMalloc((void **) &d_result, MAX_MAPPING_READS * sizeof(ResultMatchedReads)));
	/////////////////////////////////////////////////////////////////////////////

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
			HANDLE_ERROR( cudaMemcpy(d_reads, readsMatch.readsInBits, (readsMatch.nReads) * sizeof(InBits), cudaMemcpyHostToDevice));
			HANDLE_ERROR( cudaMemcpy(d_readsLen, readsMatch.readsLen, (readsMatch.nReads) * sizeof(SIZE_T), cudaMemcpyHostToDevice));
			gpu_kernel<<<BLOCKS, THREADS>>>(opt.matchOpt.nMaxMismatch, refGenome, nRefSize, hashCounter, hashIndex, d_reads, d_readsLen,
					readsMatch.nReads, d_result);
			HANDLE_ERROR( cudaMemcpy(result, d_result, MAX_MAPPING_READS * sizeof(ResultMatchedReads), cudaMemcpyDeviceToHost));
			OutPutResult(fout, result, readsMatch.nReads, readID);
			readID += readsMatch.nReads;
			readsMatch.nReads = 0;
		}
	}
	if (readsMatch.nReads > 0) {
		HANDLE_ERROR( cudaMemcpy(d_reads, readsMatch.readsInBits, (readsMatch.nReads) * sizeof(InBits), cudaMemcpyHostToDevice));
		HANDLE_ERROR( cudaMemcpy(d_readsLen, readsMatch.readsLen, (readsMatch.nReads) * sizeof(SIZE_T), cudaMemcpyHostToDevice));
		gpu_kernel<<<BLOCKS, THREADS>>>(opt.matchOpt.nMaxMismatch, refGenome, nRefSize, hashCounter, hashIndex, d_reads, d_readsLen,
				readsMatch.nReads, d_result);
		HANDLE_ERROR( cudaMemcpy(result, d_result, MAX_MAPPING_READS * sizeof(ResultMatchedReads), cudaMemcpyDeviceToHost));
		OutPutResult(fout, result, readsMatch.nReads, readID);
		readID += readsMatch.nReads;
		readsMatch.nReads = 0;
	}

	fclose(fout);
}

void Matching(const Option & opt, const InBits * refGenome, SIZE_T nRefSize, const SIZE_T * hashCounter, const SIZE_T * hashIndex) {
	/* device variables memory allocation start*/
	SIZE_T * d_hashCounter;
	SIZE_T * d_hashIndex;
	InBits * d_refGenome;

	//printf(“%s\n”, cudaGetErrorString( cudaGetLastError() ) );
	/**************************************************************************/
	SIZE_T nRefSizeInWordSize = (nRefSize - 1) / wordSize + 1;
	SIZE_T size = (nRefSizeInWordSize + 1) * sizeof(InBits);
	HANDLE_ERROR(cudaMalloc((void **) &(d_refGenome), size));
	HANDLE_ERROR(cudaMemcpy(d_refGenome, refGenome, size, cudaMemcpyHostToDevice));
	/**************************************************************************/
	size = (NO_OF_BUCKET + 1) * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMalloc((void **) &(d_hashCounter), size));
	HANDLE_ERROR( cudaMemcpy(d_hashCounter, hashCounter, size, cudaMemcpyHostToDevice));

	size = (nRefSize + 1) * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMalloc((void **) &(d_hashIndex), size));
	HANDLE_ERROR( cudaMemcpy(d_hashIndex, hashIndex, size, cudaMemcpyHostToDevice));

	/**************************************************************************/
	ReadReads(opt, d_refGenome, nRefSize, d_hashCounter, d_hashIndex);

	cudaFree(d_hashCounter);
	cudaFree(d_hashIndex);
	cudaFree(d_refGenome);
}
