#include "allocate.h"

void MemoryAlloocate(RefGenome * d_refGenome, HashTable * d_hashTable,
		ReadsMatch * d_readsMatch, const RefGenome * refGenome,
		const HashTable * hashTable) {
	d_refGenome->nRefSize = refGenome->nRefSize;
	d_refGenome->nRefSizeInWordSize = refGenome->nRefSizeInWordSize;
	int size = refGenome->nRefSizeInWordSize * sizeof(InBits);
	HANDLE_ERROR(cudaMalloc((void **) &d_refGenome->refInBits, size));
	HANDLE_ERROR(
			cudaMemcpy(d_refGenome->refInBits, refGenome->refInBits, size, cudaMemcpyHostToDevice));
	d_hashTable->NO_OF_BUCKET = hashTable->NO_OF_BUCKET;
	size = (hashTable->NO_OF_BUCKET + 1) * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMalloc((void **) &d_hashTable->counter, size));
	HANDLE_ERROR(
			cudaMemcpy(d_hashTable->counter, hashTable->counter, size, cudaMemcpyHostToDevice));

	size = (refGenome->nRefSize + 1) * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMalloc((void **) &d_hashTable->index, size));
	HANDLE_ERROR(
			cudaMemcpy(d_hashTable->index, hashTable->index, size, cudaMemcpyHostToDevice));

	size = (MAX_MAPPING_READS + 1) * sizeof(InBits);
	HANDLE_ERROR(cudaMalloc((void **) &d_readsMatch->readsInBits, size));
	size = (MAX_MAPPING_READS + 1) * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMalloc((void **) &d_readsMatch->readsLenAretN, size));
}

void OUTPUT(FILE * fout, const ResultMatchedReads * result,
		const ReadsMatch * readsMatch) {
	for (SIZE_T i = 0; i < readsMatch->nReads; i++) {
		fprintf(fout, "read %d:", i + 1);
		for (SIZE_T j = 0; j < readsMatch->readsLenAretN[i]; j++) {
			fprintf(fout, " (%u, %u)", result[i].nStartPos[j],
					result[i].nMismatch[j]);
		}
		fprintf(fout, "\n");
	}
}

void run_kernel(FILE * fout, ReadsMatch * readsMatch,
		ResultMatchedReads * result, ReadsMatch * d_readsMatch,
		MatchOpt * d_matchOpt, RefGenome * d_refGenome, HashTable * d_hashTable,
		ResultMatchedReads * d_result) {
	d_readsMatch->nReads = readsMatch->nReads;
	SIZE_T size = readsMatch->nReads * sizeof(InBits);
	HANDLE_ERROR(cudaMemset(d_readsMatch->readsInBits, 0x00, size));
	HANDLE_ERROR(
			cudaMemcpy(d_readsMatch->readsInBits, readsMatch->readsInBits, size, cudaMemcpyHostToDevice));
	size = readsMatch->nReads * sizeof(SIZE_T);
	HANDLE_ERROR(cudaMemset(d_readsMatch->readsLenAretN, 0x00, size));
	HANDLE_ERROR(
			cudaMemcpy(d_readsMatch->readsLenAretN, readsMatch->readsLenAretN, size, cudaMemcpyHostToDevice));

	CUDA_Mapping<<<BLOCKS, THREADS>>>(d_matchOpt, d_refGenome, d_hashTable,
			d_readsMatch, d_result);
	//cudaDeviceSynchronize();
	size = (MAX_MAPPING_READS + 1) * sizeof(ResultMatchedReads);
	memset(result, 0x00, size);
	HANDLE_ERROR(cudaMemcpy(result, d_result, size, cudaMemcpyDeviceToHost));
	OUTPUT(fout, result, readsMatch);
	readsMatch->nReads = 0;
}

void Matching(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable) {

	if (opt->matchOpt.nIsGPU == 0) {
		CPU_Matching(opt, refGenome, hashTable);
		return;
	}
	/* device variables memory allocation start*/
	MatchOpt d_matchOpt;
	RefGenome d_refGenome;
	HashTable d_hashTable;
	ReadsMatch d_readsMatch;
	ResultMatchedReads * d_result;
	ResultMatchedReads * result;

	FILE * fout = fopen("test/test_out.txt", "wb");

	/* CUDA memory allocation start*/
	SIZE_T size = (MAX_MAPPING_READS + 1) * sizeof(ResultMatchedReads);
	HANDLE_ERROR(cudaMalloc((void **) &d_result, size));
	HANDLE_ERROR(cudaMemset(d_result, 0x00, size));

	result = (ResultMatchedReads *) malloc(size);
	if (result == NULL)
		MEMORY_ALLOCATE_ERROR
	;
	memset(result, 0x00, size);

	d_matchOpt = opt->matchOpt;
	MemoryAlloocate(&d_refGenome, &d_hashTable, &d_readsMatch, refGenome,
			hashTable);

	/* CUDA memory allocation end*/
	/* Input reads and mapping*/
	ReadsMatch readsMatch;
	SIZE_T readsLen = ReadWholeFile(opt->readsFile, &readsMatch.strReads);
	readsMatch.readsInBits = (InBits *) malloc(
			sizeof(InBits) * (MAX_MAPPING_READS + 1));
	if (readsMatch.readsInBits == NULL)
		MEMORY_ALLOCATE_ERROR
	;

	readsMatch.readsLenAretN = (SIZE_T *) malloc(
			sizeof(SIZE_T) * (MAX_MAPPING_READS + 1));
	if (readsMatch.readsLenAretN == NULL)
		MEMORY_ALLOCATE_ERROR
	;

	/*allocate read to gpu and mapping*/
	readsMatch.nReads = 0;
	char strRead[MAX_LINE_LEN];
	for (int i = 0; i < readsLen; i++) {
		int len = GetLineFromString(&readsMatch.strReads[i], strRead);
		i += len;
		if (strRead[0] == '>')
			continue;
		EncodeRead(strRead, &readsMatch.readsInBits[readsMatch.nReads], len);
		readsMatch.readsLenAretN[readsMatch.nReads] = len;
		readsMatch.nReads++;
		if (readsMatch.nReads >= MAX_MAPPING_READS) {
			run_kernel(fout, &readsMatch, result, &d_readsMatch, &d_matchOpt,
					&d_refGenome, &d_hashTable, d_result);
		}
	}
	if (readsMatch.nReads > 0) {
		run_kernel(fout, &readsMatch, result, &d_readsMatch, &d_matchOpt,
				&d_refGenome, &d_hashTable, d_result);
	}
	fclose(fout);
	printf("END HERE\n");
	/* host memory release */
	free(readsMatch.strReads);
	free(readsMatch.readsInBits);
	free(hashTable->index);
	free(hashTable->counter);
	free(refGenome->refInBits);
	/* CUDA memory release */
	cudaFree(d_refGenome.refInBits);
	cudaFree(d_hashTable.counter);
	cudaFree(d_hashTable.index);
	cudaFree(d_readsMatch.readsInBits);
	cudaFree(d_result);
}
