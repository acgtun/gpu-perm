#include "alloc_kernel.h"

void cpu_kernel(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable, ReadsMatch * readsMatch,
		ResultMatchedReads * result) {
	ResultMatchedReads oneResult;
	for (SIZE_T i = 0; i < readsMatch->nReads; i++) {
		MappingOneRead(opt->matchOpt, refGenome, hashTable,
				readsMatch->readsInBits[i], readsMatch->readsLen[i],
				&oneResult);
		result[i] = oneResult;
	}
}

void ReadReads(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable) {
	FILE * fout = fopen("test_result_out.txt", "wb");
	ResultMatchedReads * result;
	SIZE_T size = (MAX_MAPPING_READS + 1) * sizeof(ResultMatchedReads);
	result = (ResultMatchedReads *) malloc(size);
	if (result == NULL)
		MEMORY_ALLOCATE_ERROR;
	memset(result, 0x00, size);

	ReadsMatch readsMatch;
	SIZE_T readsLen = ReadWholeFile(opt->readsFile, &readsMatch.strReads);
	readsMatch.readsInBits = (InBits *) malloc(
			sizeof(InBits) * (MAX_MAPPING_READS + 1));
	if (readsMatch.readsInBits == NULL)
		MEMORY_ALLOCATE_ERROR;
	readsMatch.readsLen = (SIZE_T *) malloc(
			sizeof(SIZE_T) * (MAX_MAPPING_READS + 1));
	if (readsMatch.readsLen == NULL)
		MEMORY_ALLOCATE_ERROR;

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
			cpu_kernel(opt, refGenome, hashTable, &readsMatch, result);
			OutPutResult(fout, result, readsMatch.nReads, readID);
			readID += readsMatch.nReads;
			readsMatch.nReads = 0;
		}
	}
	if (readsMatch.nReads > 0) {
		cpu_kernel(opt, refGenome, hashTable, &readsMatch, result);
		OutPutResult(fout, result, readsMatch.nReads, readID);
		readID += readsMatch.nReads;
		readsMatch.nReads = 0;
	}

	fclose(fout);
}

void Matching(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable) {
	ReadReads(opt, refGenome, hashTable);
}
