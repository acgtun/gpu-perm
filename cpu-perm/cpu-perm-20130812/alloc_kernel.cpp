#include "alloc_kernel.h"

void cpu_kernel(const MapOpt mapOpt, const CReference refGenome, const CHashTable hashTable, CReadArray reads, CResult * result) {
	for (SIZE_T i = 0; i < reads.nReadsNum; i++) {
		//cout << i << " " << reads.nReadsNum << endl;
		MappingOneRead(mapOpt, &refGenome, &hashTable, reads.reads[i], &(result[i]));
	}
}

//void Reverse_Kernel(CReadArray reads, CReadArray reads_rev) {
//	SIZE_T i = threadIdx.x + blockDim.x * blockIdx.x;
//	SIZE_T stride = blockDim.x * gridDim.x;
//	while (i < reads.nReadsNum) {
//		char strRead[MAX_READ_LEN];
//		reverseCompliment(strRead, reads.reads[i].readInBits, reads.reads[i].readLen);
//		EncodeRead(strRead, &reads_rev.reads[i].readInBits, reads.reads[i].readLen);
//		reads_rev.reads[i].readLen = reads.reads[i].readLen;
//		i += stride;
//	}
//}

void run_kernel(FILE * fout, const MapOpt & mapOopt, const CReference * refGenome, const CHashTable * hashTable, const CReadArray * reads,
		CResult * result, const SIZE_T & readID) {
	cpu_kernel(mapOopt, *refGenome, *hashTable, *reads, result);
	OutPutResult(fout, result, reads->nReadsNum, readID);
}

void ReadReads(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {
	FILE * fout = fopen(opt.outputFile, "wb");
	CReadArray reads;
	MEMORY_ALLOCATE_CHECK(reads.reads = (CRead *) malloc(sizeof(CRead) * MAX_MAPPING_READS));
	/************************************************************************************/
	CResult * result;
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
			cout << "run+kernel" << endl;
			run_kernel(fout, opt.mapOpt, refGenome, hashTable, &reads, result, readID);
			readID += reads.nReadsNum;
			reads.nReadsNum = 0;
		}
	}

	fclose(fout);
	free(reads.reads);
	free(strReads);
}

void Matching(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {

	ReadReads(opt, refGenome, hashTable);

	/* free memory*/
	free(refGenome->refInBits);
	free(hashTable->counter);
	free(hashTable->index);
}
