#include "alloc_kernel.h"

void Match_Kernel(const CReference * refGenome, const CHashTable * hashTable, const CReadInStr * readInBits, const CReadInStr * readInBits_rev,
		CRegion * d_reads_region, const SIZE_T & nReadsNum) {
	SIZE_T NUM = nReadsNum * TOTAL_SHIFT_rev;
	for (SIZE_T i = 0; i < NUM; i++) {
		if (i % TOTAL_SHIFT_rev < TOTAL_SHIFT) {
			Match(refGenome, hashTable, &readInBits[i / TOTAL_SHIFT_rev], i % TOTAL_SHIFT, &(d_reads_region[i]));
		} else {
			Match(refGenome, hashTable, &readInBits_rev[i / TOTAL_SHIFT_rev], i % TOTAL_SHIFT, &(d_reads_region[i]));
		}
	}
}

SIZE_T NumOfDiff(const MapOpt & mapOpt, const CReadInStr * oneRead, SIZE_T index, const CReference * refGenome) {

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
	CReadInStr * reads, *reads_rev;
	CRegion * region;
	char * strReads;

	MEMORY_ALLOCATE_CHECK( reads = (CReadInStr *) malloc( MAX_MAPPING_READS * sizeof(CReadInStr)));
	MEMORY_ALLOCATE_CHECK( reads_rev = (CReadInStr *) malloc( MAX_MAPPING_READS * sizeof(CReadInStr)));
	MEMORY_ALLOCATE_CHECK( region = (CRegion *) malloc(MAX_MAPPING_READS * TOTAL_SHIFT_rev * sizeof(CRegion)));

	/* read reads from the file*/
	INFO("read reads from", opt.readsFile);
	SIZE_T readsLen = ReadWholeFile(opt.readsFile, &strReads);

	char strRead[MAX_LINE_LEN];
	SIZE_T nReadsNum = 0;
	SIZE_T readID = 0;
	map<SIZE_T, SIZE_T> mapPosCount;

	for (SIZE_T i = 0; i < readsLen; i++) {
		SIZE_T len = GetLineFromString(&strReads[i], strRead);
		i += len;
		if (strRead[0] != '>' && len != 0) {
			CHECK_READ_LEN(len, nReadsNum);
			strcpy(reads[nReadsNum].readInStr, strRead);
			reads[nReadsNum].readLen = len;
			nReadsNum++;
		}

		if (nReadsNum == MAX_MAPPING_READS || (nReadsNum > 0 && i >= readsLen - 1)) {
			Reverse_Kernel(reads, reads_rev, nReadsNum);
			TIME_INFO(Match_Kernel(refGenome, hashTable, reads, reads_rev, region, nReadsNum), "match time");
			SIZE_T NUM = TOTAL_SHIFT_rev * nReadsNum;
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

	fclose(fout);
	free(strReads);
	free(reads);
	free(region);
}

void Matching(const Option & opt, const CReference * refGenome, const CHashTable * hashTable) {

	ReadReads(opt, refGenome, hashTable);
	/* free memory*/
	free(refGenome->refInBits);
	free(hashTable->counter);
	free(hashTable->index);
}
