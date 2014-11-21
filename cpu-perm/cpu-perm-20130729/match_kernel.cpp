#include "match_kernel.h"
#include "seed.h"

FILE * fsee;

InBits GetShiftRead(SIZE_T shift, InBits read) {
	InBits r;
	r.ub = read.ub >> shift;
	r.lb = read.lb >> shift;
	return r;
}

int cuSTRCMP(char * str1, int len1, char * str2, int len2) {
	int i = 0;
	while (i < len1 && i < len2) {
		if (str1[i] < str2[i])
			return -1;
		else if (str1[i] > str2[i])
			return 1;
		i++;
	}
	if (i == len1 && i == len2)
		return 0;
	if (i == len2)
		return 1;
	else
		return -1;
}

int CMP(char * strRead, int len, SIZE_T index, const RefGenome * d_refGenome) {
	InBits r, v;
	SIZE_T s = GetKmer(d_refGenome, index, wordSize, &r);
	int ss = GetF2SeedForBits(r, s, &v);

	char strRef64[MAX_READ_LEN];
	DecodeReadReverse(strRef64, ss, &v);

	return cuSTRCMP(strRead, len, strRef64, ss);
}

SIZE_T LowerBound(SIZE_T low, SIZE_T high, char * strRead, int s,
		const RefGenome * d_refGenome, const HashTable * d_hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high) / 2;
		if (CMP(strRead, s, d_hashTable->index[mid], d_refGenome) <= 0)
			high = mid;
		else
			low = mid + 1;
	}
	return low;
}

SIZE_T UpperBound(SIZE_T low, SIZE_T high, char * strRead, int s,
		const RefGenome * d_refGenome, const HashTable * d_hashTable) {
	SIZE_T mid = 0;
	while (low < high) {
		mid = (low + high + 1) / 2;
		if (CMP(strRead, s, d_hashTable->index[mid], d_refGenome) >= 0)
			low = mid;
		else
			high = mid - 1;
	}
	return low;
}

void MappingOneRead(const MatchOpt * matchOpt, const RefGenome * refGenome,
		const HashTable * hashTable, InBits readInBits, int readLen,
		ResultMatchedReads * oneResult) {
	for (int i = 0; i <= NUMBER_OF_SHIFT; i++) {
		InBits read = GetShiftRead(i, readInBits);
		SIZE_T hashValue = GetHashValue(read);
		fprintf(fsee, "hashValue = %d ", hashValue);
		InBits ret;
		int len = readLen - i;
		int s = GetF2SeedForBits(read, len, &ret);
		char strRead[MAX_READ_LEN];
		DecodeReadReverse(strRead, s, &ret);

		SIZE_T l = hashTable->counter[hashValue];
		SIZE_T u = hashTable->counter[hashValue + 1] - 1;

		SIZE_T lower = LowerBound(l, u, strRead, s, refGenome, hashTable);
		SIZE_T upper = UpperBound(l, u, strRead, s, refGenome, hashTable);
		fprintf(fsee, " (%d,%d)\n", lower, upper);
		for (SIZE_T j = lower; j <= upper; j++) {
			int s = GetKmer(refGenome, hashTable->index[j] - i, readLen, &ret);
			if (s != readLen)
				continue;
			SIZE_T nDiff = bitsStrNCompare(ret, readInBits, readLen);
			if (nDiff <= matchOpt->nMaxMismatch) {
				oneResult->nMismatch[oneResult->nRet] = nDiff;
				oneResult->nStartPos[oneResult->nRet] = hashTable->index[j] - i;
				//printf("readID = %d %d %d\n", hashTable->index[j], i);
				oneResult->nRet++;
			}
		}
	}
}

void CPU_Matching(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable) {
	ReadsMatch readsMatch;
	ResultMatchedReads oneResult;
	InBits readInBits, readInBits_rev;
	vector<ResultMatchedReads> result;

	fsee = fopen("see.txt", "wb");

	SIZE_T readsLen = ReadWholeFile(opt->readsFile, &readsMatch.strReads);
	char strRead[MAX_LINE_LEN];
	int readID = 0;
	for (SIZE_T i = 0; i < readsLen; i++) {
		//cout << "i = " << i << "readsLen = " << readsLen << endl;
		int len = GetLineFromString(&readsMatch.strReads[i], strRead);
		i += len;
		if (strRead[0] == '>')
			continue;
		EncodeRead(strRead, &readInBits, len);
		//cout << "strRead = " << strRead << endl;
		//printWORD(readInBits.ub, len);
		//printWORD(readInBits.lb, len);
		fprintf(fsee, "readID = %d ", ++readID);
		fprintf(fsee, "%s\n", strRead);
		oneResult.nRet = 0;
		MappingOneRead(&opt->matchOpt, refGenome, hashTable, readInBits, len,
				&oneResult);
		reverseCompliment(strRead, len);
		EncodeRead(strRead, &readInBits_rev, len);
		MappingOneRead(&opt->matchOpt, refGenome, hashTable, readInBits_rev,
				len, &oneResult);
		fprintf(fsee, "\n");
		result.push_back(oneResult);
		for (size_t j = 0; j < oneResult.nRet; j++) {
			fprintf(fsee, " (%u, %u)", oneResult.nStartPos[j],
					oneResult.nMismatch[j]);
		}
		fprintf(fsee, "\n");
	}
	FILE * fout = fopen("test_result_out.txt", "wb");
	fprintf(fout, "%u\n", result.size());
	for (size_t i = 0; i < result.size(); i++) {
		fprintf(fout, "read %d:", i);
		for (size_t j = 0; j < result[i].nRet; j++) {
			fprintf(fout, " (%u, %u)", result[i].nStartPos[j],
					result[i].nMismatch[j]);
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
	fclose(fsee);
}
