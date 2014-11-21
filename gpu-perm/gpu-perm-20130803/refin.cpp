#include "refin.h"

int RemoveNonACGTNBase(char * strRef, int refLen) {
	/* This function removes all non-ACGTN characters */
	FILE * frand = fopen("rand1.txt", "wb");
	int cnt = 0;
	char strRet[MAX_LINE_LEN];
	int j = 0;
	srand(time(NULL));
	for (int i = 0; i < refLen; i++) {
		if (strRef[i] == '>') {
			i += GetLineFromString(&strRef[i], strRet);
		} else if (isACGT(strRef[i])) {
			strRef[j++] = toupper(strRef[i]);
		} else {
			int r = rand() % 4;
			fprintf(frand, "%c ", strRef[i]);
			strRef[j++] = getNT(r);
			fprintf(frand, "%d\n", r);
			cnt++;
		}
	}
	fprintf(frand, "cnt = %d\n", cnt);
	fprintf(frand, "j = %d\n", j);
	fclose(frand);
	strRef[j] = 0;
	return j;
}

void genTestData(char * strRef, int len) {
	FILE * fread = fopen("testread_chr1.fa", "wb");
	FILE * fans = fopen("ans.txt", "wb");

	srand(time(NULL));

	for (int i = 0; i < 1000000; i++) {
		int r = rand() % len;
		fprintf(fread, ">read%d\n", i);
		for (int j = r, n = 0; j < len && n < 50; n++, j++) {
			fprintf(fread, "%c", strRef[j]);
		}
		fprintf(fread, "\n");
		fprintf(fans, "read%d %d\n", i, r);
	}
	fclose(fread);
	fclose(fans);
}

void RefEncodeToBits(CReference * refGenome, char * strRef) {
	/* This function transfers the genome to binary format. The binary format genome
	 * is stored as many InBits, which has two WORD_SIZE. For each character A(00),C(01),
	 * G(10),T(11), the upper bit stores in ub, and lower bit stores in lb.
	 * */
	refGenome->nRefSizeInWordSize = (refGenome->nRefSize - 1) / wordSize + 1;
	MEMORY_ALLOCATE_CHECK(refGenome->refInBits = (InBits * ) malloc(sizeof(InBits) * (refGenome->nRefSizeInWordSize + 1)));
	char strReads[wordSize + 1];
	for (SIZE_T i = 0; i < refGenome->nRefSizeInWordSize - 1; i++) {
		memcpy(&strReads, &(strRef[i * wordSize]), wordSize);
		strReads[wordSize] = 0;
		EncodeRead(strReads, &(refGenome->refInBits[i]), wordSize);
	}
	int codesize = (refGenome->nRefSizeInWordSize - 1) * wordSize;
	int remSize = refGenome->nRefSize - codesize;
	memcpy(strReads, &(strRef[codesize]), (SIZE_T) remSize);
	strReads[remSize] = 0;
	EncodeRead(strReads, &(refGenome->refInBits[refGenome->nRefSizeInWordSize - 1]), remSize);
	free(strRef);
}

void GetReference(CReference * refGenome, const Option & opt) {
	LOG_INFO;
	char * strRef;
	SIZE_T refLen = ReadWholeFile(opt.refFile, &strRef);
	refGenome->nRefSize = RemoveNonACGTNBase(strRef, refLen);
	genTestData(strRef, refGenome->nRefSize);
	RefEncodeToBits(refGenome, strRef);
}
