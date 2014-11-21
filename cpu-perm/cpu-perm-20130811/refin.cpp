#include "refin.h"

SIZE_T SEED_SPACE[] = { 0, 1, 2, 4, 7, 8, 9, 11, 14, 15, 16, 18, 21 };

int RemoveNonACGTNBase(char * strRef, SIZE_T refLen,
		set<SIZE_T> & setConsectiveN) {
	/* This function removes all non-ACGTN characters */
	char strRet[MAX_LINE_LEN];
	SIZE_T j = 0;
	for (SIZE_T i = 0; i < refLen; i++) {
		if (strRef[i] == '>') {
			i += GetLineFromString(&strRef[i], strRet);
		} else if (isACGT(strRef[i]) || strRef[i] == 'N' || strRef[i] == 'n') {
			strRef[j++] = toupper(strRef[i]);
		}
	}
	strRef[j] = 0;

	/* delete the index whose seed pattern has 'N' */
	for (SIZE_T i = 0; i < j; i++) {
		if (i % NUMBER_OF_SPACE != 0)
			continue;
		for (SIZE_T k = 0; k < 13; k++) {
			if (strRef[i + SEED_SPACE[k]] == 'N') {
				setConsectiveN.insert(i);
				break;
			}
		}
	}

	/* change all 'N' to A,C,G,T with the same probability*/
	srand(time(NULL));
	int r = 0;
	for (SIZE_T i = 0; i < j; i++) {
		if (strRef[i] == 'N') {
			r = rand() % 4;
			strRef[i] = getNT(r);
			cout << "strRef[i] = " << strRef[i] << endl;
		}
	}

	return j;
}

void genTestData(char * strRef, SIZE_T len) {
	FILE * fread = fopen("testread_chr1.fa", "wb");
	FILE * fans = fopen("ans.txt", "wb");

	srand(time(NULL));

	for (SIZE_T i = 0; i < 1000000; i++) {
		int r = rand() % len;
		fprintf(fread, ">read%d\n", i);
		for (SIZE_T j = r, n = 0; j < len && n < 50; n++, j++) {
			fprintf(fread, "%c", strRef[j]);
		}
		fprintf(fread, "\n");
		fprintf(fans, "read%d %d\n", i, r);
	}
	fclose(fread);
	fclose(fans);
}

void RefEncodeToBits(CReference * refGenome, const char * strRef) {
	/* This function transfers the genome to binary format. The binary format genome
	 * is stored as many InBits, which has two WORD_SIZE. For each character A(00),C(01),
	 * G(10),T(11), the upper bit stores in ub, and lower bit stores in lb.
	 * */
	cout << "1" << endl;
	refGenome->nRefSizeInWordSize = (refGenome->nRefSize - 1) / wordSize + 1;
	MEMORY_ALLOCATE_CHECK(
			refGenome->refInBits = (InBits * ) malloc(sizeof(InBits) * refGenome->nRefSizeInWordSize));
	cout << "4" << endl;
	char strReads[wordSize + 1];
	for (SIZE_T i = 0; i < refGenome->nRefSizeInWordSize - 1; i++) {
		memcpy(&strReads, &(strRef[i * wordSize]), wordSize);
		strReads[wordSize] = 0;
		EncodeRead(strReads, &(refGenome->refInBits[i]), wordSize);
	}
	cout << "2" << endl;
	SIZE_T codesize = (refGenome->nRefSizeInWordSize - 1) * wordSize;
	SIZE_T remSize = refGenome->nRefSize - codesize;
	memcpy(strReads, &(strRef[codesize]), (SIZE_T) remSize);
	strReads[remSize] = 0;
	EncodeRead(strReads,
			&(refGenome->refInBits[refGenome->nRefSizeInWordSize - 1]),
			remSize);
	cout << "3" << endl;
}

void GetReference(CReference * refGenome, Option & opt) {
	LOG_INFO;
	char * strRef;
	cout << "hah1" << endl;
	SIZE_T refLen = ReadWholeFile(opt.refFile, &strRef);
	cout << "hah2" << endl;
	cout << "refLen = " << refLen << endl;
	refGenome->nRefSize = RemoveNonACGTNBase(strRef, refLen,
			opt.setConsectiveN);
	cout << refGenome->nRefSize << endl;
	//genTestData(strRef, refGenome->nRefSize);
	RefEncodeToBits(refGenome, strRef);
	cout << " opt.setConsectiveN.size() = " << opt.setConsectiveN.size()
			<< endl;
	cout << "size = " << refGenome->nRefSize << endl;
	free(strRef);
}
