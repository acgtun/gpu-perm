#include "refin.h"

int RemoveNonACGTNBase(char * strRef, int refLen) {
	char strRet[MAX_LINE_LEN];
	int j = 0;
	for (int i = 0; i < refLen; i++) {
		if (strRef[i] == '>') {
			i += GetLineFromString(&strRef[i], strRet);
		} else if (isACGT(strRef[i]) || strRef[i] == 'N' || strRef[i] == 'n') {
			strRef[j++] = toupper(strRef[i]);
		}
	}
	strRef[j] = 0;
	return j;
}

void RefEncodeToBits(RefGenome * refGenome) {
	refGenome->nRefSizeInWordSize = (refGenome->nRefSize - 1) / wordSize + 1;
	refGenome->refInBits = (InBits *) malloc(
			sizeof(InBits) * (refGenome->nRefSizeInWordSize + 1));
	cout << "siz1e =  " << refGenome->nRefSizeInWordSize * sizeof(InBits)
			<< endl;
	if (refGenome->refInBits == NULL)
		MEMORY_ALLOCATE_ERROR;

	char strReads[wordSize + 1];
	for (SIZE_T i = 0; i < refGenome->nRefSizeInWordSize - 1; i++) {
		memcpy(&strReads, &(refGenome->strRef[i * wordSize]), wordSize);
		strReads[wordSize] = 0;
		EncodeRead(strReads, &refGenome->refInBits[i], wordSize);
	}
	int codesize = (refGenome->nRefSizeInWordSize - 1) * wordSize;
	int remSize = refGenome->nRefSize - codesize;
	memcpy(strReads, &(refGenome->strRef[codesize]), (SIZE_T) remSize);
	strReads[remSize] = 0;
	EncodeRead(strReads,
			&(refGenome->refInBits[refGenome->nRefSizeInWordSize - 1]),
			remSize);
}

void GetReference(RefGenome * refGenome, const Option * opt) {
	SIZE_T refLen = ReadWholeFile(opt->refFile, &refGenome->strRef);
	refGenome->nRefSize = RemoveNonACGTNBase(refGenome->strRef, refLen);
	FILE * fout = fopen("reference_out.txt", "wb");
	cout << "len = " << " " << refGenome->nRefSize << endl;
	for (size_t i = 0; i < refGenome->nRefSize; i++) {

		fprintf(fout, "%09d %c\n", i + 1, refGenome->strRef[i]);
	}
	fclose(fout);
	RefEncodeToBits(refGenome);
	free(refGenome->strRef);
}
