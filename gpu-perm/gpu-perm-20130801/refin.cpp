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

void RefEncodeToBits(InBits ** refGenome, SIZE_T nRefSize, char * strRef) {
	SIZE_T nRefSizeInWordSize = (nRefSize - 1) / wordSize + 1;
	*refGenome = (InBits *) malloc(sizeof(InBits) * (nRefSizeInWordSize + 1));

	if (refGenome == NULL)
		MEMORY_ALLOCATE_ERROR;

	char strReads[wordSize + 1];
	for (SIZE_T i = 0; i < nRefSizeInWordSize - 1; i++) {
		memcpy(&strReads, &(strRef[i * wordSize]), wordSize);
		strReads[wordSize] = 0;
		EncodeRead(strReads, &((*refGenome)[i]), wordSize);
	}
	int codesize = (nRefSizeInWordSize - 1) * wordSize;
	int remSize = nRefSize - codesize;
	memcpy(strReads, &(strRef[codesize]), (SIZE_T) remSize);
	strReads[remSize] = 0;
	EncodeRead(strReads, &((*refGenome)[nRefSizeInWordSize - 1]), remSize);
}

void GetReference(InBits ** refGenome, SIZE_T * nRefSize, const Option & opt) {
	LOG_INFO
	char * strRef;
	SIZE_T refLen = ReadWholeFile(opt.refFile, &strRef);

	*nRefSize = RemoveNonACGTNBase(strRef, refLen);
	RefEncodeToBits(refGenome, *nRefSize, strRef);
	free(strRef);
}
