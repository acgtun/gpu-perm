#include "iofile.h"

int GetLineFromString(const char * strVal, char * strRet) {
	int i;
	for (i = 0; strVal[i] != 0; i++) {
		if (0xa == strVal[i] || 0xd == strVal[i])
			break;
	}
	memcpy(strRet, strVal, (SIZE_T) i);
	strRet[i] = 0;
	return i;
}

SIZE_T ReadWholeFile(const char * fileName, char **strRet) {
	FILE * fin = fopen(fileName, "rb");
	if (fin == NULL) {
		FILE_OPEN_ERROR(fileName);
		fclose(fin);
		return 0;
	} else {
		fseek(fin, 0, SEEK_END);
		SIZE_T size = ftell(fin);
		*strRet = (char*) malloc(sizeof(char) * (size + 1));
		if (*strRet == NULL)
			MEMORY_ALLOCATE_ERROR;
		fseek(fin, 0, SEEK_SET);
		SIZE_T nsize = fread(*strRet, 1, size, fin);
		if (nsize != size) {
			perror("read file error!");
		}
		fclose(fin);
		return size;
	}
}

void WriteWholeFile(const char * fileName, char * strRet, SIZE_T len) {
	FILE * fout = fopen(fileName, "wb");
	if (fout == NULL) {
		FILE_OPEN_ERROR(fileName);
		fclose(fout);
	}
	fwrite(strRet, 1, len, fout);
	fclose(fout);
}

void OutPutResult(FILE * fout, ResultMatchedReads * result, SIZE_T nReads,
		int readID) {
	for (SIZE_T i = 0; i < nReads; i++) {
		fprintf(fout, "read %d:", readID + i);
		if (0 == result[i].nRet) {
			cout << "now result\n" << endl;
		}
		for (SIZE_T j = 0; j < result[i].nRet; j++) {
			fprintf(fout, " (%u, %u)", result[i].nStartPos[j],
					result[i].nMismatch[j]);
		}
		fprintf(fout, "\n");
	}
}
