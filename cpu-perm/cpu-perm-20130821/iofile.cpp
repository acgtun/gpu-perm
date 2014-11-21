#include "iofile.h"

int GetLineFromString(const char * strVal, char * strRet) {
	SIZE_T i;
	for (i = 0; strVal[i] != 0; i++) {
		//printf("%d %c\n", strVal[i], strVal[i]);
		if (0xA == strVal[i] || 0xD == strVal[i])
			break;
	}
	//cout << "i = " << i << endl;
	memcpy(strRet, strVal, i);
	strRet[i] = 0;
	//cout << "strRet = " << strRet << endl;
	return i;
}

SIZE_T ReadWholeFile(const string & fileName, char **strRet) {
	LOG_INFO;
	FILE * fin = fopen(fileName.c_str(), "rb");
	FILE_OPEN_CHECK(fin);
	LOG_INFO;
	fseek(fin, 0, SEEK_END);
	SIZE_T size = ftell(fin);
	MEMORY_ALLOCATE_CHECK(*strRet = (char* ) malloc(sizeof(char) * (size + 1)));
	fseek(fin, 0, SEEK_SET);
	FREAD_CHECK(fread(*strRet, 1, size, fin), size);
	fclose(fin);
	return size;
}

void WriteWholeFile(const string & fileName, char * strRet, SIZE_T len) {
	FILE * fout = fopen(fileName.c_str(), "wb");
	FILE_OPEN_CHECK(fout);
	fwrite(strRet, 1, len, fout);
	fclose(fout);
}

//void OutPutResult(FILE * fout, CResult * result, SIZE_T nReads, SIZE_T readID) {
//	for (SIZE_T i = 0; i < nReads; i++) {
//		fprintf(fout, "read %d:", readID + i);
//		if (0 == result[i].nRet) {
//			//cout << "no result\n" << endl;
//		}
//		for (SIZE_T j = 0; j < result[i].nRet; j++) {
//			fprintf(fout, " (%u, %u)", result[i].nStartPos[j], result[i].nMismatch[j]);
//		}
//		fprintf(fout, "\n");
//	}
//}
