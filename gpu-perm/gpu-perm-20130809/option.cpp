#include "option.h"

void PrintSynopsis() {
	printf("The input command is incorrect.\nFor more info, please check: http://code.google.com/p/perm/\n");
}

int GetIntVal(int argc, const char** argv, const char * str, SIZE_T * nVal) {
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], str) == 0 && i + 1 < argc - 1) {
			if (argv[i + 1][0] != '-') {
				*nVal = (SIZE_T) atoi(argv[i + 1]);
				return 1;
			}
		}
	}
	return 0;
}

int GetStrVal(int argc, const char** argv, const char* str, char * strVal) {
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], str) == 0 && i + 1 < argc - 1) {
			if (argv[i + 1][0] != '-') {
				strcpy(strVal, argv[i + 1]);
				return 1;
			}
		}
	}
	return 0;
}

int ChkStrExist(int argc, const char** argv, const char* str) {
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], str) == 0) {
			return 1;
		}
	}
	return 0;
}

void GetNameAfterDot(char * strFile, char * fileName) {
	int size = strlen(strFile);
	int i = 0, j = 0;
	for (int i = size - 1; i >= 0; i--) {
		if (strFile[i] == '.')
			break;
	}
	if (i <= 0) {
		strcpy(fileName, strFile);
	}
	int k = 0;
	for (j = i; j < size; j++) {
		fileName[k] = strFile[j];
		k++;
	}
	fileName[k] = 0;
}

void GetNameBeforeDot(char * strFile, char * fileName) {
	int size = strlen(strFile);
	int i = 0, j = 0;
	for (i = size - 1; i >= 0; i--) {
		if (strFile[i] == '.')
			break;
	}
	if (i <= 0) {
		strcpy(fileName, strFile);
	}
	for (j = 0; j < i; j++) {
		fileName[j] = strFile[j];
	}
	fileName[j] = 0;
}

int CheckIndexExsit(char * strFile) {
	//.index
	char fileName[MAX_FILEPATH_LEN] = {0};
	GetNameAfterDot(strFile, fileName);
	if (strcmp(fileName, ".index") == 0)
		return 1;
	return 0;
}

void GetParameter(int argc, const char* argv[], Option * opt) {
	LOG_INFO;
	if (argc <= 2 || argv == NULL)
		PrintSynopsis();
	else {
		strcpy(opt->refFile, argv[1]);
		opt->bIndexExist = 0;
		if (CheckIndexExsit(opt->refFile)) {
			opt->bIndexExist = 1;
		}
		strcpy(opt->readsFile, argv[2]);
		if (GetIntVal(argc, argv, "-v", &opt->mapOpt.nMaxMismatch) == 0)
			opt->mapOpt.nMaxMismatch = 2;
		if (GetStrVal(argc, argv, "-o", opt->outputFile) == 0) {
			char fileName[MAX_FILEPATH_LEN] = {0};
			GetNameBeforeDot(opt->readsFile, fileName);
			strcpy(opt->outputFile, fileName);
			strcat(opt->outputFile, ".sam");
		}
		opt->bSaveIndex = 0;
		if (ChkStrExist(argc, argv, "-s")) {
			opt->bSaveIndex = 1;
		}
		if (GetStrVal(argc, argv, "-s", opt->indexFile) == 0 && opt->bSaveIndex == 1) {
			char fileName[MAX_FILEPATH_LEN] = {0};
			GetNameBeforeDot(opt->refFile, fileName);
			strcpy(opt->indexFile, fileName);
			strcat(opt->indexFile, ".index");
		}
	}
}
