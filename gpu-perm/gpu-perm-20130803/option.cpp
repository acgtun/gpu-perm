#include "option.h"

void PrintSynopsis() {
	printf("The input command is incorrect.\nFor more info, please check: http://code.google.com/p/perm/\n");
}

int GetIntVal(int argc, const char** argv, const char * str, SIZE_T * nVal) {
	for (int i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], str) == 0) {
			if (argv[i + 1][0] != '-') {
				*nVal = (SIZE_T) atoi(argv[i + 1]);
				return 1;
			}
		}
	}
	return 0;
}

void GetStrVal(int argc, const char** argv, const char* str, char * strVal) {
	for (int i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], str) == 0) {
			if (argv[i + 1][0] != '-') {
				strcpy(strVal, argv[i + 1]);
			}
			return;
		}
	}
}

int ChkStrExist(int argc, const char** argv, const char* str) {
	for (int i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], str) == 0) {
			return 1;
		}
	}
	return 0;
}

void GetParameter(int argc, const char* argv[], Option * opt) {
	LOG_INFO;
	if (argc <= 2 || argv == NULL)
		PrintSynopsis();
	else {
		strcpy(opt->refFile, argv[1]);
		strcpy(opt->readsFile, argv[2]);
		if (GetIntVal(argc, argv, "-v", &opt->mapOpt.nMaxMismatch) == 0)
			opt->mapOpt.nMaxMismatch = 2;
		//cout << "mismatch = " << opt->matchOpt.nMaxMismatch << endl;
		//system("pause");
	}
}
