#include "option.h"

void PrintSynopsis() {
  printf(
      "The input command is incorrect.\nFor more info, please check: http://code.google.com/p/perm/\n");
}

int GetIntVal(int argc, const char** argv, const char * str, SIZE_T & nVal) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], str) == 0 && i + 1 < argc - 1) {
      if (argv[i + 1][0] != '-') {
        nVal = (SIZE_T) atoi(argv[i + 1]);
        return 1;
      }
    }
  }
  return 0;
}

int GetStrVal(int argc, const char** argv, const char* str, string & strVal) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], str) == 0 && i + 1 < argc - 1) {
      if (argv[i + 1][0] != '-') {
        strVal = argv[i + 1];
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

//void GetNameAfterDot(char * strFile, char * fileName) {
//	int size = strlen(strFile);
//	int i = 0, j = 0;
//	for (int i = size - 1; i >= 0; i--) {
//		if (strFile[i] == '.')
//			break;
//	}
//	if (i <= 0) {
//		strcpy(fileName, strFile);
//	}
//	int k = 0;
//	for (j = i; j < size; j++) {
//		fileName[k] = strFile[j];
//		k++;
//	}
//	fileName[k] = 0;
//}

void GetNameBeforeDot(const string & strFile, string & fileName) {
  int pos = strFile.find_last_of('.');
  fileName = strFile.substr(0, pos);
}

void GetParameter(int argc, const char* argv[], Option & opt) {
  LOG_INFO;
  if (argc <= 2 || argv == NULL)
    PrintSynopsis();
  else {
    opt.refFile = argv[1];
    cout << "The reference file is " << opt.refFile << endl;
    cout << opt.refFile.length() << endl;
    cout << opt.refFile.size() << endl;
    opt.bIndexExist = 0;
    cout << opt.refFile.size() << endl;
    cout << opt.refFile.substr(opt.refFile.size() - 6) << endl;
    cout << "hah" << endl;
    if (opt.refFile.size() > 6
        && opt.refFile.substr(opt.refFile.size() - 6) == ".index") {
      opt.bIndexExist = 1;
    }
    opt.readsFile = argv[2];
    cout << "The reads file is " << opt.readsFile << endl;

    if (GetIntVal(argc, argv, "-v", opt.mapOpt.nMaxMismatch) == 0)
      opt.mapOpt.nMaxMismatch = 2;

    if (GetStrVal(argc, argv, "-o", opt.outputFile) == 0) {
      string fileName;
      GetNameBeforeDot(opt.readsFile, fileName);
      opt.outputFile = fileName;
      opt.outputFile += ".sam";
    }
    opt.bSaveIndex = 0;
    if (ChkStrExist(argc, argv, "-s")) {
      opt.bSaveIndex = 1;
    }
    if (GetStrVal(argc, argv, "-s", opt.indexFile) == 0
        && opt.bSaveIndex == 1) {
      string fileName;
      GetNameBeforeDot(opt.refFile, fileName);
      opt.indexFile = fileName;
      opt.indexFile += ".index";
    }
  }
}
