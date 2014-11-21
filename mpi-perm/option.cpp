#include "option.h"

void PrintSynopsis() {
  printf(
      "The input command is incorrect. please reference the following command. \n "
      "./MPImapping [reference file or index file] [reads file]\n");
}

int GetStrVal(int argc, char** argv, const char* str, string & strVal) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], str) == 0 && i + 1 <= argc - 1) {
      if (argv[i + 1][0] != '-') {
        strVal = argv[i + 1];
        return 1;
      }
    }
  }
  return 0;
}

int ChkStrExist(int argc, char** argv, const char* str) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], str) == 0) {
      return 1;
    }
  }
  return 0;
}

void GetNameBeforeDot(const string & strFile, string & fileName) {
  int pos = strFile.find_last_of('.');
  fileName = strFile.substr(0, pos);
}

void GetReadLength(Option & opt) {
  ifstream fin(opt.readsFile.c_str());
  if (!fin.good()) {
    printf("--ERROR INFO-- reads file open error. %s\n", opt.readsFile.c_str());
    exit (EXIT_FAILURE);
  }
  opt.nNumOfreads = 0;
  string strLine;
  while (getline(fin, strLine)) {
    if (strLine[0] == '>') {
      opt.nNumOfreads++;
      continue;
    }
    opt.readLen = strLine.size();
    break;
  }

  while (getline(fin, strLine)) {
    if (strLine[0] == '>') {
      opt.nNumOfreads++;
      continue;
    }
    if (opt.readLen != strLine.size()) {
      printf("--ERROR INFO-- The lengths of the reads are not identical.\n");
      exit (EXIT_FAILURE);

    }
  }
  fin.close();
  INFO("There are", opt.nNumOfreads, "reads in the reads file");
  INFO("The length of the read is", opt.readLen);
}

bool isBuildIndex(const string strVal, const string strFind) {
  if (strVal.find(strFind) != string::npos)
    return true;
  return false;
}

void GetParameter(int argc, char* argv[], Option & opt) {
  cout << "--INFO-- Input command:";
  for (int i = 0; i < argc; i++) {
    cout << " " << argv[i];
  }
  cout << endl;

  /* reference file name */
  opt.refFile = argv[1];
  cout << "--INFO-- The reference file is " << opt.refFile << endl;

  /* determine the index file exist of not */
  opt.bIndexExist = 0;
  string indexsf = ".sdindex";
  if (opt.refFile.size() > indexsf.size()
      && opt.refFile.substr(opt.refFile.size() - indexsf.size()) == indexsf) {
    opt.bIndexExist = 1;
  }
  if (!isBuildIndex(argv[0], "buildindex")) {
    opt.readsFile = argv[2];
    cout << "--INFO-- The reads file is " << opt.readsFile << "." << endl;
  }
  /* setting the output file name */
  if (GetStrVal(argc, argv, "-o", opt.outputFile) == 0) {
    string fileName;
    GetNameBeforeDot(opt.readsFile, fileName);
    opt.outputFile = fileName;
    opt.outputFile += ".out";
  }

  /* determine whether save the index file or not*/
  opt.bSaveIndex = 0;
  if (ChkStrExist(argc, argv, "-s")) {
    opt.bSaveIndex = 1;
  }
  if ((GetStrVal(argc, argv, "-s", opt.indexFile) == 0 && opt.bSaveIndex == 1)
      || isBuildIndex(argv[0], "buildindex")) {
    string fileName;
    opt.bSaveIndex = 1;
    GetNameBeforeDot(opt.refFile, fileName);
    opt.indexFile = fileName;
    opt.indexFile += indexsf;
  }
  if (opt.bSaveIndex != 1)
    GetReadLength(opt);
}

