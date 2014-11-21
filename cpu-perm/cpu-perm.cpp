/*
 * GPU-PerM
 * Haifeng Chen
 * University of Southern California
 * July 1, 2013
 *
 * Compiler C99
 */

/*
 * File name: all the letters are lower case
 * Function name: the first letter of every word are upper case, \
 *   the others are lower case
 * Variable name: the first letter of each word except \
 *   the first word are upper case, the others are lower case
 *
 * make it as simple as possible
 */

#include "option.h"
#include "bitscode.h"
#include "hash.h"
#include "refin.h"
#include "alloc_kernel.h"

int main(int argc, const char* argv[]) {
  /* (1) Get parameters */
  Option opt;
  GetParameter(argc, argv, opt);

  CReference refGenome;
  CHashTable hashTable;

  cout << "opt.bIndexExist = " << opt.bIndexExist << endl;
  if (opt.bIndexExist) {
    /* (2) Read Reference and Hash Table from index */
    TIME_INFO(ReadIndexAndRef(&refGenome, &hashTable, opt),
              "read reference and hash table");
  } else {
    /* (3) Read Reference */
    TIME_INFO(GetReference(&refGenome, opt), "encode reference");
    cout << "hf2" << endl;
    /* (4) Build Reference Index */
    TIME_INFO(MakeHashTable(&refGenome, &hashTable, opt), "build hash table");
  }

  /* (5) Build Reference Index */
  TIME_INFO(Matching(opt, &refGenome, &hashTable), "mapping");

  return 0;
}
