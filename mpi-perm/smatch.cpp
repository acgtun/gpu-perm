#include "sdk.h"
#include "iofile.h"
#include "option.h"
#include "hash.h"
#include "match.h"

int main(int argc, char* argv[]) {

  /* (1) Get parameters */
  Option opt;
  GetParameter(argc, argv, opt);

  CReference refGenome;
  CHashTable hashTable;

  vector <string> vReadsSeq;

  /* (1) Read Reference and Hash Table from index */
  TIME_INFO(ReadIndexAndRef(opt, &refGenome, &hashTable),
            "Read reference and hash table");

  /* (2) Read the reads sequences */
  ReadReads(opt, vReadsSeq);

  /* (3) Serial Matching */
  TIME_INFO(Serial_Matching(opt, &refGenome, &hashTable, vReadsSeq),
            "Serial matching");

  /* (4) release memory*/
  free(refGenome.refSeq);
  free(hashTable.counter);
  free(hashTable.index);

  return 0;
}
