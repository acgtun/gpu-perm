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
	GetParameter(argc, argv, &opt);

	/* (2) Read Reference Index */
	CReference refGenome;
	TIME_INFO(GetReference(&refGenome, opt), "encode reference");
	test("1.txt");
	/* (3) Build Reference Index */
	CHashTable d_hashTable;
	CReference d_refGenome;
	TIME_INFO(MakeHashTable(&refGenome, &d_refGenome, &d_hashTable), "build hash table");
	test("2.txt");
	/* (4) Build Reference Index */
	TIME_INFO(Matching(opt, &d_refGenome, &d_hashTable), "mapping");

	/* (5) free memory*/
	cudaFree(d_refGenome.refInBits);
	cudaFree(d_hashTable.counter);
	cudaFree(d_hashTable.index);
	return 0;
}
