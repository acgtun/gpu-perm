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
	printf("------haifengc------\n");
	Option opt;
	RefGenome refGenome;
	HashTable hashTable;
	GetParameter(argc, argv, &opt);

	/* (2) Build Reference Index */
	TIME_INFO(GetReference(&refGenome, &opt), "encode reference");
	TIME_INFO(MakeHashTable(&refGenome, &hashTable), "build hash table");
	TIME_INFO(Matching(&opt, &refGenome, &hashTable), "mapping");

	/* (3) mapping */
	//return 0;
}
