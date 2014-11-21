#include "sdk.h"
#include "iofile.h"
#include "option.h"
#include "hash.h"

int main(int argc, char* argv[]) {
	Option opt;
	GetParameter(argc, argv, opt);
	
	CReference refGenome;
	CHashTable hashTable;

	TIME_INFO(BuildIndex(opt, &refGenome, &hashTable), "Build index");

	/* release memory*/
	free(refGenome.refSeq);
	free(hashTable.counter);
	free(hashTable.index);

	return 0;
}
