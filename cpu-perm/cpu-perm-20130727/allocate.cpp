#include "allocate.h"

void Matching(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable) {
	CPU_Matching(opt, refGenome, hashTable);
}
