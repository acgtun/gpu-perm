#ifndef ALLOCATE_H_
#define ALLOCATE_H_

#include "stdafx.h"
#include "bitscode.h"
#include "hash.h"
#include "refin.h"
#include "iofile.h"
#include "match_kernel.h"

void Matching(const Option * opt, const RefGenome * refGenome,
		const HashTable * hashTable);

#endif /* ALLOCATE_H_ */
