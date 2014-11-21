#ifndef ALLOCATE_H_
#define ALLOCATE_H_

#include "stdafx.h"
#include "bitscode.h"
#include "hash.h"
#include "refin.h"
#include "iofile.h"
#include "match.h"

void Matching(const Option & opt, const InBits * refGenome, const SIZE_T nRefSize, const SIZE_T * hashCounter, const SIZE_T * hashIndex);

#endif /* ALLOCATE_H_ */
