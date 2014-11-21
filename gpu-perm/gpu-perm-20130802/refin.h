#pragma once
#ifndef REFIN_H_
#define REFIN_H_

#include "stdafx.h"
#include "option.h"
#include "iofile.h"
#include "bitscode.h"

typedef struct {
	InBits * refInBits;
	SIZE_T nRefSize;
	SIZE_T nRefSizeInWordSize;
} CReference;

void GetReference(CReference * refGenome, const Option & opt);

#endif /* REFIN_H_ */
