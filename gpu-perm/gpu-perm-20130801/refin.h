#pragma once
#ifndef REFIN_H_
#define REFIN_H_

#include "stdafx.h"
#include "option.h"
#include "iofile.h"
#include "bitscode.h"

//typedef struct {
//	char * strRef;
//	InBits * refInBits;
//	//SIZE_T nRefSize;
//	//SIZE_T * nRefSizeInWordSize;
//} RefGenome;

void GetReference(InBits ** refGenome, SIZE_T * nRefSize, const Option & opt);

#endif /* REFIN_H_ */
