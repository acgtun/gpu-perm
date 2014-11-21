#pragma once
#ifndef PARAMETER_H_
#define PARAMETER_H_

#include "stdafx.h"
#include "bitscode.h"

typedef struct {
	/* mapping parameter */
	SIZE_T nMaxMismatch;
} MapOpt;

typedef struct {
	//input and out parameter
	char readsFile[MAX_FILEPATH_LEN];
	char refFile[MAX_FILEPATH_LEN];

	char outputFile[MAX_FILEPATH_LEN];
	char indexFile[MAX_FILEPATH_LEN];
	int bSaveIndex; // Default is false
	int bIndexExist;

	MapOpt mapOpt;

	/////////////////////////////////
	set<SIZE_T> setConsectiveN;
} Option;

void PrintSynopsis();
void GetParameter(int argc, const char* argv[], Option * opt);

#endif /* PARAMETER_H_ */
