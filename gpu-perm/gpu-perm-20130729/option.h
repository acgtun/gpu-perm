#pragma once
#ifndef PARAMETER_H_
#define PARAMETER_H_

#include "stdafx.h"
#include "bitscode.h"

typedef struct {
	//mapping parameter
	SIZE_T nMaxMismatch;
} MatchOpt;

typedef struct {
	//input and out parameter
	char readsFile[MAX_FILEPATH_LEN];
	char refFile[MAX_FILEPATH_LEN];
	int bSaveIndex; // Default is false

	char outputFile[MAX_FILEPATH_LEN];

	MatchOpt matchOpt;
} Option;

void PrintSynopsis();
void GetParameter(int argc, const char* argv[], Option * opt);
int GetIntVal(int argc, const char** argv, const char * str, SIZE_T * nVal);
void GetStrVal(int argc, const char** argv, const char* str, char * strVal);

#endif /* PARAMETER_H_ */
