#pragma once
#ifndef PARAMETER_H_
#define PARAMETER_H_

#include "sdk.h"

typedef struct {
  //input and out parameter

  string readsFile;
  uint32_t readLen;
  string refFile;

  uint32_t nNumOfreads;

  string outputFile;
  string indexFile;

  int bSaveIndex;  // Default is false
  int bIndexExist;

} Option;

void PrintSynopsis();
void GetParameter(int argc, char* argv[], Option & opt);

#endif /* PARAMETER_H_ */
