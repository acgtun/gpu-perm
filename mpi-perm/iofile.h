#pragma once
#ifndef IOFILE_H_
#define IOFILE_H_

#include "sdk.h"

uint32_t ReadWholeFile(const string & fileName, char ** strRet);
uint32_t GetLineFromString(const char * strVal, char * strRet);

#endif /* IOFILE_H_ */
