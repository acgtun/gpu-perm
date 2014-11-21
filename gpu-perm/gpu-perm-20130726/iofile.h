#pragma once
#ifndef IOFILE_H_
#define IOFILE_H_

#include "stdafx.h"

SIZE_T ReadWholeFile(const char * fileName, char ** strRet);
int GetLineFromString(const char * strVal, char * strRet);

void WriteWholeFile(const char * fileName, char * strRet, SIZE_T len);

#endif /* IOFILE_H_ */
