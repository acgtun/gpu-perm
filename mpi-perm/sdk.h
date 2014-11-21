#pragma once
#ifndef STDAFX_H_
#define STDAFX_H_

#include <time.h>
#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include <map>
#include <set>
#include <limits>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#define debug

using namespace std;

const double GB = 1024 * 1024 * 1024;

#define HASHLEN 13     /* use the first 13 characters as the hash value */
#define NO_OF_BUCKET 0x4000000  /* number of bucket in the hash table, 4^BITS_FOR_HASHING */
#define MAX_LINE_LEN 10024

inline void MemoryAllocateCheck(void * pointer, const char * file, int line) {
  if (pointer == NULL) {
    printf("Memory allocate error in %s at line %d\n", file, line);
    exit(EXIT_FAILURE);
  }
}

inline void FileOpenCheck(FILE * pfile, const char * file, int line) {
  if (pfile == NULL) {
    printf("File open error in %s at line %d\n", file, line);
    exit(EXIT_FAILURE);
  }
}

#define HANDLE_ERROR(err) (HandleError( err, __FILE__, __LINE__ ))
#define LOG_INFO printf("--- %s:%s:%d\n",  __FILE__, __func__, __LINE__)
#define LOG_INFOrank(rank)  printf("[rank %d] %s:%s:%d\n", rank,  __FILE__, __func__, __LINE__)
#define LOG_INFOranki(rank, i)  printf("[rank %d] i = %d: %s:%s:%d\n", rank, i,  __FILE__, __func__, __LINE__)
#define FILE_OPEN_CHECK(pfile) (FileOpenCheck( pfile, __FILE__, __LINE__))
#define CHECK_READ_LEN(len, nReadsNum) (checkReadLen(len, nReadsNum, __FILE__, __LINE__))
#define MEMORY_ALLOCATE_CHECK(pointer)  (MemoryAllocateCheck(pointer, __FILE__, __LINE__))
#define LOG_INFO_CPP cout << "-----" << __FILE__ << " " << __func__ << " " << __LINE__ << endl

#define FREAD_CHECK(func, size) { \
	uint32_t s = func; \
	if(s != size) { \
		printf("read file error. --- %s:%s:%d\n", __FILE__, __func__, __LINE__); \
		exit(EXIT_FAILURE); \
	} \
}

#define TIME_INFO(func, msg) { \
	clock_t start_t, end_t; \
	start_t = clock(); \
	func; \
	end_t = clock(); \
	printf("--INFO-- %s takes %.3lf seconds.\n", msg, (double) ((end_t - start_t) / CLOCKS_PER_SEC )); \
}

inline int isACGT(const char & nt) {
  switch (nt) {
    case 'a':
    case 'c':
    case 'g':
    case 't':
    case 'A':
    case 'C':
    case 'G':
    case 'T':
      return 1;
    default:
      return 0;
  }
}

inline void INFO(const string & msg) {
  cout << "--INFO-- " << msg << "." << endl;
}
/*
 inline void INFO(const int & rank) {
 cout << "[rank " << rank << "]";
 printf("--- %s:%s:%d\n",  __FILE__, __func__, __LINE__);
 }
 */
inline void INFO(const string & msg, const string & val) {
  cout << "--INFO-- " << msg << " " << val << "." << endl;
}
inline void INFO(const string & msg, const uint32_t &val) {
  cout << "--INFO-- " << msg << " " << val << "." << endl;
}
inline void INFO(const string & msg, const uint32_t & val1,
                 const string & val2) {
  cout << "--INFO-- " << msg << " " << val1 << " " << val2 << "." << endl;
}

inline uint32_t GetHashValue(const char * strSeed) {
  uint32_t hashValue = 0;
  for (uint32_t i = 0; i < HASHLEN; i++) {
    hashValue <<= 2;
    switch (strSeed[i]) {
      case 'A':
        hashValue += 0;
        break;
      case 'C':
        hashValue += 1;
        break;
      case 'G':
        hashValue += 2;
        break;
      case 'T':
        hashValue += 3;
        break;
    }
  }
  return hashValue;
}

inline char getNT(const int & nt) {
  switch (nt) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
  }
  return 'A';
}

#endif /* STDAFX_H_ */
