#pragma once
#ifndef STDAFX_H_
#define STDAFX_H_

#include <time.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <map>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

//#define CUDA

#ifdef CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#define BLOCKS 512
#define THREADS 512

static void HandleError(cudaError_t err, const char *file, int line) {
	if (err != cudaSuccess) {
		printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
		//exit(EXIT_FAILURE);
	} else {
		printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
	}
}
#endif

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#define printf(f, ...) ((void)(f, __VA_ARGS__),0)
#endif

typedef unsigned int SIZE_T;

#define MAX_LINE_LEN 1024
#define MIN_READ_LEN 30
#define MAX_READ_LEN 193

typedef unsigned long long WORD_SIZE; /* no matter 32-bit or 64-bit, make sure they are identical. */
#define wordSize (8 * sizeof(WORD_SIZE))

#define NO_OF_BUCKET 0x4000000  /* number of bucket in the hash table, 4^BITS_FOR_HASHING */
#define BITS_FOR_HASHING 13     /* use the first 13 characters as the hash value */
#define LEN_FOR_13_CARE 22      /* the first 22 nucleotide contains 13 care ones */
#define NUMBER_OF_SHIFT 7		/* one less than the length of pattern */
#define NUMBER_OF_SPACE 1       /* if NUMBER_OF_SPACE = 5, every 5 position index one */
#define TOTAL_SHIFT (NUMBER_OF_SHIFT * NUMBER_OF_SPACE)
#define TOTAL_SHIFT_rev (2 * TOTAL_SHIFT)

#define MAX_MAPPING_READS 1000000 /* 256 * 256 */

#define MAX_CHAR_INPUT  ((wordSize) * 3932160) //300Mb
#define MAX_CHAR_OUTPUT  ((wordSize) * 3932160) //300Mb
#define TYPE_NEW_LINE 0
#define TYPE_SPACE 1
#define TYPE_NOTHING 2

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

inline void printWORD_SIZE(WORD_SIZE word, SIZE_T len) {
	/* print the WORD_SIZE as binary format */
	for (SIZE_T i = 0; i < len; i++) {
		if ((word & 0x01) == 1) {
			printf("1");
		} else {
			printf("0");
		}
		word >>= 1;
	}
	printf("\n");
}

inline void checkReadLen(const SIZE_T & len, const SIZE_T & nReadsNum, const char * file, int line) {
	if (len > MAX_READ_LEN || len < MIN_READ_LEN) {
		cout << "The length of read " << nReadsNum + 1 << " is not between " << MIN_READ_LEN << " and " << MAX_READ_LEN << ". It will be ignored."
				<< "--- in " << file << " at line " << line << endl;
	}
}

#define INFO(msg, file) cout << msg << " " << file << endl
#define HANDLE_ERROR(err) (HandleError( err, __FILE__, __LINE__ ))
#define LOG_INFO printf("--- %s:%s:%d\n",  __FILE__, __func__, __LINE__)
#define FILE_OPEN_CHECK(pfile) (FileOpenCheck( pfile, __FILE__, __LINE__))
#define CHECK_READ_LEN(len, nReadsNum) (checkReadLen(len, nReadsNum, __FILE__, __LINE__))
#define MEMORY_ALLOCATE_CHECK(pointer)  (MemoryAllocateCheck(pointer, __FILE__, __LINE__))
#define LOG_INFO_CPP cout << "-----" << __FILE__ << " " << __func__ << " " << __LINE__ << endl

#define FREAD_CHECK(func, size) { \
	SIZE_T s = func; \
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
	printf("%s takes %.3lf seconds.\n", msg, (double) ((end_t - start_t) / CLOCKS_PER_SEC )); \
}

typedef struct {
	SIZE_T nMismatch; //it's very dangerous.
	SIZE_T nStartPos;
	char org_rev;
} CResult;

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

#endif /* STDAFX_H_ */
