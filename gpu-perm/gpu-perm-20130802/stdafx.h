#pragma once
#ifndef STDAFX_H_
#define STDAFX_H_

#include <time.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "cuda.h"
#include "cuda_runtime.h"

typedef unsigned int SIZE_T;

#define MAX_FILEPATH_LEN 256
#define MAX_LINE_LEN 2048
#define MIN_READ_LEN 13
#define MAX_READ_LEN 64

typedef unsigned long long WORD_SIZE; //no matter 32-bit or 64-bit, make sure they are coincident.
#define wordSize (8 * sizeof(WORD_SIZE))

#define NO_OF_BUCKET 0x4000000
#define BITS_FOR_HASHING 13  /* use the first 13 characters as the hash value */
#define NUMBER_OF_SHIFT 6 /* one less than the length of pattern */

#define BLOCKS 512 
#define THREADS 512

#define MAX_MAPPING_READS ((BLOCKS) * (THREADS)) /* 256 * 256 */
#define MAX_NUMBER_OF_RESULT 200

#define MAX_CHAR_INPUT  ((wordSize) * 3932160) //300Mb
#define MAX_CHAR_OUTPUT  ((wordSize) * 3932160) //300Mb
#define TYPE_NEW_LINE 0
#define TYPE_SPACE 1
#define TYPE_NOTHING 2

#define  QuickSort_MAX_LEVELS  300

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#define printf(f, ...) ((void)(f, __VA_ARGS__),0)
#endif

static void HandleError(cudaError_t err, const char *file, int line) {
	if (err != cudaSuccess) {
		printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
		exit(EXIT_FAILURE);
	}
}

static void MemoryAllocateCheck(void * pointer, const char * file, int line) {
	if (pointer == NULL) {
		printf("Memory allocate error in %s at line %d\n", file, line);
		exit(EXIT_FAILURE);
	}
}

static void FileOpenCheck(FILE * pfile, const char * file, int line) {
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

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define MEMORY_ALLOCATE_CHECK( pointer )  (MemoryAllocateCheck(pointer, __FILE__, __LINE__))
#define FILE_OPEN_CHECK(pfile) (FileOpenCheck( pfile, __FILE__, __LINE__))
#define LOG_INFO printf("--- %s:%s:%d\n",  __FILE__, __func__, __LINE__)

#define TIME_INFO(func, msg) { \
	clock_t start_t, end_t; \
	start_t = clock(); \
	func; \
	end_t = clock(); \
	printf("%s takes %.3lf seconds.\n", msg, (double) ((end_t - start_t) )); \
}
//	printf("%s takes %.3lf seconds.\n", msg, (double) ((end_t - start_t) / CLOCKS_PER_SEC));

typedef struct {
	SIZE_T nMismatch[MAX_NUMBER_OF_RESULT]; //it's very dangerous.
	SIZE_T nStartPos[MAX_NUMBER_OF_RESULT];
	SIZE_T nRet;
} CResult;

inline int isACGT(char nt) {
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
