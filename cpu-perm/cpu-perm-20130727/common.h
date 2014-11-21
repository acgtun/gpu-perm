#pragma once
#ifndef COMMON_H_
#define COMMON_H_

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

inline void Swap(char * strVal, int len) {
	char chr;
	for (int i = 0; i < len / 2; i++) {
		chr = strVal[i];
		strVal[i] = strVal[len - i - 1];
		strVal[len - i - 1] = chr;
	}
}

inline char complimentBase(char nt) {
	switch (nt) {
	case 'a':
		return ('t');
	case 'c':
		return ('g');
	case 'g':
		return ('c');
	case 't':
		return ('a');
	case 'A':
		return ('T');
	case 'C':
		return ('G');
	case 'G':
		return ('C');
	case 'T':
		return ('A');
	default:
		return ('N');
	}
}

inline void reverseCompliment(char * strVal, int len) {
	Swap(strVal, len);
	for (int i = 0; i < len; i++) {
		strVal[i] = complimentBase(strVal[i]);
	}
}

#endif /* COMMON_H_ */
