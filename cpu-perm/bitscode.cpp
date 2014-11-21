#include "bitscode.h"

/*
 * uiReadLength must < WORD_SIZE
 * WORD_SIZE is 32 bp in 32 big machine and 64 bp in 64 bit machine
 * Each base is encoded into 2 bits: A -> 00, C->01, G->10 and T->11.
 * These two digits are located in two word, for bits operation.
 * The first nucleotide is encoded as the last digit.
 */

void EncodeRead(const char * strRead, InBits * readsInBits, int len) {
  /* This function transfers the ACGT string to binary format */
  readsInBits->ub = 0;
  readsInBits->lb = 0;

  // A 00
  // C 01
  // G 10
  // T 11

  for (int i = len - 1; i >= 0; i--) {
    if (strRead[i] == 'A' || strRead[i] == 'a') {

    } else if (strRead[i] == 'C' || strRead[i] == 'c') {
      readsInBits->lb++;
    } else if (strRead[i] == 'G' || strRead[i] == 'g') {
      readsInBits->ub++;
    } else if (strRead[i] == 'T' || strRead[i] == 't') {
      readsInBits->ub++;
      readsInBits->lb++;
    } else {
      //printf("Not A, C, G, T\n");
      //anything else as A
    }
    if (i != 0) {
      readsInBits->ub <<= 1;  //left shift 1
      readsInBits->lb <<= 1;
    }
  }
}

void DecodeRead(char * strReads, int readLen, const InBits * readsInBits) {
  /* This function transfers binary format to ACGT string */
  WORD_SIZE UpperBits = readsInBits->ub;
  WORD_SIZE LowerBits = readsInBits->lb;
  int strReadsl = 0;
  for (int i = 0; i < readLen; i++) {
    WORD_SIZE c = (UpperBits & 0x01) << 1 | (LowerBits & 0x01);
    switch (c) {
      case 0x00:
        strReads[strReadsl++] = 'A';
        break;
      case 0x01:
        strReads[strReadsl++] = 'C';
        break;
      case 0x02:
        strReads[strReadsl++] = 'G';
        break;
      case 0x03:
        strReads[strReadsl++] = 'T';
        break;
      default:
        strReads[strReadsl++] = 'N';
        break;
    }
    LowerBits >>= 1;
    UpperBits >>= 1;
  }
  strReads[strReadsl] = 0;
}

void DecodeReadReverse(char * strRead, int readLen,
                       const InBits * readsInBits) {
  DecodeRead(strRead, readLen, readsInBits);
  Swap(strRead, readLen);
}

SIZE_T bitsSetCount(WORD_SIZE bits) {
  /* magic function to calculate how many ones are there */
  SIZE_T c;  // c accumulates the total bits set in v
  for (c = 0; bits; c++) {
    bits &= bits - 1;  // clear the least significant bit set
  }
  return c;
}

char getNT(const int & nt) {
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

SIZE_T bitsStrNCompare(InBits r1, InBits r2, SIZE_T N) {
  /* compare only the last N bases (bits) */
  WORD_SIZE bits = (r1.ub ^ r2.ub) | (r1.lb ^ r2.lb);
  bits <<= (wordSize - N);
  return ((SIZE_T) bitsSetCount(bits));
}

void Swap(char * strVal, int len) {
  char chr;
  for (int i = 0; i < len / 2; i++) {
    chr = strVal[i];
    strVal[i] = strVal[len - i - 1];
    strVal[len - i - 1] = chr;
  }
}

char complimentBase(char nt) {
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
