#include "seed.h"

/*
 * seed (111*1**)(111*1**)(111*1**)(111*1**)
 * hash (111*1**)(111*1**)(111*1**)(1)
 */

SIZE_T GetKmer(const CReference * refGenome, SIZE_T nRefStart, SIZE_T kmerLen, InBits * r) {
	/* This function get kmerLen characters from the genome, start from position nRefstart*/
	if (nRefStart + kmerLen > refGenome->nRefSize) {
		kmerLen = refGenome->nRefSize - nRefStart;
	}

	r->ub = 0;
	r->lb = 0;

	SIZE_T indexInWords = nRefStart / wordSize;
	SIZE_T bitsShift = nRefStart % wordSize;
	/* the bitsShift bits are delete from the refInBits[indexInWords]
	 * get (WORD_SIZE - bitsShift) bits */
	r->ub = refGenome->refInBits[indexInWords].ub >> bitsShift;
	r->lb = refGenome->refInBits[indexInWords].lb >> bitsShift;

	/* kmer in two WORDSIZE, here kmerLen should less than WORD_SIZE */
	if (bitsShift != 0) {
		/* delete the high  (wordSize - bitsShift) bits, and get (bitsShift) bits */
		r->ub |= (refGenome->refInBits[indexInWords + 1].ub << (wordSize - bitsShift));
		r->lb |= (refGenome->refInBits[indexInWords + 1].lb << (wordSize - bitsShift));
	}

	SIZE_T elimatedBits = wordSize - kmerLen;
	r->ub <<= elimatedBits;
	r->lb <<= elimatedBits;
	r->ub >>= elimatedBits;
	r->lb >>= elimatedBits;

	return kmerLen;
}

int GetF2Seed(WORD_SIZE encode, SIZE_T len, WORD_SIZE * seed) {
	/* This function is based on the seed pattern (111*1**)(111*1**)(111*1**)1.
	 * To save the number of CPU instructions, the value encoded bits that doesn't follow the order of reads
	 * It still gets the first 3 digit and the fifth digit.. and so on to get total 13 digits
	 * The first 22 base are used as hash key (111*1**)(111*1**)(111*1**)1
	 * 111*1**111*1**111*1**1
	 * 1234567890123456789012
	 *
	 * 1001011100101110010111001011100101110010111001011100101110010111
	 * 0xE9D3A74E9D3A74E9
	 * */
	int nBits = 0;
	WORD_SIZE code = 0x972E5CB972E5CB97;
	*seed = 0;
	for (SIZE_T i = 1; i <= len; i++) {
		if (code & 0x01) {
			nBits++;
			*seed <<= 1;
			*seed += encode & 0x01;
		}
		encode >>= 1;
		code >>= 1;
	}
	return nBits;
}

int GetF2SeedForBits(InBits r, SIZE_T len, InBits * ret) {
	int s = GetF2Seed(r.ub, len, &ret->ub);
	GetF2Seed(r.lb, len, &ret->lb);
	return s;
}

SIZE_T GetHashValue(InBits r) {
	WORD_SIZE hashValueU, hashValueL;
	GetF2Seed(r.ub, 22, &hashValueU);
	GetF2Seed(r.lb, 22, &hashValueL);
	WORD_SIZE hashValue = (hashValueU << BITS_FOR_HASHING) + hashValueL;
	return ((SIZE_T) hashValue);
}
