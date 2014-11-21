#include "seed.h"

/*
 * seed (111*1**)(111*1**)(111*1**)(111*1**)
 * hash (111*1**)(111*1**)(111*1**)(1)
 */

void GetF2Seed(WORD_SIZE encode, const SIZE_T & len, WORD_SIZE & seed) {
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
	//WORD_SIZE code = 0x972E5CB972E5CB97;
	//cout << "F2SEEDPATTERN = " << F2SEEDPATTERN << endl;
	//cout << strlen(F2SEEDPATTERN) << endl;
	seed = 0;
	//int cnt = 0;
	//cout << "len = " << len << endl;
	for (SIZE_T i = 0; i < len; i++) {
		if (F2SEEDPATTERN[i] == '1') {
			//cnt++;
			seed <<= 1;
			seed += encode & 0x01;
		}
		encode >>= 1;
	}
	//cout << "seed = " << seed << endl;
	//cout << "cnt = " << cnt << endl;
}

SIZE_T GetHashValue(const InBits & r) {
	WORD_SIZE hashValueU = 0, hashValueL = 0;
	//cout << "xfasdf" << endl;
	GetF2Seed(r.ub, LEN_FOR_13_CARE, hashValueU);
	GetF2Seed(r.lb, LEN_FOR_13_CARE, hashValueL);
	//cout << "gegegege" << endl;
	//cout << hashValueU << endl;
	//cout << hashValueL << endl;
	//cout << "------------------------" << endl;
	return (SIZE_T) ((hashValueU << BITS_FOR_HASHING) + hashValueL);
}

SIZE_T GetHashValue(const char * strVal, const SIZE_T & nStartPos) {
	WORD_SIZE hashValueU = 0, hashValueL = 0;
	for (SIZE_T i = 0; i < LEN_FOR_13_CARE; i++) {
		if (F2SEEDPATTERN[i] == '1') {
			hashValueU <<= 1;
			hashValueL <<= 1;
			if (strVal[nStartPos + i] == 'A') {
				//hashValueU += 0;
				//hashValueL += 0;
			} else if (strVal[nStartPos + i] == 'C') {
				//hashValueU += 0;
				hashValueL += 1;
			} else if (strVal[nStartPos + i] == 'G') {
				hashValueU += 1;
				//hashValueL += 0;
			} else { // (strVal[nStartPos + i] == 'T')
				hashValueU += 1;
				hashValueL += 1;
			}
		}
	}
	return (SIZE_T) ((hashValueU << BITS_FOR_HASHING) + hashValueL);
}

SIZE_T GetKmer(const CReference * refGenome, const SIZE_T & nRefStart, SIZE_T kmerLen, InBits * r) {
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
