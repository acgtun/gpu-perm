#ifndef MATCHING_H_
#define MATCHING_H_

#include "sdk.h"
#include "hash.h"
#include "option.h"

inline void NoRegion(pair<uint32_t, uint32_t> & ret) {
  ret.first = 1;
  ret.second = 0;
}

void Serial_Matching(const Option & opt, const CReference * refGenome,
                     const CHashTable * hashTable,
                     const vector<string> & vReadsSeq);
void ReadReads(const Option & opt, vector<string> & vReadsSeq);
void MappingOneRead(const CReference * refGenome, const CHashTable * hashTable,
                    const char * querySeq, pair<uint32_t, uint32_t> & ret,
                    const uint32_t & startPos);

#endif /* MATCHING_H_ */
