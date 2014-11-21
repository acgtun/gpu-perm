#include "match.h"

int StrCMPproteinSeq(const char * querySeq, const uint32_t & index, const CReference * refGenome) {
	uint32_t i = 0, j = index, len = strlen(querySeq), size = refGenome->nRefSize;
	//cout << querySeq << " size = " << size <<  "len " <<  len << " index " << index << endl;
	while (i < len && j < size) {
		if (querySeq[i] > refGenome->refSeq[j])
			return 1;
		else if (querySeq[i] < refGenome->refSeq[j])
			return -1;
		i++;
		j++;
	}

	if (i == len)
		return 0;
	if (j == size)
		return 1;

	return 0;
}

uint32_t UpperBound(uint32_t low, uint32_t high, const char * querySeq, const CReference * refGenome,
		const CHashTable * hashTable) {
	uint32_t mid = 0;
	while (low < high) {
		mid = (low + high + 1) / 2;
		if (StrCMPproteinSeq(querySeq, hashTable->index[mid], refGenome) >= 0)
			low = mid;
		else
			high = mid - 1;
	}
	return low;
}

uint32_t LowerBound(uint32_t low, uint32_t high, const char * querySeq, const CReference * refGenome,
		const CHashTable * hashTable) {
	uint32_t mid = 0;
	while (low < high) {
		mid = (low + high) / 2;
		if (StrCMPproteinSeq(querySeq, hashTable->index[mid], refGenome) <= 0)
			high = mid;
		else
			low = mid + 1;
	}
	return low;
}

void MappingOneRead(const CReference * refGenome, const CHashTable * hashTable,
		const char * querySeq, pair<uint32_t, uint32_t> & ret, const uint32_t & startPos) {
	uint32_t hashValue = GetHashValue(querySeq) - startPos;
	//cout << GetHashValue(querySeq)<< " " << startPos << " " << querySeq <<  " "  << hashValue + 1 << " "  << hashTable->nSizeCounter << endl;
	
	uint32_t l = hashTable->counter[hashValue];
	if (hashTable->counter[hashValue + 1] == 0)
		return NoRegion(ret);
	uint32_t u = hashTable->counter[hashValue + 1] - 1;
	if (l > u)
		return NoRegion(ret);
	l = LowerBound(l, u, querySeq, refGenome, hashTable);
	u = UpperBound(l, u, querySeq, refGenome, hashTable);
	if (l > u)
		return NoRegion(ret);
	if (l == u) {
		if (StrCMPproteinSeq(querySeq, hashTable->index[l], refGenome) != 0)
			return NoRegion(ret);
	}

	ret.first = l;
	ret.second = u;
}

void ReadReads(const Option & opt, vector<string> & vReadsSeq) {
	INFO("Read read sequences from", opt.readsFile);
	char * strReads;
	uint32_t readsLen = ReadWholeFile(opt.readsFile, &strReads);
	char strRead[MAX_LINE_LEN];
	uint32_t readLen;
	for (uint32_t i = 0; i < readsLen; i++) {
		readLen = GetLineFromString(&strReads[i], strRead);
		i += readLen;
		if (strRead[0] == '>') {
			//sscanf(strRead, ">%s", readName);
			//vReadsName.push_back(readName);
			continue;
		} else {
			vReadsSeq.push_back(strRead);
		}
	}
	free(strReads);
}

void Serial_Matching(const Option & opt, const CReference * refGenome, const CHashTable * hashTable,
		const vector<string> & vReadsSeq) {
	ofstream fout((opt.outputFile + "_serial").c_str());
	pair<uint32_t, uint32_t> ret;
	for (uint32_t i = 0; i < vReadsSeq.size(); i++) {
		MappingOneRead(refGenome, hashTable, vReadsSeq[i].c_str(), ret, 0);

		fout << ">" << vReadsSeq[i] << endl;
		for (uint32_t j = ret.first; j <= ret.second; j++) {
			fout << hashTable->index[j] << " ";
		}
		fout << endl;
	}
	fout.close();
}
