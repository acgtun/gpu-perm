#include "sdk.h"
#include "hash.h"
#include "option.h"
#include "match.h"
#include <mpi.h>

void spiltHashTable(const CReference * refGenome, const CHashTable * hashTable, int numprocs,
		vector<uint32_t> & segStartID) {
	uint32_t avg = refGenome->nRefSize / (numprocs - 1);
	uint32_t j = avg;
	segStartID.push_back(0);
	for (uint32_t i = 0; i < hashTable->nSizeCounter; i++) {
		if (hashTable->counter[i] > j) {
			segStartID.push_back(i);
			j += avg;
		}
	}
	segStartID.push_back(hashTable->nSizeCounter);
#ifdef debug1
	cout << "segStartID size = " << segStartID.size() << " :";
	for (uint32_t i = 0; i < segStartID.size(); i++) {
		cout << " " << segStartID[i];
	}
	cout << endl;
#endif
}

uint32_t lower_bound(const vector<uint32_t> & vInt, const uint32_t & val) {
	uint32_t low = 0, high = vInt.size();
	while (low < high) {
		if (low + 1 == high) {
			if (val >= vInt[low + 1])
				return low + 1;
			return low;
		}
		uint32_t mid = (low + high) / 2;
		if (val < vInt[mid]) {
			high = mid - 1;
		} else
			low = mid;
	}
	return low;
}

void Parallel_Matching(int argc, char* argv[]) {

	int rank, namelen, numprocs;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	MPI_Status status;
	vector < uint32_t > segStartID;
	uint32_t nCounterStart, nCounterSize, nIndexSize;

	CReference refGenome;
	CHashTable hashTable;

	uint32_t nReads;
	uint32_t readLen;

	char outputFilename[100];
	uint32_t outputFilenamelen;

	vector < string > vReadsSeq;

	/* distribute index hash table to workers */
	if (rank == 0) {
		Option opt;
		GetParameter(argc, argv, opt);
		strcpy(outputFilename, opt.outputFile.c_str());
		outputFilenamelen = strlen(outputFilename);
		ReadIndexAndRef(opt, &refGenome, &hashTable);
		readLen = opt.readLen;
		ReadReads(opt, vReadsSeq);
		
		cout << "read len " << readLen << endl;
		//char rmcommand[100];
		//sprintf(rmcommand, "/bin/rm %s_parallel_rank*", outputFilename);
		//system(rmcommand);

		spiltHashTable(&refGenome, &hashTable, numprocs, segStartID);

		for (int i = 1; i < numprocs; i++) {
			nCounterStart = segStartID[i - 1];
			nCounterSize = segStartID[i] - segStartID[i - 1] + 1;
			nIndexSize = hashTable.counter[segStartID[i]] - hashTable.counter[segStartID[i - 1]] + 1;
			MPI_Send(&nCounterStart, 1, MPI_UINT32_T, i, i, MPI_COMM_WORLD);
			MPI_Send(&nCounterSize, 1, MPI_UINT32_T, i, i, MPI_COMM_WORLD);
			MPI_Send(&nIndexSize, 1, MPI_UINT32_T, i, i, MPI_COMM_WORLD);
			MPI_Send(&(hashTable.counter[segStartID[i - 1]]), nCounterSize, MPI_UINT32_T, i, i, MPI_COMM_WORLD);
			MPI_Send(&(hashTable.index[hashTable.counter[segStartID[i - 1]]]), nIndexSize, MPI_UINT32_T, i, i,
					MPI_COMM_WORLD);

			MPI_Send(&refGenome.nRefSize, 1, MPI_UINT32_T, i, i, MPI_COMM_WORLD);
			MPI_Send(refGenome.refSeq, refGenome.nRefSize, MPI_CHAR, i, i, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&nCounterStart, 1, MPI_UINT32_T, 0, rank, MPI_COMM_WORLD, &status);
		MPI_Recv(&nCounterSize, 1, MPI_UINT32_T, 0, rank, MPI_COMM_WORLD, &status);
		MPI_Recv(&nIndexSize, 1, MPI_UINT32_T, 0, rank, MPI_COMM_WORLD, &status);
		hashTable.counter = (uint32_t *) malloc(nCounterSize * sizeof(uint32_t));
		hashTable.index = (uint32_t *) malloc(nIndexSize * sizeof(uint32_t));
		hashTable.nSizeCounter = nCounterSize;
		hashTable.nSizeIndex = nIndexSize;

		MPI_Recv(hashTable.counter, nCounterSize, MPI_UINT32_T, 0, rank, MPI_COMM_WORLD, &status);
		MPI_Recv(hashTable.index, nIndexSize, MPI_UINT32_T, 0, rank, MPI_COMM_WORLD, &status);

		uint32_t tmp = hashTable.counter[0];
		for(uint32_t i = 0;i < nCounterSize;i++) {
			hashTable.counter[i] -= tmp;
		}
		/************************	
		  char chrtest[40];	
		  sprintf(chrtest, "ranktest%d.txt", rank);
		  ofstream fout(chrtest);
		  for(uint32_t i = 0;i < nCounterSize;i++) {
		  fout << i << " " <<  hashTable.counter[i] << " " << hashTable.nSizeIndex << endl;
		  }

		//ifor(uint32_t i =0;i < nIndexSize;i++) {
		//	fout << i << " " << hashTable.index[i] << endl;
		//}
		fout.close();	
		 ***************************************/
		MPI_Recv(&refGenome.nRefSize, 1, MPI_UINT32_T, 0, rank, MPI_COMM_WORLD, &status);
		refGenome.refSeq = (char *) malloc((refGenome.nRefSize + 1) * sizeof(char));
		MPI_Recv(refGenome.refSeq, refGenome.nRefSize, MPI_CHAR, 0, rank, MPI_COMM_WORLD, &status);
		refGenome.refSeq[refGenome.nRefSize] = 0;
	}

	MPI_Bcast(&readLen, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
	MPI_Bcast(&outputFilenamelen, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

	char * readsSeq;
	/* distribute read sequences table to workers */
	if (rank == 0) {
		uint32_t hashValue, procID;
		vector < vector<string> > reads(numprocs - 1);
		for (uint32_t i = 0; i < vReadsSeq.size(); i++) {
			hashValue = GetHashValue(vReadsSeq[i].substr(0, HASHLEN).c_str());
			procID = lower_bound(segStartID, hashValue);
			reads[procID].push_back(vReadsSeq[i]);
		}

		for (int i = 1; i < numprocs; i++) {
			nReads = reads[i - 1].size();
			MPI_Send(&nReads, 1, MPI_UINT32_T, i, i, MPI_COMM_WORLD);
			readsSeq = (char *) malloc(sizeof(char) * nReads * readLen);
			int k = 0;
			for (uint32_t j = 0; j < nReads; j++) {
				for (uint32_t l = 0; l < readLen; l++) {
					readsSeq[k++] = reads[i - 1][j][l];
				}
			}
			MPI_Send(readsSeq, nReads * readLen, MPI_CHAR, i, i, MPI_COMM_WORLD);
			MPI_Send(outputFilename, outputFilenamelen, MPI_CHAR, i, i, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&nReads, 1, MPI_UINT32_T, 0, rank, MPI_COMM_WORLD, &status);
		readsSeq = (char *) malloc(sizeof(char) * nReads * readLen);
		MPI_Recv(readsSeq, nReads * readLen, MPI_CHAR, 0, rank, MPI_COMM_WORLD, &status);
		MPI_Recv(outputFilename, outputFilenamelen, MPI_CHAR, 0, rank, MPI_COMM_WORLD, &status);
		outputFilename[outputFilenamelen] = 0;

		char * oneRead;
		oneRead = (char *) malloc(sizeof(char) * (readLen + 1));
		char outPut[100];	
		sprintf(outPut, "%s_parallel_np%d_rank%d", outputFilename, numprocs, rank);
		ofstream fout(outPut);

		pair < uint32_t, uint32_t > ret;
		int k = 0;
		for (uint32_t i = 0; i < nReads; i++) {
			uint32_t l;
			for (l = 0; l < readLen; l++) {
				oneRead[l] = readsSeq[k++];
			}
			oneRead[l] = 0;
			MappingOneRead(&refGenome, &hashTable, oneRead, ret, nCounterStart);
			fout << ">" << oneRead << endl;
			for(uint32_t j = ret.first; j <= ret.second;j++) {
				fout << hashTable.index[j] << " ";	
			}
			fout << endl;
		}
		fout.close();
		free(readsSeq);
	}
	free(refGenome.refSeq);
	free(hashTable.counter);
	free(hashTable.index);
	if(rank == 0) {
		char command[100];
		sprintf(command, "cat %s_parallel_np%d_rank*  >%s_parallel_np%d_all", outputFilename, numprocs, outputFilename, numprocs);
		system(command);
		//sprintf(rmcommand, "/bin/rm %s_parallel_rank*", outputFilename);
		//system(rmcommand);
	}
}

int main(int argc, char* argv[]) {
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double start_time = MPI_Wtime();
	Parallel_Matching(argc, argv);

	if(rank == 0) printf("%.2lfs time\n", MPI_Wtime() - start_time);
	MPI_Finalize();

	return 0;
}
