/*
 *the main function
 *
 */

#include "chemFP.h"

/**
 * get the name of the file where data is stored
 */
char * _getFname (char *numBits, int targetFile) {
	char *fName = (char*) helper_malloc(200*sizeof(char),
			__FILE__, __LINE__);

	if (targetFile)
		sprintf(fName, "%s%s%s%s",G_DATA_DIR,"CHEMFP_benchmark/benchmark/targets_",numBits,".fps");
	else
		sprintf(fName, "%s%s%s%s",G_DATA_DIR,"CHEMFP_benchmark/benchmark/queries_",numBits,".fps");
	return fName;
}

/**
 * count the number of molecules in the file
 */
int _dataSetSize (char *fName) {
	char *buffer = NULL;
	size_t bufsiz = 0;
	ssize_t nbytes;
	int n=0;

	FILE *fp = fopen(fName, "r");
	helper_errFP (fp , fName, __FILE__,__LINE__);

	while ((nbytes = getline(&buffer, &bufsiz, fp)) != -1) {
		//this is not a comment line
	    if (buffer[0] != '#')
	    	n++;
	}
	free(buffer);

	fclose(fp);

	return n;
}

/**
 * Load the data into a common format
 */
void _loadDataSmall (char *fName, DataSmall *data, int numBits, int sLimit) {
	int numMolecules = _dataSetSize(fName);
	int sizeI = ceil ((double)numBits/NUMBITS_fingerPrintIntType); //the length of the fingerprintI array

	//allocate memory
	data_init (data, numMolecules, sizeI);

	FILE *fp = fopen(fName, "r");

	helper_errFP (fp , fName, __FILE__,__LINE__);

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	while ((read = getline(&line, &len, fp)) != -1) {
		//this is a comment line
		if (line[0] == '#')
			continue;

		if (sLimit >= 0) {
			if (sLimit-- == 0)
				break;
		}

		//break the line into 2 fingerprint and id
		char *fingerPrintS 	= strtok (line, 		"\t");
		char *id 			= strtok ((char *)NULL, "\n");

		//printf ("%s : %s\n", fingerPrintS, id);
		//printf ("%s,", fingerPrintS);

		//load the fingerprint into data
		data_addFingerprint(data, fingerPrintS, id);
	}
	printf ("\n");

	free(line);
	fclose(fp);
}

/*****************************************************************/

struct _solutionStat {
	int intersection;
	int unionSize;
	int length;
	double sim;
	int molID;
};

struct _solutionStat* _new_solutionStat (int intersection, int unionSize, int length,
	double sim, int molID) {
	struct _solutionStat *stat = (struct _solutionStat*) helper_malloc(sizeof(struct _solutionStat), __FILE__, __LINE__);

	stat->intersection 	= intersection;
	stat->unionSize		= unionSize;
	stat->length		= length;
	stat->sim			= sim;
	stat->molID			= molID;

	return stat;
}

void printResultStat (void *element, double key, int printHeader) {

	if (printHeader) {
		printf ("KEY, Similarity, L_q, len, intersection, union, molID\n");
	} else {
		struct _solutionStat *stat = (struct _solutionStat*) element;
		int L_q = stat->unionSize+stat->intersection-stat->length;
		printf ("%f,%f,%3d,%3d,%3d,%3d,%d\n", key, stat->sim, L_q, stat->length,
				stat->intersection, stat->unionSize, stat->molID);
	}
}

void runResultDistribution (DataSmall *data, DataSmall *dataQueries, int numExp) {
	int dSize = data->size;
	int qSize = dataQueries->size;

	int fpLen = data->fingerPrintLen;

	int *featureIndex = (int *)helper_malloc((fpLen*NUMBITS_fingerPrintIntType+1)*sizeof(int),__FILE__, __LINE__);

	minHeaptype *solutionHeap = newMinHeaptype(100, printResultStat);
	//run over each query
	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		fingerPrintIntType *queryFP = data_getFingerPrint(dataQueries, queryIndex);


		int L_q = data_getIndexOf1Bits (featureIndex, queryFP, fpLen);

		printf ("%d,",queryIndex);

		//run through each fingerprints
		for (int dataIndex=0; dataIndex<dSize; dataIndex++) {
			fingerPrintIntType *dataFP = data_getFingerPrint(data, dataIndex);

			int l = data_getIndexOf1Bits (featureIndex, dataFP, fpLen);

			double intersection = data_getIntersection(queryFP, dataFP, fpLen);
			int u = l+L_q-intersection;
			double sim = intersection/u;

			struct _solutionStat *stat = _new_solutionStat(intersection, u, l, sim, dataIndex);

			int notadded = minHeap_addIfKeyHigher(solutionHeap, sim, stat);

			if (notadded)
				free (stat);
		}

		//we have the top-K results
		minHeap_print(solutionHeap);
	}
	free (featureIndex);
}

void chemFP_run_public () {
	printf ("What is the dimension of the fingerprint : ");
	u_long numBits;
	scanf ("%lu", &numBits);

	printf ("You have entered the dimension to be %lu\n",numBits);

	char fNameTarget[100], fNameQueries[100];
	sprintf (fNameTarget, "%sBinary/targets_%lu.fps",G_DATA_DIR,numBits);
	sprintf (fNameQueries,"%sBinary/queries_%lu.fps",G_DATA_DIR,numBits);

	printf ("The target  file is %s\n", fNameTarget);
	printf ("The queries file is %s\n", fNameQueries);

	//data
	DataSmall data, dataQueries;

	_loadDataSmall (fNameQueries, &dataQueries, numBits, -1);
	_loadDataSmall (fNameTarget,  &data		  , numBits, -1);
	data_sort(&data);
	//data_printHex (&dataQueries);

	char rNamePartial[200];
	sprintf (rNamePartial,"%sChemFP_%lu_",G_RESULT_DIR,numBits);

	experiments_work_small (&data, &dataQueries, rNamePartial);

}
void chemFP_run(char *fileID, int numBits) {
	//data
	DataSmall data, dataQueries;

	//read the data from the file
	int tartgetFile = 1;
	char *fNameTarget  = _getFname (fileID, tartgetFile);
	tartgetFile = 0;
	char *fNamequeries = _getFname (fileID, tartgetFile);

	//load the data from file
	//int numBits = atoi(numBitsString);

	_loadDataSmall (fNamequeries, &dataQueries, numBits, -1);
	_loadDataSmall (fNameTarget,  &data		  , numBits, -1);
	data_sort(&data);
	//data_printHex (&dataQueries);

	char rNamePartial[200];
	sprintf (rNamePartial,"%s%s_",G_RESULT_DIR,fileID);

	experiments_work_small (&data, &dataQueries, rNamePartial);
}
/************************************************************************************/

arrayListtype* _runRangeHelper_Linear(void *index, fingerPrintIntType *queryFP , double r,
		double *pruned, double *unPruned) {
	DataSmall *data = (DataSmall *)index;
	u_long dSize = data->size;
	int fpLen = data->fingerPrintLen;

	arrayListtype *solutionLinear = new_arrayListtype(1000);
	for (int i=0; i<dSize; i++) {
		fingerPrintIntType *dataFP = data_getFingerPrint(data, i);
		double sim = data_tanimotoSimilarity (queryFP, dataFP, fpLen);

		if (sim >= r) {
			arrayList_put(solutionLinear, i);
		}
	}
	return solutionLinear;
}

minHeapInttype* _runTopKHelper_Linear (void *index, fingerPrintIntType *queryFP , int k,
		double *pruned, double *unPruned) {
	DataSmall *data = (DataSmall *)index;
	u_long dSize = data->size;
	int fpLen = data->fingerPrintLen;

	minHeapInttype *solutionHeapLinear = newMinHeapInttype(k);
	for (u_long i=0; i<dSize; i++) {
		fingerPrintIntType *dataFP = data_getFingerPrint(data, i);
		double sim = data_tanimotoSimilarity (queryFP, dataFP, fpLen);

		minHeapInt_addIfKeyHigher(solutionHeapLinear, sim, i);
	}

	return solutionHeapLinear;
}

void* _linear_init_index (DataSmall *data) {
	return data;
}

void _linear_free_index  (void *index) {
	//DO nothing
}

workerFunctions_type_small linear_getWorkerFunction () {
	workerFunctions_type_small worker;

	worker.dataDistribution = notImplemented;
	worker.init_index 	= _linear_init_index;
	worker.free_index	= _linear_free_index;
	worker.runRange		= _runRangeHelper_Linear;
	worker.runTopK		= _runTopKHelper_Linear;
	return worker;
}
