/*
 * nonBinaryFingerPrint.c
 *
 *  Created on: 30-Dec-2018
 *      Author: jithin
 */

#include "nonBinaryFingerPrint.h"

/****************************************************************/
/********************  Helper   *********************************/
/****************************************************************/
char* _getFname_NBF (enum options_datasetTypes d, int dataFile) {
	char *fName = (char*) helper_malloc(200*sizeof(char), __FILE__, __LINE__);

	if (dataFile == 1) {
		switch (d) {
		case _START_data:
		case _EXIT_data:
		case chemFP_166:
		case chemFP_881:
		case chemFP_1024:
		case chemFP_2048:
		case chemFP_ZINC_400K:
		case _END_data:
		case DUD:
		case pubChem:
		case ZINC:
		case pubChem_RDKit:
		case ZINC_RDKit:
		case chemFP_DUD:
			printf("This dataset does not have a 3D fingerprint\n");
			break;
		case pubChem3DE3fp_14:
			sprintf(fName, "%s%s",G_DATA_DIR,"pubChem3D_e3fp_14/data.e3fp");
			break;
		case pubChem3DLOOB:
			sprintf(fName, "%s%s",G_DATA_DIR,"LOOB/pubChem3D.loob");
			break;
		}
	} else if (dataFile == 0){
		switch (d) {
		case _START_data:
		case _EXIT_data:
		case chemFP_166:
		case chemFP_881:
		case chemFP_1024:
		case chemFP_2048:
		case chemFP_ZINC_400K:
		case _END_data:
		case DUD:
		case pubChem:
		case ZINC:
		case pubChem_RDKit:
		case ZINC_RDKit:
		case chemFP_DUD:
			printf("This dataset does not have a 3D fingerprint\n");
			break;
		case pubChem3DE3fp_14:
			sprintf(fName, "%s%s",G_DATA_DIR,"pubChem3D_e3fp_14/Aux/test.e3fp");
			break;
		case pubChem3DLOOB:
			sprintf(fName, "%s%s",G_DATA_DIR,"LOOB/Aux/pubChem3D_Test");
			break;
		}
	} else if (dataFile == -1){
		//Feature file
		switch (d) {
		case _START_data:
		case _EXIT_data:
		case chemFP_166:
		case chemFP_881:
		case chemFP_1024:
		case chemFP_2048:
		case chemFP_ZINC_400K:
		case _END_data:
		case DUD:
		case pubChem:
		case ZINC:
		case pubChem_RDKit:
		case ZINC_RDKit:
		case chemFP_DUD:
			printf("This dataset does not have a 3D fingerprint\n");
			break;
		case pubChem3DE3fp_14:
			sprintf(fName, "%s%s",G_DATA_DIR,"pubChem3D_e3fp_14/Aux/features.ice");
			break;
		case pubChem3DLOOB:
			sprintf(fName, "%s%s",G_DATA_DIR,"LOOB/Aux/Pubchem3D_features.ice");
			break;
		}
	}
	return fName;
}

u_long _NBF_featureLoader_helper(char fName[100], u_long n,
		minHeapPlain_Type** minHeap) {
	char* buffer = NULL;
	size_t bufsiz = 0;
	ssize_t nbytes;
	FILE* fp = fopen(fName, "r");
	helper_errFP(fp, fName, __FILE__, __LINE__);
	while ((nbytes = getline(&buffer, &bufsiz, fp)) != -1) {
		if (buffer[0] != '#') {
			//its not a valid fingerprint
			continue;
		}
		//count the number of molecules available
		n++;
		helper_printProgress(n);
		//load features to the min heap,
		char* str = strtok(buffer, " ");
		str = strtok((char*) NULL, " "); //first value is molID
		do {
			u_long feature;
			int val = -1;
			sscanf(str, "%lu:%d", &feature, &val);
			if (val > 0) {
				minHeapPlain_add_WO_fail_unique(&*minHeap, feature);
			}
			str = strtok((char*) NULL, " ");
		} while (str);
	}
	free(buffer);
	fclose(fp);
	return n;
}

void _loadData_NBF_helper (char *fName, arrayListtype *featureIds, DataNonBinary *data) {
	char* buffer = NULL;
	size_t bufsiz = 0;
	ssize_t nbytes;
	FILE* fp = fopen(fName, "r");
	helper_errFP(fp, fName, __FILE__, __LINE__);

	arrayList_Iterator_type *iterator = arrayList_getIterator(featureIds);

	int n = 0;
	while ((nbytes = getline(&buffer, &bufsiz, fp)) != -1) {
		if (buffer[0] != '#') {
			//its not a valid fingerprint
			continue;
		}
		//count the number of molecules available
		n++;
		helper_printProgress(n);

		dataNonBinary_addMolecule(data, buffer, iterator, featureIds);
	}
	free(buffer);
	fclose(fp);
}

u_long _dataSize_NBF (char *fName) {
	char* buffer = NULL;
	size_t bufsiz = 0;
	ssize_t nbytes;
	FILE* fp = fopen(fName, "r");
	helper_errFP(fp, fName, __FILE__, __LINE__);

	u_long n = 0;

	while ((nbytes = getline(&buffer, &bufsiz, fp)) != -1) {
		if (buffer[0] != '#') {
			//its not a valid fingerprint
			continue;
		}
		//count the number of molecules available
		n++;
	}
	free(buffer);
	fclose(fp);
	return n;
}

void _loadData_NBF (char *dName, arrayListtype *featureIds,
		DataNonBinary **data_org, int n, int loadFromDir) {
	if (n <= 0)
		n = _dataSize_NBF(dName);

	//initialize the data
	DataNonBinary *data = new_dataNonBinary(n, featureIds->size);

	if (loadFromDir) {
		DIR *directory;
		struct dirent* file;

		directory = opendir(dName);
		if (directory == NULL) {
			helper_perror(__FILE__, __LINE__, dName);
			exit(2);
		}

		while ((file=readdir(directory)) != NULL) {
			int len = strlen(file->d_name);
			if (len > 6) {
				char *extn = &(file->d_name[len-5]);

				if (strcmp(extn, ".e3fp") == 0) {
					char fName[100];
					sprintf (fName, "%s/%s",dName,file->d_name);

					_loadData_NBF_helper (fName, featureIds, data);
				}
			}
		}

		closedir(directory);
	} else {
		_loadData_NBF_helper (dName, featureIds, data);
	}

	*data_org = data;
}
/**
 * return a sorted list of featureids used in the dataset
 * num : the number of molecules in the dataset
 */
arrayListtype* _loadFeatureIds_NBF (char *dName, char *fNameFeatures, int *num, int loadFromDir) {
	unsigned int size = 100000;

	//read the feature from file if exists
	FILE *fp_f = fopen(fNameFeatures, "r");
	if (fp_f) {
		// the feature file exits
		char *buffer = NULL;
		size_t bufsiz = 0;
		ssize_t nbytes;

		int n;
		{
			int x = fscanf(fp_f, "size_data %d\n",&n);
			if (x < 1)
				helper_perror(__FILE__, __LINE__, "");
		}
		//populate the feature in to the array list
		arrayListtype *featureIds = new_arrayListtype(size);
		u_long feature;
		while ((nbytes = getline(&buffer, &bufsiz, fp_f)) != -1) {
			sscanf(buffer, "%lu",&feature);
			arrayList_put(featureIds, feature);
		}
		*num = n;
		free(buffer);
		fclose(fp_f);
		return featureIds;
	} else {
		minHeapPlain_Type *minHeap = new_minHeapPlain_Type(size);
		u_long n = 0;

		if (loadFromDir) {
			DIR *directory;
			struct dirent* file;

			directory = opendir(dName);
			if (directory == NULL) {
				helper_perror(__FILE__, __LINE__, dName);
				exit(2);
			}

			while ((file=readdir(directory)) != NULL) {
				int len = strlen(file->d_name);
				if (len > 6) {
					char *extn = &(file->d_name[len-5]);

					if (strcmp(extn, ".e3fp") == 0) {
						char fName[100];
						sprintf (fName, "%s/%s",dName,file->d_name);

						n = _NBF_featureLoader_helper(fName, n, &minHeap);
					}
				}
			}

			closedir(directory);
		} else {
			n = _NBF_featureLoader_helper(dName, n, &minHeap);
		}
		printf ("\n");
		*num = n;

		fprintf(stderr,	"Do not terminate now the program will not function correctly\n");
		fprintf(stderr,	"If the program terminates/crashes in between kindly delete %s\n", fNameFeatures);
		//write featureIds to file
		fp_f = fopen(fNameFeatures, "w");
		fprintf(fp_f, "size_data %lu\n",n);

		//populate the feature in to the array list
		arrayListtype *featureIds = new_arrayListtype(size);
		while (minHeap->size > 0) {
			u_long value;
			int count;
			minHeapPlain_popMultiple(minHeap, &value, &count);
			arrayList_put(featureIds, value);
			fprintf(fp_f, "%lu\n",value);
		}
		fclose(fp_f);
		fprintf(stderr,	"It is safe to terminate the program if you wish to do so\n");
		minHeapPlain_free(minHeap);
		return featureIds;
	}
}

void _generate_Test (char* dName, char *tName) {
	FILE *fp_o = fopen (tName,"w");

	DIR *directory;
	struct dirent* file;

	directory = opendir(dName);
	if (directory == NULL) {
		helper_perror(__FILE__, __LINE__, dName);
		exit(2);
	}

	srand(time(NULL));   // Initialization, should only be called once.

	int n = 0;
	while ((file=readdir(directory)) != NULL) {
		int len = strlen(file->d_name);
		if (len > 6) {
			char *extn = &(file->d_name[len-5]);

			if (strcmp(extn, ".e3fp") == 0) {
				char fName[100];
				sprintf (fName, "%s/%s",dName,file->d_name);

				do {
					char* buffer = NULL;
					size_t bufsiz = 0;
					ssize_t nbytes;
					FILE* fp = fopen(fName, "r");
					helper_errFP(fp, fName, __FILE__, __LINE__);

					while ((nbytes = getline(&buffer, &bufsiz, fp)) != -1) {
						if (buffer[0] != '#') {
							//its not a valid fingerprint
							continue;
						}
						n++;
						helper_printProgress(n);

						int r = rand()%10000;
						if (r == 1) {
							fprintf(fp_o, "%s\n",buffer);
						}
					}
					free(buffer);
					fclose(fp);

				} while (0);
			}
		}
	}

	closedir(directory);
	fclose(fp_o);
	exit (0);
}

/****************************************************************/
/****************************************************************/
void nonBinaryFingerPrint_run (enum options_datasetTypes d) {
	DataNonBinary *data, *dataQueries;

	//read the data from the file
	int dataFile = 1;
	char *fNameData  	= _getFname_NBF (d, dataFile);
	dataFile = 0;
	char *fNamequeries 	= _getFname_NBF (d, dataFile);
	dataFile = -1;
	char *fNameFeatures	= _getFname_NBF (d, dataFile);

	int n;
	int loadFromDir =0;

	arrayListtype *featureIds = _loadFeatureIds_NBF (fNameData, fNameFeatures, &n, loadFromDir); //sorted list of featureids used in the dataset

	printf ("size : %d\n",n);

	_loadData_NBF (fNamequeries, featureIds, &dataQueries, 0, 0);
	_loadData_NBF (fNameData,  	featureIds, &data, n, loadFromDir);
	dataNonBinary_sort (data);

	char rNamePartial[200];
	sprintf (rNamePartial,"%s%s_",G_RESULT_DIR,options_getDataSetName(d));

	experiments_nonBinary_work (data, dataQueries, rNamePartial);
}

void nonBinaryFingerPrint_run_public () {
	printf ("What is the dataset name : ");
	char dataName[100];
	scanf ("%s", dataName);

	char fNameTarget[100], fNameQueries[100], fNameFeatures[100];
	sprintf (fNameTarget,  "%sNonBinary/%s",G_DATA_DIR,dataName);
	sprintf (fNameQueries, "%sNonBinary/%s_test",G_DATA_DIR,dataName);
	sprintf (fNameFeatures,"%sNonBinary/Aux/%s_features",G_DATA_DIR,dataName);

	printf ("The target  file is %s\n", fNameTarget);
	printf ("The queries file is %s\n", fNameQueries);

	DataNonBinary *data, *dataQueries;

	int n;
	int loadFromDir =0;

	arrayListtype *featureIds = _loadFeatureIds_NBF (fNameTarget, fNameFeatures, &n, loadFromDir); //sorted list of featureids used in the dataset

	printf ("size : %d\n",n);

	_loadData_NBF (fNameQueries,featureIds, &dataQueries, 0, 0);
	_loadData_NBF (fNameTarget,	featureIds, &data, n, loadFromDir);
	dataNonBinary_sort (data);

	char rNamePartial[200];
	sprintf (rNamePartial,"%s_Non_Binary_%s_",G_RESULT_DIR,dataName);

	experiments_nonBinary_work (data, dataQueries, rNamePartial);
}

/*****************************************************************/
/*********              Linear                  ******************/
/*****************************************************************/
void* _linear_NBFinit_index (DataNonBinary *data) {
	return data;
}

void _linear_NBFfree_index  (void *index) {
	//DO nothing
}


arrayListtype* _runRangeHelper_Linear_NBF(void *index, MolNonBinary *queryFP , double r,
		double *pruned, double *unPruned) {
	DataNonBinary *data = (DataNonBinary *)index;
	arrayListtype *solutionLinear = new_arrayListtype(1000);
	for (int i=0; i<data->size; i++) {
		MolNonBinary *dataFP = dataNonBinary_getFingerPrint(data, i);
		double sim = dataNonBinary_tanimotoSimilarity (queryFP, dataFP);

		if (sim >= r) {
			arrayList_put(solutionLinear, i);
		}
	}
	return solutionLinear;
}

minHeapInttype* _runTopKHelper_Linear_NBF (void *index, MolNonBinary *queryFP , int k,
		double *pruned, double *unPruned) {
	DataNonBinary *data = (DataNonBinary *)index;
	minHeapInttype *solutionHeapLinear = newMinHeapInttype(k);
	for (int i=0; i<data->size; i++) {
		MolNonBinary *dataFP = dataNonBinary_getFingerPrint(data, i);
		double sim = dataNonBinary_tanimotoSimilarity (queryFP, dataFP);

		minHeapInt_addIfKeyHigher(solutionHeapLinear, sim, i);
	}

	return solutionHeapLinear;
}

workerFunctions_NB_type linear_NBF_getWorkerFunction () {
	workerFunctions_NB_type worker;

	worker.dataDistribution = notImplemented;
	worker.init_index 	= _linear_NBFinit_index;
	worker.free_index	= _linear_NBFfree_index;
	worker.runRange		= _runRangeHelper_Linear_NBF;
	worker.runTopK		= _runTopKHelper_Linear_NBF;
	return worker;
}

/*****************************************************************/
/*********          Inverted Index              ******************/
/*****************************************************************/
typedef struct{
	int lengthMin;
	int lengthMax;
	arrayListtype_flexible *molIDs;
	arrayListtype_flexible *cummulativeFeatureValue;
	arrayListtype_flexible *featureValue;
	int *_sizeIndexStart;
	int *_sizeIndexEnd;
	int cvMin; 					//Cumulative value start and end within each len
	int cvMax;
}_moleculeList_NB;

typedef struct{
	_moleculeList_NB **mList;
	DataNonBinary *data;
}_invertedIndex_NB;

_moleculeList_NB* _new_moleculeList_NB (int dSize, int lengthMax) {
	_moleculeList_NB *mList 		= (_moleculeList_NB *)helper_calloc(1, sizeof(_moleculeList_NB), __FILE__, __LINE__);
	mList->molIDs 					= new_arrayListtype_flexible(10, dSize);
	mList->cummulativeFeatureValue 	= new_arrayListtype_flexible(10, dSize);
	mList->featureValue 			= new_arrayListtype_flexible(10, dSize);
	mList->_sizeIndexStart 	= (int *)helper_calloc(lengthMax+1, sizeof(int), __FILE__, __LINE__);
	mList->_sizeIndexEnd 	= (int *)helper_calloc(lengthMax+1, sizeof(int), __FILE__, __LINE__);

	for (int i=0; i<=lengthMax; i++)
		mList->_sizeIndexStart[i]=-1;

	mList->lengthMin 		= INT_MAX;
	mList->cvMin 			= INT_MAX;
	return mList;
}

void _add2mlist (u_long molID, u_long cVal, u_long val, u_long len,  _moleculeList_NB *mList) {
	int n = mList->molIDs->size;

	arrayList_flexible_put(&(mList->cummulativeFeatureValue), cVal);
	arrayList_flexible_put(&(mList->featureValue), val);
	arrayList_flexible_put(&(mList->molIDs), molID);

	if (len < mList->lengthMin)
		mList->lengthMin = len;
	if (len > mList->lengthMax)
		mList->lengthMax = len;

	if (cVal < mList->cvMin)
		mList->cvMin = cVal;
	if (cVal > mList->cvMax)
		mList->cvMax = cVal;

	if (mList->_sizeIndexStart[len] < 0)
		mList->_sizeIndexStart[len] = n;
	mList->_sizeIndexEnd[len] = n;
}

void* _II_NBFinit_index (DataNonBinary *data) {
	if (!data->_isSorted) {
		helper_err(__FILE__, __LINE__, "data is not sorted");
	}

	u_long numFeatures 	= data->numFeatures;
	u_long dSize 		= data->size;

	_invertedIndex_NB *index =(_invertedIndex_NB *)helper_calloc(1, 		 sizeof(_invertedIndex_NB), __FILE__, __LINE__);
	_moleculeList_NB **mLists=(_moleculeList_NB **)helper_calloc(numFeatures,sizeof(_moleculeList_NB*), __FILE__, __LINE__);

	for (int i=0; i<numFeatures; i++) {
		mLists[i] = _new_moleculeList_NB(dSize, data->maxLen);
	}
	index->mList = mLists;

	//add molecules
	u_long minLen = data->minLen;
	u_long maxLen = data->maxLen;

	for (int len=minLen; len<=maxLen; len++) {
		int indexStart = data->_sizeIndexStart[len];
		if (indexStart >= 0) {
			int indexEnd = data->_sizeIndexEnd[len];

			for (u_long molID=indexStart; molID<=indexEnd; molID++) {
				 MolNonBinary *fp = dataNonBinary_getFingerPrint(data, molID);

				 u_long cVal 	 = fp->cummulativeFeatureVal;
				 u_long *features= fp->features->data;
				 u_long *values	 = fp->values->data;

				 for (u_long f=0; f<len; f++) {
					 u_long feature = features[f];
					 u_long value 	= values[f];

					 _moleculeList_NB *mList = mLists[feature];

					 _add2mlist (molID, cVal, value, len, mList);
				 }
			}
		}
	}

	index->data = data;
	return index;
}


void _dataDistribution_NBF(void *source, char*fName) {
	FILE *fp = (FILE *) fopen (fName, "w");

	if (!fp)
		return;

	_invertedIndex_NB *index= (_invertedIndex_NB *)source;
	DataNonBinary *data 	= index->data;

	//length distribution
	fprintf (fp,"len, #mol\n");

	for (int len=1; len<=data->maxNumFeatures; len++) {
		int startIndex  	= data->_sizeIndexStart[len];
		if (startIndex >= 0) {
			int endIndex 	= data->_sizeIndexEnd[len];
			int size = endIndex - startIndex + 1;
			fprintf (fp,"%d, %d\n", len, size);
		}
	}

	//feature list distribution
	fprintf (fp,"\nfeatureId, #mol\n");
	u_long dSize = data->size;
	u_long lenLimit = dSize/10;

	arrayListtype *hevayHitters = new_arrayListtype(1000);
	u_long numFeatures = data->numFeatures;
	for (u_long f=0;f<numFeatures;f++) {
		_moleculeList_NB *mlist = index->mList[f];

		u_long size = mlist->molIDs->size;
		fprintf (fp,"%lu, %lu\n",f,size);


		if (size >= lenLimit)
			arrayList_put(hevayHitters, f);
	}

	fclose (fp);

	//print information about heavy hitters

	arrayList_Iterator_type *iterator = arrayList_getIterator(hevayHitters);

	double totalHH = 0;
	int maxHH = 0;

	for (int molId=0; molId < dSize; molId++) {
		MolNonBinary *mol =  dataNonBinary_getFingerPrint (data, molId);

		int n = mol->features->size;

		int numHH=0;
		arrayList_flexible_Iterator_type *fiterator = arrayList_flexible_getIterator(mol->features);
		for (int i=0; i<n; i++){
			int feature = arrayList_flexible_IteratorGetNext(fiterator);

			int found = 0;
			int notenoughelements = arrayList_IteratorFind(iterator, feature, &found);
			if (notenoughelements)
				break;

			if (found) {
				numHH++;
			}
		}
		arrayList_flexible_IteratorFree(fiterator);
		if (numHH > maxHH)
			maxHH = numHH;

		totalHH += numHH;

		arrayList_resetIterator(iterator, hevayHitters);
	}

	printf ("dataset  num HH : %lu\n", hevayHitters->size);
	printf ("mean num HH : %f\n", totalHH/dSize);
	printf ("max  num HH : %d\n", maxHH);

	arrayList_IteratorFree(iterator);
}


void _II_NBFfree_index  (void *index) {
	//DO nothing
}

void _sortFeatures_NBF (u_long *featureIndex, u_long *featureValues,
		MolNonBinary *queryFP,_invertedIndex_NB *index, int L_q) {
	u_long len[L_q];

	u_long *qFeatures= queryFP->features->data;
	u_long *qValuess = queryFP->values->data;
	for (u_long f=0; f<L_q; f++) {
		u_long feature = qFeatures[f];

		_moleculeList_NB *mList = index->mList[feature];
		u_long mSize 	= mList->molIDs->size;

		featureIndex [f] = feature;
		featureValues[f] = qValuess[f];
		len[f]			 = mSize;
	}

	//sort
	for (int i=0;i<L_q; i++) {
		for (int j=i+1;j<L_q; j++) {
			if (len[i] > len[j]) {
				//swap
				u_long temp;

				temp  = len[i];
				len[i]= len[j];
				len[j]= temp;

				temp  			= featureIndex[i];
				featureIndex[i]	= featureIndex[j];
				featureIndex[j]	= temp;

				temp  			= featureValues[i];
				featureValues[i]= featureValues[j];
				featureValues[j]= temp;
			}
		}
	}
}

arrayListtype* _runRangeHelper_II_NBF_1(void *source, MolNonBinary *queryFP , double r,
		double *pruned, double *unPruned) {
	_invertedIndex_NB *index = (_invertedIndex_NB *)source;

	DataNonBinary *data = index->data;
	u_long dSize = data->size;

	int *candidates = (int *)helper_calloc(dSize, sizeof(int), __FILE__, __LINE__);

	u_long L_q = queryFP->features->size;

	u_long featureIndex [L_q];
	u_long featureValues[L_q];
	_sortFeatures_NBF (featureIndex, featureValues, queryFP, index, L_q);
	arrayListtype *solution = new_arrayListtype(1000);

#ifdef G_GET_PRUNING_RATE
	int *allC  = (int *) helper_calloc(data->size, sizeof(int),	__FILE__, __LINE__);
	for (int i=0; i<L_q; i++) {
		u_long feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList_NB *mList= index->mList[feature];
		u_long *molid = mList->molIDs->data;
		u_long msize  = mList->molIDs->size;
		for (int i=0;i<msize; i++) {
			allC[molid[i]]=1;
		}
	}

	for (int i=0; i<data->size; i++) {
		if (allC[i])
			*unPruned = *unPruned + 1;
	}

	free(allC);
#endif


	//The cumulative value of fingerprint
	double Sq = queryFP->cummulativeFeatureVal;

	int ll = 0; 		//lower limit on length of a molecule
	double sl = r * Sq; //lower limit on cumulative value
	//find the lower limit on length
	do {
		//sort the feature values
		u_long featureValues_sorted[L_q];
		for (int i=0;i<L_q; i++) {
			featureValues_sorted[i] = featureValues[i];
		}

		for (int i=0;i<L_q; i++) {
			for (int j=i+1;j<L_q; j++) {
				if (featureValues_sorted[i] < featureValues_sorted[j]) {
					u_long temp = featureValues_sorted[i];
					featureValues_sorted[i] = featureValues_sorted[j];
					featureValues_sorted[j] = temp;
				}
			}
		}

		//Partial Cumulative feature value
		double pQ = 0;

		int i;
		for (i=0; pQ < sl; i++) {
			pQ += featureValues_sorted[i];
		}
		i--;

		ll = i;
	}while(0);


	double lu = floor ((double)Sq/r - Sq +L_q);
	u_long su = floor((double)Sq*(r+1)/r - Sq);
	sl = ceil(sl);

	for (int m=0; m<L_q; m++) {
		u_long feature = featureIndex[m];

		//explore mList
		_moleculeList_NB *mList = index->mList[feature];

		int lMin = mList->lengthMin;
		if (lMin < ll)
			lMin = ll;

		int lMax = mList->lengthMax;
		if (lMax > lu)
			lMax = lu;

		u_long *molIDs 	= mList->molIDs->data;
		u_long *cvs     = mList->cummulativeFeatureValue->data;
		u_long *values  = mList->featureValue->data;

		for (int len=lMin; len<=lMax; len++) {
			int startIndex = mList->_sizeIndexStart[len];
			if (startIndex >= 0) {
				int endIndex = mList->_sizeIndexEnd[len];

				int i=startIndex;

				for (; i<=endIndex; i++){
					if (cvs[i] >= sl)
						break;
				}
				for (; endIndex>=i; endIndex--){
					if (cvs[endIndex] <= su)
						break;
				}

				u_long qFval = featureValues[m];
				for (; i<=endIndex; i++){
					u_long molID =molIDs[i];
					u_long fValue;
					if (qFval < values[i])
						fValue=qFval;
					else
						fValue=values[i];

					candidates[molID]+=fValue;
				}
			}
		}
	}

	do {
		u_long Q = queryFP->cummulativeFeatureVal;
		double alpha = r/(r+1);
		int PRECISION = 1000000;
		alpha = ((double)floor(alpha*PRECISION))/PRECISION;

		for (u_long molID=0; molID<dSize; molID++) {
			if (candidates[molID]){
#ifdef G_GET_PRUNING_RATE
				*pruned = *pruned + 1;
#endif
				MolNonBinary *dataFP = dataNonBinary_getFingerPrint(data, molID);

				double minVal = (alpha * (Q+dataFP->cummulativeFeatureVal));
				if (candidates[molID] >= minVal) {
					arrayList_put(solution, molID);
				}
			}
		}
	} while (0);

	//check the solution
	if (CHECK_SOLUTION_G) {
		static int c=0;
		printf ("%d : comparing the 2 solutions : \n", c++);
		arrayListtype *solutionLinear = new_arrayListtype(1000);
		for (int i=0; i<dSize; i++) {
			MolNonBinary *dataFP = dataNonBinary_getFingerPrint(data, i);
			double sim = dataNonBinary_tanimotoSimilarity (queryFP, dataFP);

			if (sim >= r) {
				arrayList_put(solutionLinear, i);
			}
		}

		//compare the 2 solutions
		u_long sSize = solution->size;

		if (sSize != solutionLinear->size) {
			printf ("II : %lu,%lu Linear; the size of solutions do not match\n",sSize, solutionLinear->size);
			//helper_err(__FILE__, __LINE__, "the size of solutions do not match");
		//} else {
			arrayList_Iterator_type *iteratorI = arrayList_getIterator(solution);
			arrayList_Iterator_type *iteratorL = arrayList_getIterator(solutionLinear);

			for (int i=0; i<sSize; i++) {
				int a = arrayList_IteratorGetNext(iteratorI);
				int b = arrayList_IteratorGetNext(iteratorL);

				if (a != b) {
					printf ("%d,%d; ",a,b);
					helper_err(__FILE__, __LINE__, "Solutions do not match");
				}
			}
			printf("\n");
			arrayList_IteratorFree(iteratorI);
			arrayList_IteratorFree(iteratorL);
		}

		free_arrayListtype(solutionLinear);
	}

	return solution;
}

void _topK_NBF_helper_helper(int len, double sl, u_long su,
		u_long *featureValues, int m, _moleculeList_NB* mList, u_long* cvs,
		u_long* molIDs, u_long* values, int* candidates) {
	int startIndex = mList->_sizeIndexStart[len];
	if (startIndex >= 0) {
		int endIndex = mList->_sizeIndexEnd[len];
		int i = startIndex;
		for (; i <= endIndex; i++) {
			if (cvs[i] >= sl)
				break;
		}
		for (; endIndex >= i; endIndex--) {
			if (cvs[endIndex] <= su)
				break;
		}
		u_long qFval = featureValues[m];
		for (; i <= endIndex; i++) {
			u_long molID = molIDs[i];
			u_long fValue;
			if (qFval < values[i])
				fValue = qFval;
			else
				fValue = values[i];

			candidates[molID] += fValue;
		}
	}
}

void _process_candidates(int len, u_long Q, DataNonBinary* data,
		int* candidates, minHeapInttype* solutionHeap, double *pruned) {
	//process the candidates
	int startIndex = data->_sizeIndexStart[len];
	if (startIndex >= 0) {
		int endIndex = data->_sizeIndexEnd[len];
		for (int molID = startIndex; molID <= endIndex; molID++) {
			int c = candidates[molID];
			if (c) {
#ifdef G_GET_PRUNING_RATE
				*pruned = *pruned + 1;
#endif

				MolNonBinary* dataFP = dataNonBinary_getFingerPrint(data,
						molID);
				double sim = (double) c
						/ (Q + dataFP->cummulativeFeatureVal - c);
				minHeapInt_addIfKeyHigher(solutionHeap, sim, molID);
			}
		}
	}
}

minHeapInttype* _runTopKHelper_II_NBF_1 (void *source, MolNonBinary *queryFP , int k,
		double *pruned, double *unPruned) {
	_invertedIndex_NB *index = (_invertedIndex_NB *)source;

	DataNonBinary *data = index->data;
	u_long dSize = data->size;

	int *candidates = (int *)helper_calloc(dSize, sizeof(int), __FILE__, __LINE__);
	double defaultSim = 1.0/data->maxCV;

	minHeapInttype *solutionHeap = newMinHeapInttype_withDefault(k, defaultSim);

	u_long L_q = queryFP->features->size;

	u_long featureIndex [L_q];
	u_long featureValues[L_q];
	_sortFeatures_NBF (featureIndex, featureValues, queryFP, index, L_q);
#ifdef G_GET_PRUNING_RATE
	int *allC  = (int *) helper_calloc(data->size, sizeof(int),	__FILE__, __LINE__);
	for (int i=0; i<L_q; i++) {
		u_long feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList_NB *mList= index->mList[feature];
		u_long *molid = mList->molIDs->data;
		u_long msize  = mList->molIDs->size;
		for (int i=0;i<msize; i++) {
			allC[molid[i]]=1;
		}
	}

	for (int i=0; i<data->size; i++) {
		if (allC[i])
			*unPruned = *unPruned + 1;
	}

	free(allC);
#endif

	//The cumulative value of fingerprint
	double Sq = queryFP->cummulativeFeatureVal;

	//Partial Cumulative feature value
	double pQ = 0;
	int ll = 0; 		//lower limit on length of a molecule


	//find the lower limit on length
	//sort the feature values
	u_long featureValues_sorted[L_q];
	for (int i=0;i<L_q; i++) {
		featureValues_sorted[i] = featureValues[i];
	}

	for (int i=0;i<L_q; i++) {
		for (int j=i+1;j<L_q; j++) {
			if (featureValues_sorted[i] < featureValues_sorted[j]) {
				u_long temp = featureValues_sorted[i];
				featureValues_sorted[i] = featureValues_sorted[j];
				featureValues_sorted[j] = temp;
			}
		}
	}

	u_long Q = queryFP->cummulativeFeatureVal;

	int moveLeft=1, moveRight=1;
	int leftLen	= L_q;
	int rightLen= leftLen+1;

	while (moveLeft || moveRight){
		if (moveLeft) {
			double r = minHeapInt_peekKey_or_default(solutionHeap);
			double sl = r * Sq; //lower limit on cumulative value
			u_long su = floor((double)Sq*(r+1)/r - Sq);

			do {
				int flag=0;
				for (; pQ < sl; ll++) {
					flag=1;
					pQ += featureValues_sorted[ll];
				}
				if (flag)
					ll--;
			}while (0);

			if (leftLen < ll) {
				moveLeft = 0;
			} else {
				for (int m=L_q-1; m>=0; m--) {
					u_long feature = featureIndex[m];

					//explore mList
					_moleculeList_NB *mList = index->mList[feature];

					if (leftLen < mList->lengthMin)
						continue;
					if (leftLen > mList->lengthMax)
						continue;

					u_long *molIDs 	= mList->molIDs->data;
					u_long *cvs     = mList->cummulativeFeatureValue->data;
					u_long *values  = mList->featureValue->data;

					_topK_NBF_helper_helper(leftLen, sl, su, featureValues,
							m, mList, cvs, molIDs, values, candidates);

				}

				//process the candidates
				_process_candidates(leftLen, Q, data,candidates, solutionHeap, pruned);
				leftLen--;
			}
		}

		if (moveRight) {
			double r = minHeapInt_peekKey_or_default(solutionHeap);
			double lu = floor ((double)Sq/r - Sq +L_q);
			if (lu > data->maxLen)
				lu = data->maxLen;
			double sl = r * Sq; //lower limit on cumulative value
			u_long su = floor((double)Sq*(r+1)/r - Sq);

			if (rightLen > lu) {
				moveRight = 0;
			} else {
				for (int m=L_q-1; m>=0; m--) {
					u_long feature = featureIndex[m];

					//explore mList
					_moleculeList_NB *mList = index->mList[feature];

					if (rightLen > mList->lengthMax)
						continue;
					if (rightLen < mList->lengthMin)
						continue;

					u_long *molIDs 	= mList->molIDs->data;
					u_long *cvs     = mList->cummulativeFeatureValue->data;
					u_long *values  = mList->featureValue->data;

					_topK_NBF_helper_helper(rightLen, sl, su, featureValues,
							m, mList, cvs, molIDs, values, candidates);
				}

				//process the candidates
				_process_candidates(rightLen, Q, data,candidates, solutionHeap, pruned);
				rightLen++;
			}
		}
	}
	//check the solution
	if (CHECK_SOLUTION_G) {
		static int c;
		printf ("%d comparing the 2 solutions\n",c++);
		minHeapInttype *solutionHeapLinear = newMinHeapInttype(k);
		for (int i=0; i<dSize; i++) {
			MolNonBinary *dataFP = dataNonBinary_getFingerPrint(data, i);
			double sim = dataNonBinary_tanimotoSimilarity (queryFP, dataFP);

			minHeapInt_addIfKeyHigher(solutionHeapLinear, sim, i);
		}

		//compare the 2 solutions
		double linearIndex[k];
		double riscIndex  [k];

		//int linearIndex[k];
		//int riscIndex  [k];
		//minHeapInt_getValues(solutionHeap, 		riscIndex);
		//minHeapInt_getValues(solutionHeapLinear,linearIndex);

		minHeapInt_getKeys(solutionHeap, 		riscIndex);
		minHeapInt_getKeys(solutionHeapLinear,linearIndex);

		//sort
		for (int i=0; i<k; i++) {
			for (int j=i+1; j<k; j++) {
				if (riscIndex[i] > riscIndex[j]) {
					double t = riscIndex[i];
					riscIndex[i] = riscIndex[j];
					riscIndex[j] = t;
				}
			}
		}
		for (int i=0; i<k; i++) {
			for (int j=i+1; j<k; j++) {
				if (linearIndex[i] > linearIndex[j]) {
					double t = linearIndex[i];
					linearIndex[i] = linearIndex[j];
					linearIndex[j] = t;
				}
			}
		}
		for (int i=0; i<k; i++) {
			if (linearIndex[i] != riscIndex[i]) {
				for (int j=0; j<k; j++)
					printf ("%f,%f;\n", linearIndex[j], riscIndex[j]);

				printf ("\n\nthe linear heap\n");
				minHeapInt_print (solutionHeapLinear);
				helper_err(__FILE__ , __LINE__, "Solution not matching");
			}
		}

		minHeapInt_free(solutionHeapLinear);
	}

	free (candidates);
	return solutionHeap;
}

workerFunctions_NB_type invertedIndex_NBF_getWorkerFunction () {
	workerFunctions_NB_type worker;

	worker.dataDistribution = _dataDistribution_NBF;
	worker.init_index 	= _II_NBFinit_index;
	worker.free_index	= _II_NBFfree_index;

	worker.runRange		= _runRangeHelper_II_NBF_1;
	worker.runTopK		= _runTopKHelper_II_NBF_1;
	return worker;
}
