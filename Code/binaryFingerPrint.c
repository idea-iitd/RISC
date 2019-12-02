/*
 * binaryFingerPrint.c
 *
 *  Created on: 16-Dec-2018
 *      Author: jithin
 */

#include "binaryFingerPrint.h"

/****************************************************************/
/********************  Helper   *********************************/
/****************************************************************/
char* _getFname_BF (enum options_datasetTypes d, int dataFile) {
	char *fName = (char*) helper_malloc(500*sizeof(char),
			__FILE__, __LINE__);
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
		case pubChem3DLOOB:
		case pubChem3DE3fp_14:
		case chemFP_DUD:
			printf("This dataset does not have a Molprint2D fingerprint\n");
			break;
		case DUD:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/","dud2006");
			break;
		case pubChem:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/","pubChem");
			break;
		case ZINC:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/","zinc_total");
			break;
		case pubChem_RDKit:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/RDKIT/","pubChem");
			break;
		case ZINC_RDKit:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/RDKIT/","ZINC");
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
		case pubChem3DLOOB:
		case pubChem3DE3fp_14:
		case chemFP_DUD:
			printf("This dataset does not have a Molprint2D fingerprint\n");
			break;
		case DUD:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/Aux/","dud2006/dud2006Test");
			break;
		case pubChem:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/Aux/","pubChem/pubChemTest");
			break;
		case ZINC:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/Aux/","zinc/zincTest");
			break;
		case pubChem_RDKit:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/RDKIT/Aux/","pubChem_Test");
			break;
		case ZINC_RDKit:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/RDKIT/Aux/","ZINC_Test");
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
		case pubChem3DLOOB:
		case pubChem3DE3fp_14:
		case chemFP_DUD:
			printf("This dataset does not have a Molprint2D fingerprint\n");
			break;
		case DUD:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/Aux/","dud2006/features.ice");
			break;
		case pubChem:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/Aux/","pubChem/features.ice");
			break;
		case ZINC:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/molprint2d/Aux/","zinc/features.ice");
			break;
		case pubChem_RDKit:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/RDKIT/Aux/","pubChem_features.ice");
			break;
		case ZINC_RDKit:
			sprintf(fName, "%s%s%s",G_DATA_DIR,"/RDKIT/Aux/","ZINC_features.ice");
			break;
		}
	}

	return fName;
}

int _dataSize_BF (char *fName) {
	char *buffer = NULL;
	size_t bufsiz = 0;
	ssize_t nbytes;
	int n=0;

	FILE *fp = fopen(fName, "r");
	helper_errFP (fp , fName, __FILE__,__LINE__);

	//count the number of molecules available
	while ((nbytes = getline(&buffer, &bufsiz, fp)) != -1) {
		n++;
	}
	free(buffer);
	fclose(fp);
	return n;
}

void _loadData_BF (char *fName, arrayListtype *featureIds,
		DataBinary **data_org, int n, int alterLoad) {
	char *buffer = NULL;
	size_t bufsiz = 0;
	ssize_t nbytes;

	if (n <= 0)
		n = _dataSize_BF(fName);

	//initialize the data
	DataBinary *data = new_dataBinary(n, featureIds->size);

	FILE *fp = fopen(fName, "r");
	helper_errFP (fp , fName, __FILE__,__LINE__);

	arrayList_Iterator_type *iterator = arrayList_getIterator(featureIds);
	minHeapPlain_Type *minHeap_new = new_minHeapPlain_Type(100);
	minHeapPlain_Type *minHeap_old = new_minHeapPlain_Type(100);

	int i=1;
	while ((nbytes = getline(&buffer, &bufsiz, fp)) != -1) {
		if (i > n) {
			char msg[100];
			sprintf(msg, "Not enough memory allocated i=%u n=%u",i,n);
			helper_err(__FILE__, __LINE__, msg);
		}

		u_long feature;
		int val=-1;
		if(alterLoad) {
			if (buffer[0] != '#')
				continue;
			char *str 	= strtok (buffer," ");
			str = strtok ((char *)NULL, " "); //first value is molID

			do {
				u_long feature;
				val = -1;
				sscanf(str, "%lu:%d",&feature,&val);

				if (val > 0) {
					//we want to store the features as sorted feature
					minHeapPlain_add_WO_fail_unique(&minHeap_old, feature);
				}

				str = strtok ((char *)NULL, " ");
			} while (str);

		} else {
			if (buffer[0] != '{')
				continue;
			char *str 	= strtok (buffer,"{");
			str = strtok (str,"}");
			str = strtok (str,"\n");
			str = strtok (str,",");

			do {
				val=-1;
				sscanf(str, "%lu:%d",&feature,&val);
				if (val > 0) {

					//we want to store the features as sorted feature
					minHeapPlain_add_WO_fail_unique(&minHeap_old, feature);
				}

				str = strtok ((char *)NULL, ",");
			} while (str);

		}
		helper_printProgress(i++);

		while (minHeap_old->size) {
			u_long featureOld;
			int count;
			minHeapPlain_popMultiple(minHeap_old, &featureOld, &count);

			if (count > 1) {
				printf("featureOld:%lu, count:%d - \n",featureOld,count);
				helper_err(__FILE__, __LINE__, "invalid feature values");
			}

			int found;
			u_long newFeatureId = arrayList_IteratorFind_Binary (iterator, featureOld, &found);

			if (!found) {
				helper_err(__FILE__, __LINE__, "invalid feature values, the fingerprint has a feature not originally present in the indexing phase");
			}

			minHeapPlain_add_WO_fail_unique(&minHeap_new, newFeatureId);
		}

		//move the iterator back
		arrayList_resetIterator(iterator, featureIds);

		//add the features as a entry in the data
		dataBinary_addMolecule(data, minHeap_new);
	}

	minHeapPlain_free(minHeap_new);
	minHeapPlain_free(minHeap_old);

	printf("\n");
	free(buffer);
	fclose(fp);
	*data_org = data;
}

/**
 * return a sorted list of featureids used in the dataset
 * num : the number of molecules in the dataset
 */
arrayListtype* _loadFeatureIds (char *fName, char *fNameFeatures, int *num) {
	unsigned int size = 10000;

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
		int feature;
		while ((nbytes = getline(&buffer, &bufsiz, fp_f)) != -1) {
			sscanf(buffer, "%d",&feature);
			arrayList_put(featureIds, feature);
		}
		*num = n;
		free(buffer);
		fclose(fp_f);
		return featureIds;
	} else {
		minHeapPlain_Type *minHeap = new_minHeapPlain_Type(size);

		char *buffer = NULL;
		size_t bufsiz = 0;
		ssize_t nbytes;
		int n = 0;

		FILE *fp = fopen(fName, "r");
		helper_errFP (fp , fName, __FILE__,__LINE__);

		while ((nbytes = getline(&buffer, &bufsiz, fp)) != -1) {
			//count the number of molecules available
			n++;
			helper_printProgress(n);

			if (buffer[0] != '{')
				continue;
			//load features to the min heap,
			char *str 	= strtok (buffer,"{");
			str = strtok (str,"}");
			str = strtok (str,"\n");
			str = strtok (str,",");
			do {
				int feature, val=-1;
				sscanf(str, "%d:%d",&feature,&val);
				if (val > 0) {
					minHeapPlain_add_WO_fail_unique(&minHeap, feature);
				}

				str = strtok ((char *)NULL, ",");
			} while (str);
		}

		free(buffer);
		fclose(fp);
		printf ("\n");

		*num = n;

		//write featureIds to file
		fprintf(stderr,	"Do not terminate the program now, it will not function correctly in future\n");
		fprintf(stderr,	"If the program terminates/crashes in between kindly delete %s\n", fNameFeatures);
		fp_f = fopen(fNameFeatures, "w");
		fprintf(fp_f, "size_data %d\n",n);

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

/************************* Linear ************************************/

minHeapInttype* _runTopKHelper_Linear_BF (void *index, arrayListtype *queryFP , int k,
		double *pruned, double *unPruned) {
	DataBinary *data = (DataBinary *)index;
	minHeapInttype *solutionHeapLinear = newMinHeapInttype(k);
	for (int i=0; i<data->size; i++) {
		arrayListtype *dataFP = dataBinary_getFingerPrint(data, i);
		double sim = dataBinary_tanimotoSimilarity_fast (queryFP, dataFP);

		minHeapInt_addIfKeyHigher(solutionHeapLinear, sim, i);
	}

	return solutionHeapLinear;
}


void _runTopK_Linear_BF (void *source, DataBinary *dataQueries, int k, int numExp, char *resultFname) {
	DataBinary *data  = (DataBinary *)source;
	int qSize = dataQueries->size;

	//run over each query

	long long totalTime=0;
	double r=0;
	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		arrayListtype *queryFP = dataBinary_getFingerPrint(dataQueries, queryIndex);

		struct timespec tstart={0,0}, tend={0,0};

		clock_gettime(CLOCK_MONOTONIC, &tstart);

		minHeapInttype *solutionHeap = _runTopKHelper_Linear_BF (data, queryFP , k, NULL, NULL);

		clock_gettime(CLOCK_MONOTONIC, &tend);

		r += minHeapInt_peekKey(solutionHeap);
		minHeapInt_free(solutionHeap);

		long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
		//printf ("%d, %ld,\n", queryIndex, nsec);

		totalTime += nsec;
	}
	printf ("K : %d :",k);
	//printf("Average time : %.2lf n.sec\n", (double)totalTime/numExp);
	printf(" %.2lf m.sec, ", (double)totalTime/numExp/1000/1000);
	printf(" %.2f : Average similarity threshold\n", r/numExp);
}

arrayListtype* _runRangeHelper_Linear_BF(void *index, arrayListtype *queryFP , double r,
		double *pruned, double *unPruned) {
	DataBinary *data = (DataBinary *)index;
	arrayListtype *solutionLinear = new_arrayListtype(1000);
	for (int i=0; i<data->size; i++) {
		arrayListtype *dataFP = dataBinary_getFingerPrint(data, i);
		double sim = dataBinary_tanimotoSimilarity_fast (queryFP, dataFP);

		if (sim >= r) {
			arrayList_put(solutionLinear, i);
		}
	}
	return solutionLinear;
}

void _runRangeSmall_Linear_BF (void *source, DataBinary *dataQueries, double r, int numExp, char *resultFname) {
	DataBinary *data = (DataBinary *) source;
	int qSize = dataQueries->size;

	//run over each query

	long long totalTime=0;
	double numSolutions=0;

	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		arrayListtype *queryFP = dataBinary_getFingerPrint(dataQueries, queryIndex);

		struct timespec tstart={0,0}, tend={0,0};

		clock_gettime(CLOCK_MONOTONIC, &tstart);

		arrayListtype *solutionList =  _runRangeHelper_Linear_BF(data, queryFP , r, NULL, NULL);

		numSolutions += solutionList->size;

		clock_gettime(CLOCK_MONOTONIC, &tend);

		free_arrayListtype(solutionList);

		long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
		//printf ("%d, %ld,\n", queryIndex, nsec);

		totalTime += nsec;
	}
	printf ("range : %.2f :",r);
	//printf("Average time : %.2lf n.sec\n", (double)totalTime/numExp);
	printf(" %.2lf m.sec", (double)totalTime/numExp/1000/1000);
	printf(" %.2lf average solution size\n", numSolutions/numExp);
}

void runLinear_BF (DataBinary *data, DataBinary *dataQueries, char *rnamePartial) {

	void *source = data;

	while (1){
		//Which experiment do you want to run
		enum options_ExperimentsTypes  experiment = options_getExperimentsTypes ();

		char rFname[200];

		switch (experiment) {
		case TopK :
		{
			for (int j=kStart_G; j<kEnd_G; j++) {
				sprintf (rFname,"%s%s_%d.csv", rnamePartial, options_getExperminetnName(experiment),k_G[j]);
				_runTopK_Linear_BF 		(source, dataQueries, k_G[j], num_exp_G, rFname);
			}
		}
		break;
		case Range :
		{
			for (int j=0; j<rSize_G; j++) {
				int x = floor(r_G[j]*100);
				sprintf (rFname,"%s%s_%d.csv", rnamePartial, options_getExperminetnName(experiment),x);
				_runRangeSmall_Linear_BF(source, dataQueries, r_G[j], num_exp_G, rFname);
			}
		}
		break;
		case _EXIT_exp :
			break;
		case _START_exp:
		case _END_exp:
			break;
		}
		if (experiment == _EXIT_exp)
			break;
	}
}

/****************************************************************/
/****************************************************************/
void binaryFingerPrint_run (enum options_datasetTypes d) {
	//
	DataBinary *data, *dataQueries;

	//read the data from the file
	int dataFile = 1;
	char *fNameData  	= _getFname_BF (d, dataFile);
	dataFile = 0;
	char *fNamequeries 	= _getFname_BF (d, dataFile);
	dataFile = -1;
	char *fNameFeatures	= _getFname_BF (d, dataFile);

	int n;
	arrayListtype *featureIds = _loadFeatureIds (fNameData, fNameFeatures, &n); //sorted list of featureids used in the dataset
	
	printf ("size : %d\n",n);

	int alterLoad=0;
	_loadData_BF (fNamequeries, featureIds, &dataQueries, 0, alterLoad);
	_loadData_BF (fNameData,  	featureIds, &data, n, alterLoad);
	dataBinary_sort (data);

	char rNamePartial[200];
	sprintf (rNamePartial,"%s%s_",G_RESULT_DIR,options_getDataSetName(d));

	experiments_work (data, dataQueries, rNamePartial);
}

void binaryFingerPrint_run_public () {
	printf ("What is the dataset name : ");
	char dataName[100];
	scanf ("%s", dataName);

	char fNameTarget[100], fNameQueries[100], fNameFeatures[100];
	sprintf (fNameTarget,  "%sBinary/%s",G_DATA_DIR,dataName);
	sprintf (fNameQueries, "%sBinary/%s_test",G_DATA_DIR,dataName);
	sprintf (fNameFeatures,"%sBinary/Aux/%s_features",G_DATA_DIR,dataName);

	printf ("The target  file is %s\n", fNameTarget);
	printf ("The queries file is %s\n", fNameQueries);

	//
	DataBinary *data, *dataQueries;

	int n;
	arrayListtype *featureIds = _loadFeatureIds (fNameTarget, fNameFeatures, &n); //sorted list of featureids used in the dataset

	printf ("size : %d\n",n);

	int alterLoad=0;
	_loadData_BF (fNameQueries, featureIds, &dataQueries, 0, alterLoad);
	_loadData_BF (fNameTarget, 	featureIds, &data, n, alterLoad);
	dataBinary_sort (data);

	char rNamePartial[200];
	sprintf (rNamePartial,"%s%s_",G_RESULT_DIR,dataName);

	experiments_work (data, dataQueries, rNamePartial);
}

void* _linear_BFinit_index (DataBinary *data) {
	return data;
}

void _linear_BFfree_index  (void *index) {
	//DO nothing
}

workerFunctions_type linear_BF_getWorkerFunction () {
	workerFunctions_type worker;

	worker.dataDistribution = notImplemented;
	worker.init_index 	= _linear_BFinit_index;
	worker.free_index	= _linear_BFfree_index;
	worker.runRange		= _runRangeHelper_Linear_BF;
	worker.runTopK		= _runTopKHelper_Linear_BF;
	return worker;
}
