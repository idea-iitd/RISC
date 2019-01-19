/*
 * experiments.c
 *
 *  Created on: 25-Nov-2018
 *      Author: jithin
 */

#include "experiments.h"

//#define __TESTING

#ifdef __TESTING

int kStart_G 	= 0;
int kEnd_G 		= 3;
int k_G[]		= {1, 50, 1000};
int rSize_G 	= 3;
double r_G[] 	= {0.4, 0.85, 0.95};

int CHECK_SOLUTION_G= 1;

#else

int kStart_G 	= 0;
int kEnd_G 		= 8;
int k_G[]		= {1, 5, 10, 20, 50, 100, 500, 1000};
int rSize_G 	= 7;
double r_G[] 	= {0.4, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95};

int CHECK_SOLUTION_G= 0;

#endif

int num_exp_G 	= 30;
int releasePublic=1;

/********************************************************************************/
void _experiments_runTopK (void *index, DataBinary *data, DataBinary *dataQueries,
		int k, int numExp, char *rFname, workerFunctions_type workerFunction) {
	G_LOCATION;
	int qSize = dataQueries->size;

	//run over each query

	long long totalTime=0;

	FILE *fp=fopen (rFname,"a");

	if (fp) {
		printf ("The results of the query will be appended to the file : %s\n",rFname);
		fprintf(fp, "###########################\n");
	} else
		perror(rFname);


	double pruned=0, unPruned=0;
	int soultion[k];
	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		arrayListtype *queryFP = dataBinary_getFingerPrint(dataQueries, queryIndex);

		struct timespec tstart={0,0}, tend={0,0};

		clock_gettime(CLOCK_MONOTONIC, &tstart);

		minHeapInttype *solutionHeap = workerFunction.runTopK (index, queryFP, k, &pruned, &unPruned);

		clock_gettime(CLOCK_MONOTONIC, &tend);

		long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);

		totalTime += nsec;

		if (fp) {
			minHeapInt_getValues(solutionHeap, soultion);
			fprintf (fp,"Solution to queryID %d\n",queryIndex);
			for (int i=0; i<k; i++) {
				u_long molID = soultion[i];

				int molName = data->molIDsORG[molID];
				fprintf (fp,"%d\n",molName);
			}
		}
		minHeapInt_free(solutionHeap);

	}
	printf ("K : %d :",k);
	//printf("Average time : %.2lf n.sec\n", (double)totalTime/numExp);
	printf(" %.2lf m.sec, \n", (double)totalTime/numExp/1000/1000);
#ifdef G_GET_PRUNING_RATE
	printf(" %.2lf pruned candidate size\n", pruned/numExp);
	printf(" %.2lf un pruned candidate size\n", unPruned/numExp);
#endif

	if (fp)
		fclose(fp);
}

void _experiments_runRange (void *index, DataBinary *data, DataBinary *dataQueries,
		double r, int numExp, char *rFname, workerFunctions_type workerFunction) {
	G_LOCATION;
	int qSize = dataQueries->size;

	//run over each query

	long long totalTime=0;
	double numSolutions=0;

	FILE *fp=fopen (rFname,"a");

	if (fp) {
		printf ("The results of the query will be appended to the file : %s\n",rFname);
		fprintf(fp, "###########################\n");
	} else
		perror(rFname);

	double pruned=0, unPruned=0;
	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		arrayListtype *queryFP = dataBinary_getFingerPrint(dataQueries, queryIndex);

		struct timespec tstart={0,0}, tend={0,0};

		clock_gettime(CLOCK_MONOTONIC, &tstart);

		arrayListtype *solutionList =  workerFunction.runRange (index, queryFP , r, &pruned, &unPruned);

		int n = solutionList->size;
		numSolutions += n;

		clock_gettime(CLOCK_MONOTONIC, &tend);

		long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
		totalTime += nsec;

		if (fp) {
			arrayList_Iterator_type *iterator = arrayList_getIterator(solutionList);

			fprintf (fp,"Solution to queryID %d\n",queryIndex);
			while (iterator->nextAvailable) {
				u_long molID = arrayList_IteratorGetNext(iterator);

				int molName = data->molIDsORG[molID];
				fprintf (fp,"%d\n",molName);
			}

			arrayList_IteratorFree(iterator);
		}

		free_arrayListtype(solutionList);
	}
	printf ("range : %.2f :",r);
	//printf("Average time : %.2lf n.sec\n", (double)totalTime/numExp);
	printf(" %.2lf m.sec\n", (double)totalTime/numExp/1000/1000);
#ifdef G_GET_PRUNING_RATE
	printf(" %.2lf pruned candidate size\n", pruned/numExp);
	printf(" %.2lf un pruned candidate size\n", unPruned/numExp);
#endif

	if (fp)
		fclose(fp);
}

 void _experiments_runTopK_NB (void *index, DataNonBinary *data, DataNonBinary *dataQueries,
		 int k, int numExp, char *rFname, workerFunctions_NB_type workerFunction) {
	G_LOCATION;
	int qSize = dataQueries->size;

	//run over each query

	long long totalTime=0;

	FILE *fp=fopen (rFname,"a");

	if (fp) {
		printf ("The results of the query will be appended to the file : %s\n",rFname);
		fprintf(fp, "###########################\n");
	} else
		perror(rFname);

	double pruned=0, unPruned=0;
	int soultion[k];

	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		MolNonBinary *queryFP = dataNonBinary_getFingerPrint(dataQueries, queryIndex);

		struct timespec tstart={0,0}, tend={0,0};

		clock_gettime(CLOCK_MONOTONIC, &tstart);

		minHeapInttype *solutionHeap = workerFunction.runTopK (index, queryFP, k, &pruned, &unPruned);

		clock_gettime(CLOCK_MONOTONIC, &tend);

		long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);

		totalTime += nsec;
		if (fp) {
			minHeapInt_getValues(solutionHeap, soultion);
			fprintf (fp,"Solution to queryID %d\n",queryIndex);
			for (int i=0; i<k; i++) {
				u_long molID = soultion[i];

				char *molName = data->ids[molID];
				fprintf (fp,"%s\n",molName);
			}
		}
		minHeapInt_free(solutionHeap);
	}
	printf ("K : %d :",k);
	//printf("Average time : %.2lf n.sec\n", (double)totalTime/numExp);
	printf(" %.2lf m.sec, \n", (double)totalTime/numExp/1000/1000);
#ifdef G_GET_PRUNING_RATE
	printf(" %.2lf pruned candidate size\n", pruned/numExp);
	printf(" %.2lf un pruned candidate size\n", unPruned/numExp);
#endif

	if (fp)
		fclose(fp);
}

void _experiments_runRange_NB (void *index, DataNonBinary *data, DataNonBinary *dataQueries,
		double r, int numExp, char *rFname, workerFunctions_NB_type workerFunction) {
	G_LOCATION;
	int qSize = dataQueries->size;

	//run over each query

	long long totalTime=0;

	FILE *fp=fopen (rFname,"a");

	if (fp) {
		printf ("The results of the query will be appended to the file : %s\n",rFname);
		fprintf(fp, "###########################\n");
	} else
		perror(rFname);

	double pruned=0, unPruned=0;
	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		MolNonBinary *queryFP = dataNonBinary_getFingerPrint(dataQueries, queryIndex);

		struct timespec tstart={0,0}, tend={0,0};

		clock_gettime(CLOCK_MONOTONIC, &tstart);

		arrayListtype *solutionList =  workerFunction.runRange (index, queryFP , r, &pruned, &unPruned);

		clock_gettime(CLOCK_MONOTONIC, &tend);

		long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
		if (fp) {
			arrayList_Iterator_type *iterator = arrayList_getIterator(solutionList);

			fprintf (fp,"Solution to queryID %d\n",queryIndex);
			while (iterator->nextAvailable) {
				u_long molID = arrayList_IteratorGetNext(iterator);

				char *molName = data->ids[molID];
				fprintf (fp,"%s\n",molName);
			}

			arrayList_IteratorFree(iterator);
		}
		free_arrayListtype(solutionList);

		totalTime += nsec;
	}
	printf ("range : %.2f :",r);
	//printf("Average time : %.2lf n.sec\n", (double)totalTime/numExp);
	printf(" %.2lf m.sec\n", (double)totalTime/numExp/1000/1000);
#ifdef G_GET_PRUNING_RATE
	printf(" %.2lf pruned candidate size\n", pruned/numExp);
	printf(" %.2lf un pruned candidate size\n", unPruned/numExp);
#endif

	if (fp)
		fclose(fp);
}

void _experiments_runTopK_small (void *index, DataSmall *data, DataSmall *dataQueries, int k,
		int numExp, char *rFname, workerFunctions_type_small workerFunction) {
	G_LOCATION;
	int qSize = dataQueries->size;

	//run over each query

	long long totalTime=0;
	double r=0;

	FILE *fp=fopen (rFname,"a");

	if (fp) {
		printf ("The results of the query will be appended to the file : %s\n",rFname);
		fprintf(fp, "###########################\n");
	} else
		perror(rFname);

	double pruned=0, unPruned=0;
	int soultion[k];

	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		fingerPrintIntType *queryFP = data_getFingerPrint(dataQueries, queryIndex);

		struct timespec tstart={0,0}, tend={0,0};

		clock_gettime(CLOCK_MONOTONIC, &tstart);

		minHeapInttype *solutionHeap = workerFunction.runTopK (index, queryFP, k, &pruned, &unPruned);

		clock_gettime(CLOCK_MONOTONIC, &tend);

		double minSim = minHeapInt_peekKey(solutionHeap);
		r += minSim;

		long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
		totalTime += nsec;

		if (fp) {
			minHeapInt_getValues(solutionHeap, soultion);
			fprintf (fp,"Solution to queryID %d\n",queryIndex);
			for (int i=0; i<k; i++) {
				u_long molID = soultion[i];

				char *molName = data->ids[molID];
				fprintf (fp,"%s\n",molName);
			}
		}
		minHeapInt_free(solutionHeap);

	}
	printf ("K : %d :",k);
	//printf("Average time : %.2lf n.sec\n", (double)totalTime/numExp);
	printf(" %.2lf m.sec, \n", (double)totalTime/numExp/1000/1000);
#ifdef G_GET_PRUNING_RATE
	printf(" %.2lf pruned candidate size\n", pruned/numExp);
	printf(" %.2lf un pruned candidate size\n", unPruned/numExp);
#endif

	if (fp)
		fclose(fp);
}

void _experiments_runRange_small (void *index, DataSmall *data, DataSmall *dataQueries, double r,
		int numExp, char *rFname, workerFunctions_type_small workerFunction) {
	G_LOCATION;
	int qSize = dataQueries->size;

	//run over each query

	long long totalTime=0;
	double numSolutions=0;

	FILE *fp=fopen (rFname,"a");

	if (fp) {
		printf ("The results of the query will be appended to the file : %s\n",rFname);
		fprintf(fp, "###########################\n");
	} else
		perror(rFname);

	double pruned=0, unPruned=0;
	for (int queryIndex=0; (queryIndex<qSize && queryIndex < numExp); queryIndex++) {
		fingerPrintIntType *queryFP = data_getFingerPrint(dataQueries, queryIndex);

		struct timespec tstart={0,0}, tend={0,0};

		clock_gettime(CLOCK_MONOTONIC, &tstart);

		arrayListtype *solutionList =  workerFunction.runRange (index, queryFP , r, &pruned, &unPruned);

		int n = solutionList->size;
		numSolutions += n;

		clock_gettime(CLOCK_MONOTONIC, &tend);


		long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
		totalTime += nsec;
		if (fp) {
			arrayList_Iterator_type *iterator = arrayList_getIterator(solutionList);

			fprintf (fp,"Solution to queryID %d\n",queryIndex);
			while (iterator->nextAvailable) {
				u_long molID = arrayList_IteratorGetNext(iterator);

				char *molName = data->ids[molID];
				fprintf (fp,"%s\n",molName);
			}

			arrayList_IteratorFree(iterator);
		}
		free_arrayListtype(solutionList);
	}
	printf ("range : %.2f :",r);
	//printf("Average time : %.2lf n.sec\n", (double)totalTime/numExp);
	printf(" %.2lf m.sec\n", (double)totalTime/numExp/1000/1000);
#ifdef G_GET_PRUNING_RATE
	printf(" %.2lf pruned candidate size\n", pruned/numExp);
	printf(" %.2lf un pruned candidate size\n", unPruned/numExp);
#endif

	if (fp)
		fclose(fp);
}

void _experiments_work_helper (DataBinary *data, DataBinary *dataQueries,
		workerFunctions_type workerFunction, char *fName_1) {

	//initialize index
	if (!workerFunction.init_index) {
		printf ("Not implemented\n");
		return;
	}

	struct timespec tstart={0,0}, tend={0,0};
	clock_gettime(CLOCK_MONOTONIC, &tstart);

	void *index = workerFunction.init_index (data);

	clock_gettime(CLOCK_MONOTONIC, &tend);
	long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
	printf("index creation time %.2lf m.sec, ", (double)nsec/1000/1000);

	while (1) {
		//Which experiment do you want to run
		enum options_ExperimentsTypes  experiment = options_getExperimentsTypes ();
		char rFname[200];

		switch (experiment) {
		case TopK :
		{
			for (int j=kStart_G; j<kEnd_G; j++) {
				sprintf (rFname,"%s%s_%d.csv", fName_1, options_getExperminetnName(experiment),k_G[j]);
				_experiments_runTopK  (index, data, dataQueries, k_G[j], num_exp_G, rFname, workerFunction);
			}
		}
		break;
		case Range :
		{
			for (int j=0; j<rSize_G; j++) {
				int x = floor(r_G[j]*100);
				sprintf (rFname,"%s%s_%d.csv", fName_1, options_getExperminetnName(experiment),x);
				_experiments_runRange (index, data, dataQueries, r_G[j], num_exp_G, rFname, workerFunction);
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

void _experiments_nonBinary_work_helper (DataNonBinary *data, DataNonBinary *dataQueries,
		workerFunctions_NB_type workerFunction, char *fName_1) {

	//initialize index
	if (!workerFunction.init_index) {
		printf ("Not implemented\n");
		return;
	}

	struct timespec tstart={0,0}, tend={0,0};
	clock_gettime(CLOCK_MONOTONIC, &tstart);

	void *index = workerFunction.init_index (data);

	clock_gettime(CLOCK_MONOTONIC, &tend);
	long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
	printf("index creation time %.2lf m.sec, ", (double)nsec/1000/1000);

	while (1) {
		//Which experiment do you want to run
		enum options_ExperimentsTypes  experiment = options_getExperimentsTypes ();
		char rFname[200];

		switch (experiment) {
		case TopK :
		{
			for (int j=kStart_G; j<kEnd_G; j++) {
				sprintf (rFname,"%s%s_%d.csv", fName_1, options_getExperminetnName(experiment),k_G[j]);
				_experiments_runTopK_NB  (index, data, dataQueries, k_G[j], num_exp_G, rFname, workerFunction);
			}
		}
		break;
		case Range :
		{
			for (int j=0; j<rSize_G; j++) {
				int x = floor(r_G[j]*100);
				sprintf (rFname,"%s%s_%d.csv", fName_1, options_getExperminetnName(experiment),x);
				_experiments_runRange_NB (index, data, dataQueries, r_G[j], num_exp_G, rFname, workerFunction);
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

void _experiments_work_helper_small (DataSmall *data, DataSmall *dataQueries,
		workerFunctions_type_small workerFunction, char *fName_1) {

	//initialize index
	if (!workerFunction.init_index) {
		printf ("Not implemented\n");
		return;
	}

	struct timespec tstart={0,0}, tend={0,0};
	clock_gettime(CLOCK_MONOTONIC, &tstart);

	void *index = workerFunction.init_index (data);

	clock_gettime(CLOCK_MONOTONIC, &tend);
	long nsec = (tend.tv_nsec + 1.0e9*tend.tv_sec) - (tstart.tv_nsec + 1.0e9*tstart.tv_sec);
	printf("index creation time %.2lf m.sec, ", (double)nsec/1000/1000);

	while (1) {
		//Which experiment do you want to run
		enum options_ExperimentsTypes  experiment = options_getExperimentsTypes ();
		char rFname[200];

		switch (experiment) {
		case TopK :
		{
			for (int j=kStart_G; j<kEnd_G; j++) {
				sprintf (rFname,"%s%s_%d.csv", fName_1, options_getExperminetnName(experiment),k_G[j]);
				_experiments_runTopK_small (index, data, dataQueries, k_G[j], num_exp_G, rFname, workerFunction);
			}
		}
		break;
		case Range :
		{
			for (int j=0; j<rSize_G; j++) {
				int x = floor(r_G[j]*100);
				sprintf (rFname,"%s%s_%d.csv", fName_1, options_getExperminetnName(experiment),x);
				_experiments_runRange_small(index, data, dataQueries, r_G[j], num_exp_G, rFname, workerFunction);
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

/********************************************************************************/
void notImplemented (void *source, char *fname) {
	fprintf (stderr, "Function not implemented\n");
}

workerFunctions_type _blank_Worker_BF () {
	workerFunctions_type worker;

	worker.dataDistribution = NULL;
	worker.init_index 	= NULL;
	worker.free_index	= NULL;
	worker.runRange		= NULL;
	worker.runTopK		= NULL;
	return worker;
}

workerFunctions_type_small _blank_Worker () {
	workerFunctions_type_small worker;

	worker.dataDistribution = NULL;
	worker.init_index 	= NULL;
	worker.free_index	= NULL;
	worker.runRange		= NULL;
	worker.runTopK		= NULL;
	return worker;
}

void experiments_work (DataBinary *data, DataBinary *dataQueries, char *fName_1) {
	if (!data->_isSorted)
			helper_err(__FILE__, __LINE__, "data is not sorted");

	while (1) {
		//which indexing scheme do you want to use
		enum options_IndexTypes  indexingScheme = options_getIndexTypes ();

		char rNamePartial[200];
		sprintf (rNamePartial,"%s%s_", fName_1, options_getIndexingSchemeName(indexingScheme));

		workerFunctions_type workerFunction;

		switch (indexingScheme) {
		case RISC:
			workerFunction = invertedIndex_BF_getWorkerFunction ();
			break;
		case Linear:
			workerFunction = linear_BF_getWorkerFunction ();
			break;
		case AOR:
			workerFunction = AOR_BF_getWorkerFunction ();
			break;
		case DivideSkip:
			workerFunction = divideSkip_BF_getWorkerFunction ();
			break;
		case _EXIT_index :
		case _START_index:
		case _END_index:
			break;
		}

		if (indexingScheme == _EXIT_index)
			break;

		_experiments_work_helper(data, dataQueries, workerFunction, rNamePartial);
	}
}

void experiments_nonBinary_work (DataNonBinary *data, DataNonBinary *dataQueries, char *fName_1) {
	if (!data->_isSorted)
			helper_err(__FILE__, __LINE__, "data is not sorted");

	while (1) {
		//which indexing scheme do you want to use
		enum options_IndexTypes  indexingScheme = options_getIndexTypes ();

		char rNamePartial[200];
		sprintf (rNamePartial,"%s%s_", fName_1, options_getIndexingSchemeName(indexingScheme));

		workerFunctions_NB_type workerFunction;

		switch (indexingScheme) {
		case RISC:
			workerFunction = invertedIndex_NBF_getWorkerFunction ();
			break;
		case Linear:
			workerFunction = linear_NBF_getWorkerFunction ();
			break;
		case AOR:
		case DivideSkip:
			printf ("This indexing scheme is not implemented\n");
			break;
		case _EXIT_index :
		case _START_index:
		case _END_index:
			break;
		}

		if (indexingScheme == _EXIT_index)
			break;

		_experiments_nonBinary_work_helper(data, dataQueries, workerFunction, rNamePartial);
	}
}

void experiments_work_small (DataSmall *data, DataSmall *dataQueries, char *fName_1) {
	if (!data->_isSorted)
			helper_err(__FILE__, __LINE__, "data is not sorted");

	while (1) {
		//which indexing scheme do you want to use
		enum options_IndexTypes  indexingScheme = options_getIndexTypes ();

		char rNamePartial[200];
		sprintf (rNamePartial,"%s%s_", fName_1, options_getIndexingSchemeName(indexingScheme));

		workerFunctions_type_small workerFunction;

		switch (indexingScheme) {
		case RISC:
			workerFunction = invertedIndex_getWorkerFunction ();
			break;
		case Linear:
			workerFunction = linear_getWorkerFunction ();
			break;
		case AOR:
			workerFunction = AOR_getWorkerFunction ();
			break;
		case DivideSkip:
			workerFunction = divideSkip_getWorkerFunction ();
			break;
		case _EXIT_index :
		case _START_index:
		case _END_index:
			break;
		}

		if (indexingScheme == _EXIT_index)
			break;

		_experiments_work_helper_small(data, dataQueries, workerFunction, rNamePartial);
	}
}
