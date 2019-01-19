/*
 * AOR.c
 *
 *  Created on: 12-Dec-2018
 *      Author: jithin
 */

#include "AOR.h"

#define _ARRAY_BUCKET_SIZE_AOR 1000

int _numBins_AOR = _G_AOR;

typedef struct {
	arrayListtype 		*molList;
	fingerPrintIntType 	*referencePoint;
}_bin_AOR;

typedef struct {
	int numBins;
	_bin_AOR **bins;
}_bucket_AOR;

typedef struct {
	int numFeatures; 		// the number of features in the index
	DataSmall *data;
	_bucket_AOR **buckets;	//the array of index entries
}_index_AOR;

/******************* Initializers ************************/
/************************ index ************************/
_bin_AOR* _new_bin_AOR (int fpLen) {
	_bin_AOR *bin 		= (_bin_AOR *) 			helper_calloc(1, 	sizeof(_bin_AOR), __FILE__, __LINE__);
	bin->referencePoint = (fingerPrintIntType *)helper_calloc(fpLen,sizeof (fingerPrintIntType),
							__FILE__, __LINE__);
	bin->molList = new_arrayListtype(_ARRAY_BUCKET_SIZE_AOR);

	return bin;
}

void _free_bin_AOR (_bin_AOR *bin) {

	free_arrayListtype(bin->molList);
	free (bin->referencePoint);
	free (bin);
}

_bucket_AOR*_new_bucket_AOR (DataSmall *data, int len) {
	_bucket_AOR *bucket = NULL;
	int startIndex  	= data->_sizeIndexStart[len];
	if (startIndex >= 0) {
		int endIndex 	= data->_sizeIndexEnd[len];

		bucket 		= (_bucket_AOR *) helper_calloc (1, 	 sizeof(_bucket_AOR),__FILE__, __LINE__);
		bucket->bins= (_bin_AOR **)   helper_calloc (_numBins_AOR,sizeof(_bin_AOR*),  __FILE__, __LINE__);

		int size = endIndex-startIndex+1;

		int n = ceil ((double)size/_numBins_AOR);

		int fpIndex	=startIndex;
		int m 		= n;
		int fpLen 	= data->fingerPrintLen;
		_bin_AOR *bin = _new_bin_AOR (fpLen);
		bucket->bins[bucket->numBins++] = bin;

		do {
			fingerPrintIntType *fingerPrint = data_getFingerPrint(data, fpIndex);
			arrayList_put(bin->molList, fpIndex);

			//update the reference point
			for (int i=0; i<fpLen; i++) {
				bin->referencePoint[i] |= fingerPrint[i];
			}

			fpIndex++;
			if (fpIndex >endIndex) {
				break;
			}

			m--;
			if (m <= 0) {
				m = n;
				bin = _new_bin_AOR (fpLen);
				bucket->bins[bucket->numBins++] = bin;
			}
		}while (1);
	}
	return bucket;
}

void _free_bucket_AOR (_bucket_AOR *bucket) {
	if (bucket) {
		for (int i=0; i<_numBins_AOR; i++) {
			_bin_AOR *bin = bucket->bins[i];

			if (bin) {
				_free_bin_AOR (bin);
			}
		}
		free (bucket->bins);
		free (bucket);
	}
}

void* _init_index_AOR (DataSmall *data) {
	_index_AOR *index = (_index_AOR *) helper_calloc (1, sizeof(_index_AOR),
			__FILE__, __LINE__);
	int numBits = data->fingerPrintLen * NUMBITS_fingerPrintIntType;
	index->numFeatures 	= numBits;
	index->data			= data;

	_bucket_AOR **buckets = (_bucket_AOR **) helper_malloc ((numBits+1) * sizeof(_bucket_AOR **),
			__FILE__, __LINE__);

	for (int len=1; len<=numBits; len++) {
		buckets[len] = _new_bucket_AOR (data, len);
	}
	index->buckets = buckets;

	return index;
}

void _free_index_AOR (void *source) {
	_index_AOR *index = (_index_AOR *) source;

	int numBits = index->numFeatures;
	for (int len=1; len<=numBits; len++) {
		_free_bucket_AOR (index->buckets[len]);
	}

	free (index->buckets);
	free (index);
}
/**********************************************/

/******************* Helpers *****************************/
arrayListtype* _runRangeSmallHelper_AOR(void *source, fingerPrintIntType *queryFP , double r,
		double *pruned, double *unPruned) {
	_index_AOR *index = (_index_AOR *) source;
	DataSmall *data   = index->data;
	int dSize 		  = data->size;
	int fpLen 		  = data->fingerPrintLen;
	int numFeatures   = index->numFeatures;

	int L_q = data_getNum1Bits(queryFP, fpLen);

	arrayListtype *solution = new_arrayListtype(_ARRAY_BUCKET_SIZE_AOR);

	//get the feature lists

	int L_l = floor (r*L_q);
	int L_u = ceil ((double)L_q/r)+1;
	if (L_u > numFeatures+1)
		L_u = numFeatures+1;

	for (int l = L_l; l< L_u; l++){
		_bucket_AOR *bucket = index->buckets[l];

		if (!bucket)
			continue;

		int B_plus_C = l + L_q;

		for (int i=0; i<bucket->numBins; i++) {
			_bin_AOR *bin 					= bucket->bins[i];
			arrayListtype *arr 				= bin->molList;
			fingerPrintIntType *referenceFP = bin->referencePoint;

			double gamma = data_getIntersection(queryFP, referenceFP, fpLen);

			gamma = gamma*(r+1);
			if (B_plus_C < gamma)
				continue; //prune
			gamma = gamma/r;

			if (B_plus_C > (gamma+0.2)) //floating point error
				continue; //prune

			int size = arr->size;
			arrayList_Iterator_type *iterator = arrayList_getIterator(arr);
			while (size--) {
				u_long molId = arrayList_IteratorGetNext(iterator);

				fingerPrintIntType *dataFP = data_getFingerPrint(data, molId);
				double sim = data_tanimotoSimilarity(queryFP, dataFP, fpLen);

				if (sim >= r) {
					arrayList_put(solution, molId);
				}
			}
			arrayList_IteratorFree(iterator);
		}
	}

	if (CHECK_SOLUTION_G) {
		int *featureIndex = (int *)helper_malloc((numFeatures+1)*sizeof(int),__FILE__, __LINE__);
		//verify the solution
		static int c;

		printf ("%3d)comparing the 2 solutions : ",c++);
		arrayListtype *solutionLinear = new_arrayListtype(_ARRAY_BUCKET_SIZE_AOR);
		for (int i=0; i<dSize; i++) {
			fingerPrintIntType *dataFP = data_getFingerPrint(data, i);
			double sim = data_tanimotoSimilarity (queryFP, dataFP, fpLen);

			if (sim >= r) {
				arrayList_put(solutionLinear, i);
			}
		}

		//compare the 2 solutions
		u_long sSize = solution->size;

		if (sSize != solutionLinear->size) {
			printf ("II : %lu,%lu Linear; the size of solutions do not match\n",sSize, solutionLinear->size);
			//helper_err(__FILE__, __LINE__, "the size of solutions do not match");
		}else {
			arrayList_Iterator_type *iteratorI = arrayList_getIterator(solution);
			arrayList_Iterator_type *iteratorL = arrayList_getIterator(solutionLinear);

			for (int i=0; i<sSize; i++) {
				u_long a = arrayList_IteratorGetNext(iteratorI);
				u_long b = arrayList_IteratorGetNext(iteratorL);



				fingerPrintIntType *dataFP;
				dataFP 		= data_getFingerPrint	(data, a);
				int n1 		= data_getIndexOf1Bits 	(featureIndex, dataFP, fpLen);
				int inter_1 = data_getIntersection	(queryFP, dataFP, fpLen);
				dataFP 		= data_getFingerPrint	(data, b);
				int n2 		= data_getIndexOf1Bits 	(featureIndex, dataFP, fpLen);
				int inter_2 = data_getIntersection	(queryFP, dataFP, fpLen);

				if (a != b) {
					printf ("Point 2 : ");
					for (int j=0; j<n2; j++) {
						printf ("%d,",featureIndex[j]);
					}
					printf ("\n");
					printf ("id:len,intersection; \n");
					printf ("%lu:%d,%d; %lu:%d,%d; \n",a,n1,inter_1,b,n2,inter_2);
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

minHeapInttype* _runTopKSmallHelper_AOR(void *source, fingerPrintIntType *queryFP , int k,
		double *pruned, double *unPruned) {
	_index_AOR *index 	= (_index_AOR *) source;
	DataSmall *data 	= index->data;
	int dSize 			= data->size;
	int fpLen 			= data->fingerPrintLen;
	int numFeatures 	= index->numFeatures;

	minHeapInttype *solutionHeap = newMinHeapInttype(k);

	int L_q = data_getNum1Bits(queryFP, fpLen);

	double r = (double)1.0/numFeatures;

	int L_l = floor (r*L_q);
	int L_u = ceil ((double)L_q/r)+1;
	if (L_u > numFeatures+1)
		L_u = numFeatures+1;

	int checkLength = 0;

	for (int l = L_q-1; l>0; l--){
		_bucket_AOR *bucket = index->buckets[l];

		if (!bucket)
			continue;

		int B_plus_C = l + L_q;


		for (int i=0; i<bucket->numBins; i++) {
			if (checkLength) {
				r = minHeapInt_peekKey(solutionHeap);
				if (r==0)
					r = (double)1.0/numFeatures;;
				checkLength = 0;

				L_l = floor (r*L_q);
				if (l < L_l) {
					l = 0; // we can stop completely
					break;
				}
			}

			_bin_AOR *bin 					= bucket->bins[i];
			fingerPrintIntType *referenceFP = bin->referencePoint;

			double gamma = data_getIntersection(queryFP, referenceFP, fpLen);

			gamma = gamma*(r+1);
			if (B_plus_C < gamma-0.2)
				continue; //prune
			gamma = gamma/r;

			if (B_plus_C > (gamma+0.2)) //floating point error
				continue; //prune

			arrayListtype *arr 				= bin->molList;
			int size = arr->size;
			arrayList_Iterator_type *iterator = arrayList_getIterator(arr);
			while (size--) {
				u_long molId = arrayList_IteratorGetNext(iterator);

				fingerPrintIntType *dataFP = data_getFingerPrint(data, molId);
				double sim = data_tanimotoSimilarity(queryFP, dataFP, fpLen);

				minHeapInt_addIfKeyHigher(solutionHeap, sim, molId);
			}
			checkLength=1;
		}
	}

	for (int l = L_q; l< L_u; l++){
		_bucket_AOR *bucket = index->buckets[l];

		if (!bucket)
			continue;

		int B_plus_C = l + L_q;


		for (int i=0; i<bucket->numBins; i++) {
			if (checkLength) {
				r = minHeapInt_peekKey(solutionHeap);
				if (r==0)
					r = (double)1.0/numFeatures;;
				checkLength = 0;
				L_u = ceil ((double)L_q/r)+1;
				if (L_u > numFeatures+1)
					L_u = numFeatures+1;

				if (l >= L_u)
					break;
			}

			_bin_AOR *bin 					= bucket->bins[i];
			fingerPrintIntType *referenceFP = bin->referencePoint;
			double gamma = data_getIntersection(queryFP, referenceFP, fpLen);

			gamma = gamma*(r+1);
			if (B_plus_C < gamma-0.2)
				continue; //prune
			gamma = gamma/r;

			if (B_plus_C > (gamma+0.2)) //floating point error
				continue; //prune

			arrayListtype *arr 				= bin->molList;
			int size = arr->size;
			arrayList_Iterator_type *iterator = arrayList_getIterator(arr);
			while (size--) {
				u_long molId = arrayList_IteratorGetNext(iterator);

				fingerPrintIntType *dataFP = data_getFingerPrint(data, molId);
				double sim = data_tanimotoSimilarity(queryFP, dataFP, fpLen);

				minHeapInt_addIfKeyHigher(solutionHeap, sim, molId);
			}
			arrayList_IteratorFree(iterator);
			checkLength=1;
		}
	}


	//check the solution
	if (CHECK_SOLUTION_G) {
		static int c;
		printf ("%d comparing the 2 solutions\n",c++);
		if (c==40) {
			c=c+1-1;
		}
		minHeapInttype *solutionHeapLinear = newMinHeapInttype(k);
		for (int i=0; i<dSize; i++) {
			fingerPrintIntType *dataFP = data_getFingerPrint(data, i);
			double sim = data_tanimotoSimilarity (queryFP, dataFP, fpLen);

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
				//minHeapInt_print (solutionHeapLinear);
				//helper_err(__FILE__ , __LINE__, "Solution not matching");
				printf ("Solution not matching\n");

				printf ("************************************\n");
				double s1 = 0.08333;
				double s2 = 0.08334;
				printf ("len, molid, sim = %f, intersection, union \n",s1);
				for (int i=0; i<dSize; i++) {
					fingerPrintIntType *dataFP = data_getFingerPrint(data, i);
					double sim = data_tanimotoSimilarity (queryFP, dataFP, fpLen);
					int in = data_getIntersection(queryFP, dataFP, fpLen);
					int un = data_getUnion		 (queryFP, dataFP, fpLen);

					if ((s1 < sim) && (sim < s2)) {
						int len = data_getNum1Bits(dataFP, fpLen);
						printf ("%d : %d  %f %d %d\n", len, i, sim, in, un);
					}
				}
				printf ("************************************\n");
			}
		}



		minHeapInt_free(solutionHeapLinear);
	}

	 return solutionHeap;
}

workerFunctions_type_small AOR_getWorkerFunction () {
	workerFunctions_type_small worker;

	worker.dataDistribution = notImplemented;
	worker.init_index 	= _init_index_AOR;
	worker.free_index	= _free_index_AOR;
	worker.runRange		= _runRangeSmallHelper_AOR;
	worker.runTopK		= _runTopKSmallHelper_AOR;
	return worker;

}
