/*
 * AOR_BF.c
 *
 *  Created on: 26-Dec-2018
 *      Author: jithin
 */

#include "AOR_BF.h"

#define _ARRAY_BUCKET_SIZE_AOR_BF 1000

int _numBins_G = _G_AOR;

typedef struct {
	arrayListtype *molList;
	arrayListtype *referencePoint;
}_bin_AOR_BF;

typedef struct {
	int numBins;
	_bin_AOR_BF **bins;
}_bucket_AOR_BF;

typedef struct {
	int numFeatures; 		// the number of features in the index
	DataBinary *data;
	_bucket_AOR_BF **buckets;	//the array of index entries
}_index_AOR_BF;

/******************* Helpers *****************************/

/******************* Initializers ************************/
/************************ index ************************/
_bin_AOR_BF* _new_bin_AOR_BF () {
	_bin_AOR_BF *bin 	= (_bin_AOR_BF *)  helper_calloc(1,sizeof(_bin_AOR_BF),  __FILE__,__LINE__);
	bin->referencePoint = new_arrayListtype(100);
	bin->molList 		= new_arrayListtype(_ARRAY_BUCKET_SIZE_AOR_BF);

	return bin;
}

void _free_bin_AOR_BF (_bin_AOR_BF *bin) {
	if (!bin)
		return;
	free_arrayListtype(bin->molList);
	free_arrayListtype(bin->referencePoint);

	free (bin);
}

_bucket_AOR_BF* _new_bucket_AOR_BF (DataBinary *data, int len, int numBins) {
	_bucket_AOR_BF *bucket = NULL;
	int startIndex  	= data->_sizeIndexStart[len];
	if (startIndex >= 0) {
		int endIndex 	= data->_sizeIndexEnd[len];

		bucket 		= (_bucket_AOR_BF *)helper_calloc(1, 	  sizeof(_bucket_AOR_BF),__FILE__,__LINE__);
		bucket->bins= (_bin_AOR_BF   **)helper_calloc(numBins,sizeof(_bin_AOR_BF*),  __FILE__,__LINE__);

		int size = endIndex-startIndex+1;

		int n = ceil ((double)size/numBins);

		int fpIndex	=startIndex;
		int m 		= n;
		_bin_AOR_BF *bin = _new_bin_AOR_BF ();
		bucket->bins[bucket->numBins++] = bin;

		minHeapPlain_Type *minHeap = new_minHeapPlain_Type(100);
		do {
			arrayListtype *fingerPrint = dataBinary_getFingerPrint(data, fpIndex);
			arrayList_put(bin->molList, fpIndex);

			//update the reference point
			arrayList_Iterator_type *iterator = arrayList_getIterator(fingerPrint);
			for (int i=0; i<len; i++) {
				u_long feature = arrayList_IteratorGetNext(iterator);
				minHeapPlain_add_WO_fail_unique(&minHeap, feature);
			}
			arrayList_IteratorFree(iterator);

			fpIndex++;
			if (fpIndex >endIndex) {
				//update the reference point
				while (minHeap->size) {
					u_long feature;
					int count;

					minHeapPlain_popMultiple(minHeap, &feature, &count);
					arrayList_put(bin->referencePoint, feature);
				}

				break;
			}

			m--;
			if (m <= 0) {
				//update the reference point
				while (minHeap->size) {
					u_long feature;
					int count;

					minHeapPlain_popMultiple(minHeap, &feature, &count);
					arrayList_put(bin->referencePoint, feature);
				}

				m = n;
				bin = _new_bin_AOR_BF ();
				bucket->bins[bucket->numBins++] = bin;
			}
		}while (1);

		minHeapPlain_free(minHeap);
	}
	return bucket;
}

void _free_bucket_AOR_BF (_bucket_AOR_BF *bucket) {
	if (!bucket)
		return;

	for (int i=0; i< _numBins_G; i++) {
		_free_bin_AOR_BF(bucket->bins[i]);
	}
	free (bucket->bins);
	free (bucket);
}

void* _init_index_AOR_BF (DataBinary *data) {
	_index_AOR_BF *index = (_index_AOR_BF *) helper_calloc(1, sizeof(_index_AOR_BF),
			__FILE__, __LINE__);
	index->numFeatures 	= data->numFeatures;
	index->data			= data;

	int maxLen 			= data->maxNumFeatures;

	_bucket_AOR_BF **buckets = (_bucket_AOR_BF **)helper_malloc((maxLen)*sizeof(_bucket_AOR_BF *),
			__FILE__, __LINE__);

	buckets[0] = NULL;
	for (int len=1; len<=maxLen; len++) {
		buckets[len] = _new_bucket_AOR_BF (data, len, _numBins_G);
	}
	index->buckets = buckets;

	return index;
}

void _free_index_AOR_BF (void *source) {
	_index_AOR_BF *index = (_index_AOR_BF *)source;

	int maxLen = index->data->maxNumFeatures;
	for (int len=1; len<=maxLen; len++) {
		_free_bucket_AOR_BF(index->buckets[len]);
	}

	free (index->buckets);
	free (index);
}

/**********************************************/

/******************* Helpers *****************************/
arrayListtype* _runRangeHelper_AOR_BF(void *source, arrayListtype *queryFP , double r,
		double *pruned, double *unPruned) {
	_index_AOR_BF *index = (_index_AOR_BF *)source;

	DataBinary *data = index->data;
	int dSize = data->size;

	int L_q = queryFP->size;

	arrayListtype *solution = new_arrayListtype(_ARRAY_BUCKET_SIZE_AOR_BF);

	//get the feature lists

	int L_l = floor (r*L_q);
	int L_u = ceil ((double)L_q/r)+1;
	if (L_u > data->maxNumFeatures)
		L_u = data->maxNumFeatures;

	for (int l = L_l; l< L_u; l++){
		_bucket_AOR_BF *bucket = index->buckets[l];

		if (!bucket)
			continue;

		int B_plus_C = l + L_q;

		for (int i=0; i<bucket->numBins; i++) {
			_bin_AOR_BF   *bin 			= bucket->bins[i];
			arrayListtype *arr 			= bin->molList;
			arrayListtype *referenceFP  = bin->referencePoint;

			double gamma = dataBinary_getIntersection_fast(queryFP, referenceFP);
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

				arrayListtype *dataFP = dataBinary_getFingerPrint(data, molId);
				double sim = dataBinary_tanimotoSimilarity(queryFP, dataFP);

				if (sim >= r) {
					arrayList_put(solution, molId);
				}
			}
			arrayList_IteratorFree(iterator);
		}
	}

	if (CHECK_SOLUTION_G) {
		//int *featureIndex = (int *)helper_malloc((numFeatures+1)*sizeof(int),__FILE__, __LINE__);
		//verify the solution
		static int c;

		printf ("%3d)comparing the 2 solutions : ",c++);
		arrayListtype *solutionLinear = new_arrayListtype(_ARRAY_BUCKET_SIZE_AOR_BF);
		for (int i=0; i<dSize; i++) {
			arrayListtype *dataFP = dataBinary_getFingerPrint(data, i);
			double sim = dataBinary_tanimotoSimilarity (queryFP, dataFP);

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

				arrayListtype *dataFP;
				dataFP 		= dataBinary_getFingerPrint	 (data, a);
				int n1 		= dataFP->size;
				int inter_1 = dataBinary_getIntersection (queryFP, dataFP);
				dataFP 		= dataBinary_getFingerPrint	 (data, b);
				int n2 		= dataFP->size;
				int inter_2 = dataBinary_getIntersection (queryFP, dataFP);

				if (a != b) {
					printf ("Point 2 : ");
					dataBinary_printFeatures(dataFP);
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

minHeapInttype* _runTopKHelper_AOR_BF(void *source, arrayListtype *queryFP , int k,
		double *pruned, double *unPruned) {
	_index_AOR_BF *index = (_index_AOR_BF *)source;
	DataBinary *data = index->data;
	int dSize = data->size;

	int L_q = queryFP->size;

	double defaultKey = 1.0/(L_q+index->data->maxNumFeatures-1);
	minHeapInttype *solutionHeap = newMinHeapInttype_withDefault(k, defaultKey);

	double r = defaultKey;
	int L_l = floor (r*L_q);
	int L_u = ceil ((double)L_q/r)+1;
	if (L_u > data->maxNumFeatures)
		L_u = data->maxNumFeatures;

	int checkLength = 0;

	for (int l = L_q; l>0; l--){
		_bucket_AOR_BF *bucket = index->buckets[l];

		if (!bucket)
			continue;

		int B_plus_C = l + L_q;

		for (int i=0; i<bucket->numBins; i++) {
			if (checkLength) {
				r = minHeapInt_peekKey_or_default(solutionHeap);
				checkLength = 0;

				L_l = floor (r*L_q);
				if (l < L_l) {
					l = 0; // we can stop completely
					break;
				}
			}

			_bin_AOR_BF *bin 			= bucket->bins[i];
			arrayListtype *referenceFP 	= bin->referencePoint;

			double gamma = dataBinary_getIntersection_fast(queryFP, referenceFP);

			gamma = gamma*(r+1);
			if (B_plus_C < gamma-0.2)
				continue; //prune
			gamma = gamma/r;

			if (B_plus_C > (gamma+0.2)) //floating point error
				continue; //prune

			arrayListtype *arr 	= bin->molList;
			int size 			= arr->size;
			arrayList_Iterator_type *iterator = arrayList_getIterator(arr);
			while (size--) {
				u_long molId = arrayList_IteratorGetNext(iterator);

				arrayListtype *dataFP = dataBinary_getFingerPrint(data, molId);
				double sim = dataBinary_tanimotoSimilarity(queryFP, dataFP);

				minHeapInt_addIfKeyHigher(solutionHeap, sim, molId);
			}
			checkLength=1;
		}
	}

	for (int l = L_q+1; l< L_u; l++){
		_bucket_AOR_BF *bucket = index->buckets[l];

		if (!bucket)
			continue;

		int B_plus_C = l + L_q;


		for (int i=0; i<bucket->numBins; i++) {
			if (checkLength) {
				r = minHeapInt_peekKey_or_default(solutionHeap);
				checkLength = 0;
				L_u = ceil ((double)L_q/r)+1;
				if (L_u > data->maxNumFeatures)
					L_u = data->maxNumFeatures;

				if (l >= L_u)
					break;
			}

			_bin_AOR_BF *bin 			= bucket->bins[i];
			arrayListtype *referenceFP 	= bin->referencePoint;
			double gamma = dataBinary_getIntersection_fast(queryFP, referenceFP);

			gamma = gamma*(r+1);
			if (B_plus_C < gamma-0.2)
				continue; //prune
			gamma = gamma/r;

			if (B_plus_C > (gamma+0.2)) //floating point error
				continue; //prune

			arrayListtype *arr 					= bin->molList;
			int size 							= arr->size;
			arrayList_Iterator_type *iterator 	= arrayList_getIterator(arr);
			while (size--) {
				u_long molId = arrayList_IteratorGetNext(iterator);

				arrayListtype *dataFP = dataBinary_getFingerPrint(data, molId);
				double sim = dataBinary_tanimotoSimilarity(queryFP, dataFP);

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
		minHeapInttype *solutionHeapLinear = newMinHeapInttype(k);
		for (int i=0; i<dSize; i++) {
			arrayListtype *dataFP = dataBinary_getFingerPrint(data, i);
			double sim = dataBinary_tanimotoSimilarity (queryFP, dataFP);

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
				printf ("Solution not matching\n");
				printf ("Len, molid, sim\n");
				{
					_intheap_element *arr = solutionHeapLinear->heap;

					for (int i=0; i<k; i++) {
						int molid = arr[i].val;
						arrayListtype *dataFP = dataBinary_getFingerPrint(data, molid);
						int len = dataFP->size;
						printf ("%d, %d, %f\n",len, molid, arr[i].key);
					}
				}
				helper_err(__FILE__ , __LINE__, "Solution not matching");

				/*
				printf ("************************************\n");
				double s1 = 0.08333;
				double s2 = 0.08334;
				printf ("len, molid, sim = %f, intersection, union \n",s1);
				for (int i=0; i<dSize; i++) {
					arrayListtype *dataFP = dataBinary_getFingerPrint(data, i);
					double sim 	= dataBinary_tanimotoSimilarity	(queryFP, dataFP);
					int in 		= dataBinary_getIntersection	(queryFP, dataFP);
					int un 		= dataBinary_getUnion		 	(queryFP, dataFP);

					if ((s1 < sim) && (sim < s2)) {
						int len = dataFP->size;
						printf ("%d : %d  %f %d %d\n", len, i, sim, in, un);
					}
				}
				printf ("************************************\n");
				*/
			}
		}

		minHeapInt_free(solutionHeapLinear);
	}

	return solutionHeap;
}

workerFunctions_type AOR_BF_getWorkerFunction () {
	workerFunctions_type worker;

	worker.dataDistribution = notImplemented;
	worker.init_index 	= _init_index_AOR_BF;
	worker.free_index	= _free_index_AOR_BF;
	worker.runRange		= _runRangeHelper_AOR_BF;
	worker.runTopK		= _runTopKHelper_AOR_BF;
	return worker;

}
