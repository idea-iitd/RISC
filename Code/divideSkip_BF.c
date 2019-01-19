/*
 * divideSkip_BF.c
 *
 *  Created on: 25-Dec-2018
 *      Author: jithin
 */

#include "divideSkip_BF.h"

#define _ARRAY_BUCKET_SIZE_DS_BF 100

const double mu_BF = 0.0085;

typedef struct {
	arrayListtype_flexible **featureList;
	int maxFeatureListLen;
}_bucket_DS_BF;

typedef struct {
	int numFeatures; 		// the number of features in the index
	int maxNumOfFeatures;   // max number of features in any molecule
	DataBinary *data;
	_bucket_DS_BF **buckets;	//the array of index entries
}_index_DS_BF;

/******************* Initializers ************************/
/************************ index ************************/
_bucket_DS_BF*_new_bucket_DS_BF (DataBinary *data, int len, int numFeatures) {
	_bucket_DS_BF *bucket = NULL;
	static int progress;
	int startIndex  	= data->_sizeIndexStart[len];
	if (startIndex >= 0) {
		int endIndex 	= data->_sizeIndexEnd[len];
		helper_printProgress(progress++);

		bucket = (_bucket_DS_BF *) helper_malloc (sizeof(_bucket_DS_BF),
				__FILE__, __LINE__);
		bucket->featureList = (arrayListtype_flexible **) helper_calloc (numFeatures,
				sizeof(arrayListtype_flexible *), __FILE__, __LINE__);

		for (int fpIndex=startIndex; fpIndex <=endIndex; fpIndex++) {
			arrayListtype *fingerPrint 		= dataBinary_getFingerPrint(data, fpIndex);
			arrayList_Iterator_type *iter	= arrayList_getIterator(fingerPrint);

			for (int f=0; f<len; f++) {
				int feature = arrayList_IteratorGetNext(iter);

				if (bucket->featureList[feature] == NULL)
					bucket->featureList[feature] = new_arrayListtype_flexible(_ARRAY_BUCKET_SIZE_DS_BF, data->size);

				arrayList_flexible_put(&(bucket->featureList[feature]), fpIndex);
			}
			arrayList_IteratorFree(iter);
		}
		int maxFeatureListLen=0;
		for (int i=0; i<numFeatures; i++) {
			arrayListtype_flexible *arr = bucket->featureList[i];

			if (arr) {
				int l = arr->size;
				if (l > maxFeatureListLen)
					maxFeatureListLen = l;
			}
		}

		bucket->maxFeatureListLen = maxFeatureListLen;
	}
	return bucket;
}

void _free_bucket_DS_BF (_bucket_DS_BF *bucket, int numFeatures) {
	if (!bucket)
		return;

	for (int i=0; i<numFeatures; i++) {
		arrayListtype_flexible *arr = bucket->featureList[i];

		if (arr) {
			free_arrayListtype_flexible(arr);
		}
	}

	free (bucket->featureList);
	free (bucket);
}

void* _init_index_DS_BF (DataBinary *data) {
	_index_DS_BF *index = (_index_DS_BF *) helper_calloc (1, sizeof(_index_DS_BF),
			__FILE__, __LINE__);
	index->numFeatures 	= data->numFeatures;
	index->data			= data;

	int maxLen = data->maxNumFeatures;
	index->maxNumOfFeatures = maxLen;

	_bucket_DS_BF **buckets = (_bucket_DS_BF **) helper_malloc ((maxLen+1) * sizeof(_bucket_DS_BF *),
			__FILE__, __LINE__);

	for (int len=1; len<=maxLen; len++) {
		buckets[len] = _new_bucket_DS_BF (data, len, data->numFeatures);
	}
	index->buckets = buckets;

	return index;
}

void _free_index_DS_BF (void *source) {
	_index_DS_BF *index = (_index_DS_BF *) source;

	int maxLen = index->data->maxNumFeatures;
	for (int len=1; len<=maxLen; len++) {
		_free_bucket_DS_BF (index->buckets[len], index->data->numFeatures);
	}
	free (index->buckets);
	free (index);
}

/**********************************************/
/******************* Helpers *****************************/
/**
 * sort features in ascending order size of mList
 */
void _sortFeatures_DS_BF (int *featureIndex, arrayListtype *queryFP, _bucket_DS_BF *bucket) {
	int L_q = queryFP->size;
	int size[L_q];

	arrayList_Iterator_type *iterator = arrayList_getIterator(queryFP);
	for (int i=0; i<L_q; i++) {
		int feature 		= arrayList_IteratorGetNext(iterator);
		arrayListtype_flexible *arr	= bucket->featureList[feature];
		featureIndex[i] 	= feature;

		if (arr)
			size[i] = arr->size;
		else
			size[i] = 0;
	}
	arrayList_IteratorFree(iterator);

	//sort
	for (int i=0; i<L_q; i++) {
		for (int k=i+1; k<L_q; k++) {
			if (size[i] > size[k]) {
				//swap len
				{
					int temp= size[i];
					size[i]	= size[k];
					size[k]	= temp;
				}

				//swap the featurelist
				{
					int temp = featureIndex[i];
					featureIndex[i] = featureIndex[k];
					featureIndex[k] = temp;
				}
			}
		}
	}
}

arrayListtype* _runRangeHelper_DS(void *source, arrayListtype *queryFP , double r,
		double *pruned, double *unPruned) {
	_index_DS_BF *index = (_index_DS_BF *) source;

	DataBinary *data = index->data;
	int dSize = data->size;
	int maxNumFeatures = index->maxNumOfFeatures;

	int L_q = queryFP->size;

	arrayListtype *solution = new_arrayListtype(_ARRAY_BUCKET_SIZE_DS_BF);
	int featureIndex[maxNumFeatures +1];

	//get the feature lists
	arrayList_flexible_Iterator_type *iterators	[maxNumFeatures];

	int L_l = floor (r*L_q);
	int L_u = ceil ((double)L_q/r)+1;
	if (L_u > maxNumFeatures+1)
		L_u = maxNumFeatures+1;

	int *processed = (int *) helper_calloc(dSize, sizeof(int), __FILE__, __LINE__);
	minHeapPlain_Type *heap = new_minHeapPlain_Type(maxNumFeatures*100000);

	for (int l=L_l; l<L_u; l++) {
		_bucket_DS_BF *bucket = index->buckets[l];

		if (!bucket)
			continue;

		_sortFeatures_DS_BF (featureIndex, queryFP, bucket);

		int numFeaturesActive = L_q;

		double T = r*(L_q+l)/(r+1)-0.01;
		int  L = floor (T/(mu_BF  * log(bucket->maxFeatureListLen) + 1));

		int l_Border= numFeaturesActive - L;
		if (l_Border < 0)
			l_Border = 0;
		int lShort 	= l_Border;
		int lLong 	= numFeaturesActive;

		for (int i=0,j=0; i<lShort; i++){
			int f = featureIndex[j++];
			arrayListtype_flexible *arr	= bucket->featureList[f];
			if (arr) {
				iterators[i] = arrayList_flexible_getIterator(arr);
			} else {
				lShort--;
				iterators[lShort] = NULL;
				i--;
			}
		}

		for (int i=l_Border,j=l_Border; i<lLong; i++){
			int f = featureIndex[j++];
			arrayListtype_flexible *arr	= bucket->featureList[f];
			if (arr) {
				iterators[i] = arrayList_flexible_getIterator(arr);
			} else {
				lLong--;
				iterators[lLong] = NULL;
				i--;
			}
		}
		numFeaturesActive = lLong;

		//initialize heap
		for (int i=0; i<lShort; i++) {
			arrayList_flexible_Iterator_type *iterator = iterators[i];

			if (!iterator->nextAvailable) {
				//the iterator has been exhausted
				lShort--;
				iterators[i] = iterators[lShort];
				i--;
			} else {
				u_long molId = arrayList_flexible_IteratorGetNext(iterator);
				//populate the minheap
				minHeapPlain_add(heap, molId,0);
			}
		}

		while (heap->size > 0) {
			u_long molId;
			int count;

			minHeapPlain_popMultiple(heap, &molId, &count);

			int n = count;
				if (count >= T-L) {
					if (count >= T) {
						arrayList_put(solution, molId);
					} else {
						//count the number of occurance in long list
						for (int i=l_Border; i<lLong; i++) {
							arrayList_flexible_Iterator_type *iterator = iterators[i];

							int found;
							int notEnoughElements = arrayList_flexible_IteratorFind(iterator, molId, &found);
							if (notEnoughElements) {
								//the iterator has been exhausted
								lLong--;
								iterators[i] = iterators[lLong];
								iterators[lLong] = iterator;
								i--;
							} else {
								if (found)
									count++;
							}
						}

						if (count >= T)
							arrayList_put(solution, molId);
					}

					for (int i=0; i<lShort; i++) {
						arrayList_flexible_Iterator_type *iterator = iterators[i];

						if (!iterator->nextAvailable) {
							//the iterator has been exhausted
							lShort--;
							iterators[i] = iterators[lShort];
							iterators[lShort] = iterator;
							i--;
						} else {
							molId = arrayList_flexible_IteratorGetNext(iterator);
							//populate the minheap
							minHeapPlain_add(heap, molId,0);
						}
					}
				}

				do {
					int popCount = T-L-1-n;
					popCount = 0;
					for (;popCount>0;popCount--) {
						minHeapPlain_pop(heap, &molId);
					}
					u_long nextMolId;
					minHeapPlain_peek(heap, &nextMolId);

					for (int i=0; i<lShort; i++) {
						arrayList_flexible_Iterator_type *iterator = iterators[i];

						arrayList_flexible_IteratorMoveTo_Binary(iterator, nextMolId);
						if (!iterator->nextAvailable) {
							//the iterator has been exhausted
							lShort--;
							iterators[i] = iterators[lShort];
							iterators[lShort] = iterator;
							i--;
						} else {
							int molId_1 = arrayList_flexible_IteratorGetNext(iterator);
							//populate the minheap
							minHeapPlain_add(heap, molId_1,0);
						}
					}
				} while (0);
			}
			for (int i=0; i<numFeaturesActive; i++) {
				if (iterators[i])
					arrayList_flexible_IteratorFree(iterators[i]);
			}
	}

	if (CHECK_SOLUTION_G) {
		//verify the solution
		static int c;

		printf ("%3d)comparing the 2 solutions : ",c++);
		arrayListtype *solutionLinear = new_arrayListtype(_ARRAY_BUCKET_SIZE_DS_BF);
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
				int a = arrayList_IteratorGetNext(iteratorI);
				int b = arrayList_IteratorGetNext(iteratorL);



				arrayListtype *dataFP;
				dataFP 		= dataBinary_getFingerPrint	(data, a);
				int n1 		= dataFP->size;
				int inter_1 = dataBinary_getIntersection	(queryFP, dataFP);
				dataFP 		= dataBinary_getFingerPrint	(data, b);
				int n2 		= dataFP->size;
				int inter_2 = dataBinary_getIntersection	(queryFP, dataFP);

				if (a != b) {
					printf ("Point 2 : ");
					for (int j=0; j<n2; j++) {
						printf ("%d,",featureIndex[j]);
					}
					printf ("\n");
					printf ("id:len,intersection; \n");
					printf ("%d:%d,%d; %d:%d,%d; \n",a,n1,inter_1,b,n2,inter_2);
					helper_err(__FILE__, __LINE__, "Solutions do not match");
				}
			}
			printf("\n");
			arrayList_IteratorFree(iteratorI);
			arrayList_IteratorFree(iteratorL);
		}

		free_arrayListtype(solutionLinear);
	}

	free (processed);

	return solution;
}

void _TopKHelperHelper_DS (int l, _bucket_DS_BF *bucket, int *featureIndex,
		int L_q, minHeapInttype *solutionHeap,
		arrayList_flexible_Iterator_type *iterators[], minHeapPlain_Type *heap,
		DataBinary *data, arrayListtype *queryFP, double r) {

	_sortFeatures_DS_BF (featureIndex, queryFP, bucket);

	int numFeaturesActive = L_q;

	double T = r*(L_q+l)/(r+1)-1;
	if (T < 0)
		T = 0;
	int  L = floor (T/(mu_BF  * log(bucket->maxFeatureListLen) + 1));

	int l_Border= numFeaturesActive - L;
	if (l_Border < 0)
		l_Border = 0;
	int lShort 	= l_Border;
	int lLong 	= numFeaturesActive;

	for (int i=0,j=0; i<lShort; i++){
		int f = featureIndex[j++];

		arrayListtype_flexible *arr	= bucket->featureList[f];
		if (arr) {
			iterators[i] = arrayList_flexible_getIterator(arr);
		} else {
			lShort--;
			iterators[lShort] = NULL;
			i--;
		}
	}

	for (int i=l_Border,j=l_Border; i<lLong; i++){
		int f = featureIndex[j++];

		arrayListtype_flexible *arr	= bucket->featureList[f];
		if (arr) {
			iterators[i] = arrayList_flexible_getIterator(arr);
		} else {
			lLong--;
			iterators[lLong] = NULL;
			i--;
		}
	}
	numFeaturesActive = lLong;

	//initialize heap
	for (int i=0; i<lShort; i++) {
		arrayList_flexible_Iterator_type *iterator = iterators[i];

		if (!iterator->nextAvailable) {
			//the iterator has been exhausted
			lShort--;
			iterators[i] = iterators[lShort];
			i--;
		} else {
			int molId = arrayList_flexible_IteratorGetNext(iterator);
			//populate the minheap
			minHeapPlain_add(heap, molId,0);
		}
	}

	int tAdjuster = 0;
	while (heap->size > 0) {
		u_long molId;
		int count;

		minHeapPlain_popMultiple(heap, &molId, &count);

		int n = count;
		if (count >= T-L-tAdjuster) {
			if (count >= T-tAdjuster) {
				arrayListtype *dataFP = dataBinary_getFingerPrint(data, molId);
				double intersection = dataBinary_getIntersection_fast(queryFP, dataFP);
				double sim = 0;
				int u = l+L_q-intersection;
				sim = intersection/u;
				minHeapInt_addIfKeyHigher(solutionHeap, sim, molId);
			} else {
				//count the number of occurance in long list
				for (int i=l_Border; i<lLong; i++) {
					arrayList_flexible_Iterator_type *iterator = iterators[i];

					int found;
					int notEnoughElements = arrayList_flexible_IteratorFind(iterator, molId, &found);
					if (notEnoughElements) {
						//the iterator has been exhausted
						lLong--;
						iterators[i] = iterators[lLong];
						iterators[lLong] = iterator;
						i--;
					} else {
						if (found)
							count++;
					}
				}

				if (count >= T-tAdjuster) {
					arrayListtype *dataFP = dataBinary_getFingerPrint(data, molId);
					double intersection = dataBinary_getIntersection_fast(queryFP, dataFP);
					double sim = 0;
					int u = l+L_q-intersection;
					sim = intersection/u;
					minHeapInt_addIfKeyHigher(solutionHeap, sim, molId);
				}
			}

			for (int i=0; i<lShort; i++) {
				arrayList_flexible_Iterator_type *iterator = iterators[i];

				if (!iterator->nextAvailable) {
					//the iterator has been exhausted
					lShort--;
					iterators[i] = iterators[lShort];
					iterators[lShort] = iterator;
					i--;
				} else {
					molId = arrayList_flexible_IteratorGetNext(iterator);
					//populate the minheap
					minHeapPlain_add(heap, molId,0);
				}
			}
		}

		do {
			int popCount = T-L-1-n - tAdjuster;
			popCount = 0;
			for (;popCount>0;popCount--) {
				minHeapPlain_pop(heap, &molId);
			}
			u_long nextMolId;
			minHeapPlain_peek(heap, &nextMolId);

			for (int i=0; i<lShort; i++) {
				arrayList_flexible_Iterator_type *iterator = iterators[i];

				arrayList_flexible_IteratorMoveTo_Binary(iterator, nextMolId);
				if (!iterator->nextAvailable) {
					//the iterator has been exhausted
					lShort--;
					iterators[i] = iterators[lShort];
					iterators[lShort] = iterator;
					i--;
				} else {
					int molId_1 = arrayList_flexible_IteratorGetNext(iterator);
					//populate the minheap
					minHeapPlain_add(heap, molId_1,0);
				}
			}
		} while (0);
	}
	for (int i=0; i<numFeaturesActive; i++) {
		if (iterators[i])
			arrayList_flexible_IteratorFree(iterators[i]);
	}
}

minHeapInttype* _runTopKHelper_DS(void *source, arrayListtype *queryFP , int k,
		double *pruned, double *unPruned) {
	_index_DS_BF *index = (_index_DS_BF *) source;

	DataBinary *data = index->data;
	int maxNumFeatures = index->maxNumOfFeatures;

	int L_q = queryFP->size;
	int featureIndex[L_q];

	double defaultKey = 1.0/(L_q+index->data->maxNumFeatures-1);
	minHeapInttype *solutionHeap = newMinHeapInttype_withDefault(k, defaultKey);

	//get the feature lists
	arrayList_flexible_Iterator_type *iterators	[maxNumFeatures];
	minHeapPlain_Type *heap = new_minHeapPlain_Type(maxNumFeatures*100000);

	int exploreRight=1,exploreLeft=1;
	int lLeft =L_q;
	int lRight=L_q+1;

	while (exploreLeft || exploreRight) {
		if (exploreLeft) {
			if (lLeft >= 0) {
				if (index->buckets[lLeft]) {
					double r = minHeapInt_peekKey_or_default(solutionHeap);

					if (lLeft < (r*L_q)) {
						exploreLeft = 0;
					} else {
						_TopKHelperHelper_DS (lLeft, index->buckets[lLeft], featureIndex,
								L_q, solutionHeap, iterators, heap, data, queryFP, r);
					}
				}
				lLeft--;
			} else {
				exploreLeft = 0;
			}
		}

		if (exploreRight) {
			if(lRight<= maxNumFeatures) {
				if (index->buckets[lRight]) {
					double r = minHeapInt_peekKey_or_default(solutionHeap);

					if (lRight > ((double)L_q/r+1)) {
						exploreRight = 0;
					} else {
						_TopKHelperHelper_DS (lRight, index->buckets[lRight], featureIndex,
								L_q, solutionHeap, iterators, heap, data, queryFP, r);
					}
				}
				lRight++;
			} else {
				exploreRight = 0;
			}
		}
	}

	//check the solution
	if (CHECK_SOLUTION_G) {
		static int c;
		printf ("%d comparing the 2 solutions\n",c++);
		minHeapInttype *solutionHeapLinear = newMinHeapInttype(k);
		for (int i=0; i<data->size; i++) {
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
				minHeapInt_print (solutionHeapLinear);
				//helper_err(__FILE__ , __LINE__, "Solution not matching");
				printf ("Solution not matching\n");
			}
		}

		minHeapInt_free(solutionHeapLinear);
	}

	return solutionHeap;
}

workerFunctions_type divideSkip_BF_getWorkerFunction () {
	workerFunctions_type worker;

	worker.dataDistribution = notImplemented;
	worker.init_index 	= _init_index_DS_BF;
	worker.free_index	= _free_index_DS_BF;
	worker.runRange		= _runRangeHelper_DS;
	worker.runTopK		= _runTopKHelper_DS;
	return worker;

}
