/*
 * divideSkip.c
 *
 *  Created on: 10-Dec-2018
 *      Author: jithin
 */

#include "divideSkip.h"

#define _ARRAY_BUCKET_SIZE_DS 100

const double mu = 0.0085;

typedef struct {
	arrayListtype_flexible **featureList;
	int maxFeatureListLen;
}_bucket_DS;

typedef struct {
	int numFeatures; 		// the number of features in the index
	DataSmall *data;
	_bucket_DS **buckets;	//the array of index entries
}_index_DS;

/******************* Initializers ************************/
/************************ index ************************/
_bucket_DS*_new_bucket_DS (DataSmall *data, int len, int numBits, int *featureIndex) {
	_bucket_DS *bucket = NULL;
	int startIndex  	= data->_sizeIndexStart[len];
	if (startIndex >= 0) {
		int endIndex 	= data->_sizeIndexEnd[len];

		bucket = (_bucket_DS *) helper_malloc (sizeof(_bucket_DS),
				__FILE__, __LINE__);
		bucket->featureList = (arrayListtype_flexible **) helper_calloc (numBits,
				sizeof(arrayListtype_flexible *), __FILE__, __LINE__);

		for (int fpIndex=startIndex; fpIndex <=endIndex; fpIndex++) {
			fingerPrintIntType *fingerPrint = data_getFingerPrint(data, fpIndex);
			int num1Bits =  data_getIndexOf1Bits (featureIndex, fingerPrint, data->fingerPrintLen);

			for (int f=0; f<num1Bits; f++) {
				int feature = featureIndex[f];

				if (bucket->featureList[feature] == NULL)
					bucket->featureList[feature] = new_arrayListtype_flexible(_ARRAY_BUCKET_SIZE_DS, data->size);

				arrayList_flexible_put(&(bucket->featureList[feature]), fpIndex);
			}
		}
		int maxFeatureListLen=0;
		for (int i=0; i<numBits; i++) {
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

void _free_bucket_DS (_bucket_DS *bucket, int numBits) {
	if (bucket) {
		for (int f=0; f<numBits; f++) {
			if (bucket->featureList[f])
				free_arrayListtype_flexible(bucket->featureList[f]);
		}

		free (bucket->featureList);
		free (bucket);
	}
}

void* _init_index_DS (DataSmall *data) {
	_index_DS *index = (_index_DS *)helper_calloc(1, sizeof(_index_DS),__FILE__, __LINE__);
	int numBits = data->fingerPrintLen * NUMBITS_fingerPrintIntType;
	index->numFeatures 	= numBits;
	index->data			= data;
	int *features 		= (int *)helper_malloc((numBits+1)*sizeof(int),__FILE__, __LINE__);

	_bucket_DS **buckets = (_bucket_DS **) helper_malloc ((numBits+1) * sizeof(_bucket_DS **),
			__FILE__, __LINE__);

	for (int len=1; len<=numBits; len++) {
		buckets[len] = _new_bucket_DS (data, len, numBits, features);
	}
	index->buckets = buckets;
	free (features);
	return index;
}

void _free_index_DS (void *source) {
	_index_DS *index = (_index_DS *) source;
	int numBits = index->numFeatures;

	for (int len=1; len<=numBits; len++) {
		_free_bucket_DS (index->buckets[len], numBits);
	}

	free (index->buckets);
	free (index);
}

/**********************************************/

/******************* Helpers *****************************/
/**
 * sort features in ascending order size of mList
 */
void _sortFeatures_DS (int *featureIndex, int num1Bits, _bucket_DS *bucket) {
	int *size = (int *)helper_malloc((num1Bits)*sizeof(int),__FILE__, __LINE__);

	for (int i=0; i<num1Bits; i++) {
		int feature 		= featureIndex[i];
		arrayListtype_flexible *arr	= bucket->featureList[feature];
		if (arr)
			size[i] = arr->size;
		else
			size[i] = 0;
	}

	//sort
	for (int i=0; i<num1Bits; i++) {
		for (int k=i+1; k<num1Bits; k++) {
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

	free (size);
}

arrayListtype* _runRangeSmallHelper_DS(void *source, fingerPrintIntType *queryFP , double r,
		double *pruned, double *unPruned) {
	_index_DS *index = (_index_DS *)source;
	DataSmall *data = index->data;
	int dSize = data->size;
	int fpLen = data->fingerPrintLen;
	int numFeatures = index->numFeatures;

	int *featureIndex = (int *)helper_malloc((numFeatures+1)*sizeof(int),__FILE__, __LINE__);
	int num1Bits = data_getIndexOf1Bits (featureIndex, queryFP, fpLen);

	int L_q = num1Bits;

	arrayListtype *solution = new_arrayListtype(_ARRAY_BUCKET_SIZE_DS);

	//get the feature lists
	arrayList_flexible_Iterator_type *iterators	[num1Bits];

	int L_l = floor (r*L_q);
	int L_u = ceil ((double)L_q/r)+1;
	if (L_u > numFeatures+1)
		L_u = numFeatures+1;

	int *processed = (int *) helper_calloc(dSize, sizeof(int), __FILE__, __LINE__);

	minHeapPlain_Type *heap = new_minHeapPlain_Type(numFeatures*10000);

	for (int l=L_l; l<L_u; l++) {
		_bucket_DS *bucket = index->buckets[l];

		if (!bucket)
			continue;

		_sortFeatures_DS (featureIndex, num1Bits, bucket);

		int numFeaturesActive = num1Bits;

		double T = r*(L_q+l)/(r+1)-0.01;
		int  L = floor (T/(mu  * log(bucket->maxFeatureListLen) + 1));

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
		arrayListtype *solutionLinear = new_arrayListtype(_ARRAY_BUCKET_SIZE_DS);
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
				int a = arrayList_IteratorGetNext(iteratorI);
				int b = arrayList_IteratorGetNext(iteratorL);



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
	free (featureIndex);
	return (solution);
}

void _TopKSmallHelperHelper_DS (int l, _bucket_DS *bucket, int *featureIndex,
		int num1Bits, int L_q, minHeapInttype *solutionHeap,
		arrayList_flexible_Iterator_type *iterators[], minHeapPlain_Type *heap,
		DataSmall *data, fingerPrintIntType *queryFP, int fpLen, double r) {
	_sortFeatures_DS (featureIndex, num1Bits, bucket);

	int numFeaturesActive = num1Bits;

	double T = r*(L_q+l)/(r+1)-1;
	if (T < 0)
		T = 0;
	int  L = floor (T/(mu  * log(bucket->maxFeatureListLen) + 1));

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
				fingerPrintIntType *dataFP = data_getFingerPrint(data, molId);
				double intersection = data_getIntersection(queryFP, dataFP, fpLen);
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
					fingerPrintIntType *dataFP = data_getFingerPrint(data, molId);
					double intersection = data_getIntersection(queryFP, dataFP, fpLen);
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

minHeapInttype* _runTopKSmallHelper_DS(void *source, fingerPrintIntType *queryFP , int k,
		double *pruned, double *unPruned) {
	_index_DS *index = (_index_DS *)source;
	DataSmall *data = index->data;
	int dSize = data->size;
	int fpLen = data->fingerPrintLen;
	int numFeatures = index->numFeatures;

	int *featureIndex = (int *)helper_malloc((index->numFeatures+1)*sizeof(int),__FILE__, __LINE__);
	int num1Bits = data_getIndexOf1Bits (featureIndex, queryFP, fpLen);

	int L_q = num1Bits;
	minHeapInttype *solutionHeap = newMinHeapInttype(k);

	double r = (double)1.0/numFeatures;

	//get the feature lists
	arrayList_flexible_Iterator_type *iterators	[num1Bits];
	minHeapPlain_Type *heap = new_minHeapPlain_Type(numFeatures*10000);

	int exploreRight=1,exploreLeft=1;
	int lLeft =L_q;
	int lRight=L_q+1;

	while (exploreLeft || exploreRight) {
		if (exploreLeft) {
			if (lLeft >= 0) {
				if (index->buckets[lLeft]) {
					r = minHeapInt_peekKey(solutionHeap);

					if (r == 0)
						r=(double)1.0/numFeatures;
					if (lLeft < (r*L_q)) {
						exploreLeft = 0;
					} else {
						_TopKSmallHelperHelper_DS (lLeft, index->buckets[lLeft], featureIndex,
								num1Bits, L_q, solutionHeap, iterators, heap, data,
								queryFP, fpLen, r);
					}
				}
				lLeft--;
			} else {
				exploreLeft = 0;
			}
		}

		if (exploreRight) {
			if(lRight<= numFeatures) {
				if (index->buckets[lRight]) {
					r = minHeapInt_peekKey(solutionHeap);

					if (r == 0)
						r=(double)1.0/numFeatures;
					if (lRight > ((double)L_q/r+1)) {
						exploreRight = 0;
					} else {
						_TopKSmallHelperHelper_DS (lRight, index->buckets[lRight], featureIndex,
								num1Bits, L_q, solutionHeap, iterators, heap, data,
								queryFP, fpLen, r);
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
				minHeapInt_print (solutionHeapLinear);
				//helper_err(__FILE__ , __LINE__, "Solution not matching");
				printf ("Solution not matching\n");
			}
		}

		minHeapInt_free(solutionHeapLinear);
	}

	free (featureIndex);

	return solutionHeap;
}

workerFunctions_type_small divideSkip_getWorkerFunction () {
	workerFunctions_type_small worker;

	worker.dataDistribution = notImplemented;
	worker.init_index 	= _init_index_DS;
	worker.free_index	= _free_index_DS;
	worker.runRange		= _runRangeSmallHelper_DS;
	worker.runTopK		= _runTopKSmallHelper_DS;
	return worker;

}
