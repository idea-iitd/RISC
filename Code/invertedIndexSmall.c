/*
 * invertedindex.c
 *
 *  Created on: 27-Nov-2018
 *      Author: jithin
 */

#include "invertedIndexSmall.h"

#define _ARRAY_BUCKET_SIZE_II 1000
//#define _STAT_PRUNING
//#define _USE_TANIMOTO

/******************* Typedefs ****************************/
typedef struct {
	arrayListtype_flexible *arr;
	int minNum1Bits;				// the size of the smallest mol
	int maxNum1Bits;				// the size of the largest mol
	int *_sizeIndexStart;			// stores the starting index of molecules of given size
	int *_sizeIndexEnd;				// stores the ending index(inclusive) of molecules of given size
}_moleculeList;

typedef struct {
	int numFeatures; 					// the number of features in the index
	DataSmall *data;
	_moleculeList **moleculeList;	//the array of index entries
}_index;


/******************* Helpers *****************************/
void _addFPIndex2moleculeList 	(_moleculeList *mList, int fpIndex, int num1Bits);

/******************* Initializers ************************/
_moleculeList* _new_moleculeList (int numBits, int dSize) {
	_moleculeList *mList = (_moleculeList *) helper_calloc(1, sizeof(_moleculeList),
			__FILE__, __LINE__);

	mList->_sizeIndexStart	= (int *) helper_malloc(sizeof(int) * (numBits + 1),
			__FILE__, __LINE__);
	mList->_sizeIndexEnd	= (int *) helper_malloc(sizeof(int) * (numBits + 1),
			__FILE__, __LINE__);

	for (int i=0; i<=numBits; i++) {
		mList->_sizeIndexStart[i]=-1;
	}

	mList->arr = new_arrayListtype_flexible(_ARRAY_BUCKET_SIZE_II, dSize);
	mList->minNum1Bits = INT_MAX;
	mList->maxNum1Bits = 0;

	return mList;
}

void _free_moleculeList (_moleculeList *mList) {

	free_arrayListtype_flexible(mList->arr);
	free (mList->_sizeIndexStart);
	free (mList->_sizeIndexEnd);
	free (mList);
}

/************************ index ************************/
void _addFingerPrint2index (fingerPrintIntType *fp, _index *index,
		int fpLen, int *featureIndex, int fpIndex) {
	int num1Bits = data_getIndexOf1Bits (featureIndex, fp, fpLen);

	for (int i=0; i<num1Bits; i++) {
		int feature = featureIndex[i];

		if (feature >= index->numFeatures)
			helper_err(__FILE__, __LINE__, "the feature value is more than the possible");

		_moleculeList *mList = index->moleculeList[feature];
		_addFPIndex2moleculeList (mList, fpIndex, num1Bits);
	}
}

void* _init_index (DataSmall *data) {
	_index *index = (_index *)helper_calloc(1, sizeof(_index),__FILE__,__LINE__);
	//allocate memory
	{
		int numBits = data->fingerPrintLen * NUMBITS_fingerPrintIntType;
		index->numFeatures 	= numBits;
		index->data			= data;

		_moleculeList **mList = (_moleculeList **) helper_malloc ((size_t)numBits * sizeof(_moleculeList **),
				__FILE__, __LINE__);

		for (int i=0; i<numBits; i++) {
			mList[i] = _new_moleculeList (numBits, data->size);
		}

		index->moleculeList = mList;
	}

	// load data
	int size 		= data->size;
	int fpLen		= data->fingerPrintLen;
	int numBits 	= fpLen*NUMBITS_fingerPrintIntType;
	int *features 	= (int *)helper_malloc((numBits+1)*sizeof(int),__FILE__, __LINE__);

	for (int i=0; i<size; i++) {
		fingerPrintIntType *fingerPrint = data_getFingerPrint(data, i);
		_addFingerPrint2index(fingerPrint, index, fpLen, features, i);
	}

	free (features);

	return index;
}

void _free_index (void *source) {
	_index *index = (_index *)source;

	int numBits = index->numFeatures;

	for (int i=0; i<numBits; i++) {
		_free_moleculeList (index->moleculeList[i]);
	}

	free (index->moleculeList);
	free (index);
}

void _dataDistribution_small(void *source, char*fName) {
	FILE *fp = (FILE *) fopen (fName, "w");

	if (!fp)
		return;

	_index *index = (_index *)source;
	DataSmall *data = index->data;

	u_long dSize = data->size;
	u_long lenLimit = dSize/10;

	//length distribution
	fprintf (fp,"len, #mol\n");

	for (int len =data->minNumFeatures; len<=data->maxNumFeatures; len++) {

		if (data->_sizeIndexStart[len] >= 0) {
			int size = data->_sizeIndexEnd[len] - data->_sizeIndexStart[len] + 1;
			fprintf (fp,"%d, %d\n", len, size);
		}
	}

	//feature list distribution
	arrayListtype *hevayHitters = new_arrayListtype(1000);

	fprintf (fp,"\nfeatureId, #mol\n");

	int numBits = data->fingerPrintLen * NUMBITS_fingerPrintIntType;
	for (int f=0;f<numBits;f++) {
		_moleculeList *mlist = index->moleculeList[f];

		u_long size = mlist->arr->size;
		fprintf (fp,"%d, %lu\n",f,size);

		if (size >= lenLimit)
			arrayList_put(hevayHitters, f);
	}

	fclose (fp);

	//print information about heavy hitters
	int fpLen = data->fingerPrintLen;
	int *featureIndex = (int *)helper_malloc((index->numFeatures+1)*sizeof(int),__FILE__, __LINE__);

	arrayList_Iterator_type *iterator = arrayList_getIterator(hevayHitters);

	double totalHH = 0;
	int maxHH = 0;
	for (int molId=0; molId < dSize; molId++) {
		fingerPrintIntType *fp =  data_getFingerPrint (data, molId);

		int n = data_getIndexOf1Bits (featureIndex, fp, fpLen);

		int numHH=0;
		for (int i=n-1; i>=0; i--){
			int feature = featureIndex[i];

			int found = 0;
			int notenoughelements = arrayList_IteratorFind(iterator, feature, &found);
			if (notenoughelements)
				break;

			if (found) {
				numHH++;
			}
		}

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
/**********************************************/

/******************* Helpers *****************************/

void _set2OneInArray (int *candidates, int index) {
	candidates[index] = 1;
}

int _getIndexRangeOfCandidates(int k, int L_q, double similarityThreshold,
		int numBits,int *startIndex,  int *endIndex, _moleculeList* mList) {
	int discard = 1;
	//what is the max similarity we can get
	double numerator = k + 1;
	double denominator = L_q - (k + 1);
	// not perfect but good enough
	int minL_m = mList->minNum1Bits;
	if (minL_m > (k + 1))
		denominator += minL_m;
	else
		denominator += k + 1;

	double Tmax = numerator / denominator;
	if (Tmax < similarityThreshold) {
		//we can discard the feature
	} else {
		do {
			//we cannot discard the entire feature
			//find the boundary
			double b = (similarityThreshold + 1) * (k + 1);
			b /= similarityThreshold;
			b -= L_q;
			//int B_u = ceil (b); // originally we only used "ceil(b)" but floating point precision is causing trouble
			int B_u = ceil(b) + 1;
			if (B_u > numBits)
				B_u = numBits;
			int B_l = floor(similarityThreshold * L_q);

			if (B_u < B_l)
				break;

			//find molecules associated with this feature of size < B_u and >= B_l
			for (int i = B_l; i <= numBits; i++) {
				int index = mList->_sizeIndexStart[i];
				if (index >= 0) {
					*startIndex = index;
					break;
				}
			}
			if (*startIndex < 0)
				break;

			for (int i = B_u - 1; i > 0; i--) {
				if (mList->_sizeIndexStart[i] >= 0) {
					*endIndex = mList->_sizeIndexEnd[i];
					break;
				}
			}
			if (*endIndex < 0)
				break;

			discard = 0;
		} while (0);
	}
	return discard;
}


int _getLengthIndexRangeOfCandidates(int k, int L_q, double similarityThreshold,
		int numBits,int *startLength,  int *endLength,
		int *startIndex,  int *endIndex,_moleculeList* mList) {
	int discard = 1;
	//what is the max similarity we can get
	double numerator = k + 1;
	double denominator = L_q - (k + 1);
	// not perfect but good enough
	int minL_m = mList->minNum1Bits;
	if (minL_m > (k + 1))
		denominator += minL_m;
	else
		denominator += k + 1;

	double Tmax = numerator / denominator;
	if (Tmax < similarityThreshold) {
		//we can discard the feature
	} else {
		do {
			//we cannot discard the entire feature
			//find the boundary
			double b = (similarityThreshold + 1) * (k + 1);
			b /= similarityThreshold;
			b -= L_q;
			//int B_u = ceil (b); // originally we only used "ceil(b)" but floating point precision is causing trouble
			int B_u = ceil(b) + 1;
			if (B_u > numBits)
				B_u = numBits;
			int B_l = floor(similarityThreshold * L_q);

			if (B_u < B_l)
				break;

			//find molecules associated with this feature of size < B_u and >= B_l
			for (int i = B_l; i <= numBits; i++) {
				int index = mList->_sizeIndexStart[i];
				if (index >= 0) {
					*startLength = i;
					*startIndex = index;
					break;
				}
			}
			if (*startIndex < 0)
				break;

			for (int i = B_u - 1; i > 0; i--) {
				if (mList->_sizeIndexStart[i] >= 0) {
					*endIndex = mList->_sizeIndexEnd[i];
					*endLength = i;
					break;
				}
			}
			if (*endIndex < 0)
				break;

			discard = 0;
		} while (0);
	}
	return discard;
}


int _getLengthRangeOfCandidates(int k, int L_q, double similarityThreshold,
		int numBits,int *startLength,  int *endLength, _moleculeList* mList) {
	int discard = 1;
	//what is the max similarity we can get
	double numerator = k + 1;
	double denominator = L_q - (k + 1);
	// not perfect but good enough
	int minL_m = mList->minNum1Bits;
	if (minL_m > (k + 1))
		denominator += minL_m;
	else
		denominator += k + 1;

	double Tmax = numerator / denominator;
	if (Tmax < similarityThreshold) {
		//we can discard the feature
	} else {
		do {
			//we cannot discard the entire feature
			//find the boundary
			double b = (similarityThreshold + 1) * (k + 1);
			b /= similarityThreshold;
			b -= L_q;
			//int B_u = ceil (b); // originally we only used "ceil(b)" but floating point precision is causing trouble
			int B_u = ceil(b) + 1;
			if (B_u > numBits)
				B_u = numBits;
			int B_l = floor(similarityThreshold * L_q);

			if (B_u < B_l)
				break;

			//find molecules associated with this feature of size < B_u and >= B_l
			for (int i = B_l; i <= numBits; i++) {
				int index = mList->_sizeIndexStart[i];
				if (index >= 0) {
					discard = 0;
					*startLength = i;
					break;
				}
			}

			for (int i = B_u - 1; i > 0; i--) {
				if (mList->_sizeIndexStart[i] >= 0) {
					*endLength = i;
					break;
				}
			}
		} while (0);
	}
	return discard;
}

void _updatesCandidatesRange(int *candidates, double similarityThreshold,
		int L_q, int k, _moleculeList *mList, int numBits) {
	int start=-1;
	int end=-1;


	int discard = _getIndexRangeOfCandidates(k, L_q, similarityThreshold,
			numBits, &start, &end, mList);

	if (!discard) {
		//update the candidates
		u_long *arr = mList->arr->data;

		for (int i=start; i<=end; i++) {
			u_long fpIndex 		= arr[i];
			candidates[fpIndex] = 1;
		}
	}
}

/**
 * sort features in descending order size of mList
 */
void _sortFeatures (int *featureIndex, int num1Bits, _index *index) {
	int *size = (int *)helper_malloc((num1Bits)*sizeof(int),__FILE__, __LINE__);

	for (int i=0; i<num1Bits; i++) {
		_moleculeList *mList= index->moleculeList[featureIndex[i]];
		size[i] = mList->arr->size;
	}

	//sort
	for (int i=0; i<num1Bits; i++) {
		for (int k=i+1; k<num1Bits; k++) {
			if (size[i] < size[k]) {
				//swap len
				{
					int temp= size[i];
					size[i]	= size[k];
					size[k]	= temp;
				}

				//swap the mlist
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


arrayListtype* _runRangeSmallHelper_1(void *source, fingerPrintIntType *queryFP , double r,
		double *pruned, double *unPruned) {
	_index *index = (_index *)source;
	DataSmall *data = index->data;
	int dSize 		= data->size;
	int fpLen 		= data->fingerPrintLen;
	int numFeatures = index->numFeatures;

	int maxNumFeatures = data->maxNumFeatures;
	int minNumFeatures = data->minNumFeatures;

	int *featureIndex = (int *)helper_malloc((numFeatures+1)*sizeof(int),__FILE__, __LINE__);
	int num1Bits = data_getIndexOf1Bits (featureIndex, queryFP, fpLen);

	_sortFeatures (featureIndex, num1Bits, index);

	int *candidates = (int *) helper_calloc(dSize, sizeof(int),__FILE__, __LINE__);
	int L_q = num1Bits;

	arrayListtype *solution = new_arrayListtype(_ARRAY_BUCKET_SIZE_II);

	//process the candidates

	int k = floor (r*L_q); //number of features to be pruned
	int L_l = ceil(r*L_q);
	if (L_l < minNumFeatures)
		L_l = minNumFeatures;
	int L_u = floor (L_q/r);
	if (L_u > maxNumFeatures)
		L_u = maxNumFeatures;

	/*
	//creat the special fingerprint out of the discarded features
	//TODO
	fingerPrintIntType *specialQueryFP;
	{
		minHeapPlain_Type *heap = new_minHeapPlain_Type(k);
		for (int i=0; i<k; i++) {
			u_long feature = featureIndex[i];
			minHeapPlain_add(heap, feature, 0);
		}
		for (int i=0; i<k; i++) {
			u_long feature;
			minHeapPlain_pop(heap, &feature);
			//TODO arrayList_put(specialQuery_fp, feature);
		}
		minHeapPlain_free(heap);
	}
	 */
#ifdef G_GET_PRUNING_RATE
	int *allC  = (int *) helper_calloc(data->size, sizeof(int),	__FILE__, __LINE__);

	for (int i=0; i<num1Bits; i++) {
		int feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList *mList= index->moleculeList[feature];
		u_long *molid = mList->arr->data;
		u_long msize  = mList->arr->size;
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

	for (int i=k; i<num1Bits; i++) {
		int feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList *mList= index->moleculeList[feature];

		int startIndex=-1;
		for (int l=L_l; l<=L_u; l++) {
			if (mList->_sizeIndexStart[l]>0) {
				startIndex = mList->_sizeIndexStart[l];
				break;
			}
		}
		if (startIndex >= 0) {
			int endIndex=-1;
			for (int l=L_u; l>=L_l; l--) {
				if (mList->_sizeIndexStart[l]>0) {
					endIndex = mList->_sizeIndexEnd[l];
					break;
				}
			}

			u_long *arr = mList->arr->data;
			for (int i = startIndex; i <= endIndex; i++) {
				u_long fpIndex = arr[i];
				candidates[fpIndex]++;
			}
		}
	}

	//process the candidates which were not processed earlier
	for (int l=L_l; l<=L_u; l++) {
		int startIndex = data->_sizeIndexStart[l];
		if (startIndex >= 0) {
			int endIndex = data->_sizeIndexEnd[l];
			double minimumIntersection = (r/(r+1)*(l+L_q));
			minimumIntersection = ((double)floor(G_PRECISION * minimumIntersection))/G_PRECISION;
			double m = minimumIntersection-k;
			for (int fpIndex = startIndex; fpIndex <= endIndex; fpIndex++) {
				int c = candidates[fpIndex];
				if (c) {

#ifdef G_GET_PRUNING_RATE
					*pruned = *pruned + 1;
#endif

					if (c >= minimumIntersection)
						arrayList_put(solution, fpIndex);
					else if (c >= m) {
						fingerPrintIntType *dataFP = data_getFingerPrint(data, fpIndex);

						//int intersection = data_getIntersection_fast(specialQuery_fp, dataFP) + c;
						int intersection = data_getIntersection(queryFP, dataFP, fpLen);

						if (intersection >= (minimumIntersection)) {
							arrayList_put(solution, fpIndex);
						}
					}
				}
			}
		}
	}

	//check the solution
	if (CHECK_SOLUTION_G) {
		printf ("comparing the 2 solutions : ");
		arrayListtype *solutionLinear = new_arrayListtype(_ARRAY_BUCKET_SIZE_II);
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
		} else {
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

	free (candidates);
	free (featureIndex);

	return solution;
}

void _runTopKSmallHelperHelper (minHeapInttype *solutionHeap,int *featureIndex, int a,
		_index *index, int featuresPruned, int L_q, int numBitsInFP, int *processed,
		DataSmall *data, fingerPrintIntType *queryFP,int fpLen) {

	double r = minHeapInt_peekKey(solutionHeap);

	int feature = featureIndex[a];

	//get the arraylist of molecules
	_moleculeList *mList= index->moleculeList[feature];

	u_long *arr = mList->arr->data;
	int startLength, endLength;

	int discard  = _getLengthRangeOfCandidates (featuresPruned, L_q, r, numBitsInFP,
			&startLength, &endLength, mList);

	if (!discard) {
		int leftLength, rightLength;
		int sizeL=0, sizeR=0;
		int li=-1, ri=-1;

		//find the valid length closest to L_q
		{
			if (L_q < endLength)
				leftLength = L_q;
			else
				leftLength = endLength;
			rightLength = leftLength+1;

			for (;leftLength >= startLength; leftLength--) {
				int val = mList->_sizeIndexStart[leftLength];

				if (val >= 0) {
					int leftIndex = mList->_sizeIndexEnd[leftLength];
					sizeL = leftIndex-val+1;
					li 		= leftIndex;
					break;
				}
			}
			for (;rightLength <= endLength; rightLength++) {
				int val = mList->_sizeIndexStart[rightLength];

				if (val >= 0) {
					int rightIndex = val;
					sizeR = mList->_sizeIndexEnd[rightLength]-val+1;
					ri 		= rightIndex;
					break;
				}
			}
		}

		int more2Explore;
		do{
			more2Explore=0;

			if (startLength<= leftLength) {
				//explore left
				r = minHeapInt_peekKey(solutionHeap);
				double adjuster = 1;
				double l_thresh = r * L_q - adjuster; //subtracted 0.5 to avoid floating point problems
				if (l_thresh > leftLength) {
					//no need to explore further
					leftLength=startLength-1;
					break;
				}
				more2Explore=1;
				double minimumIntersection = (r/(r+1)*(leftLength+L_q));
				minimumIntersection = ((double)floor(G_PRECISION * minimumIntersection))/G_PRECISION;

				for (;sizeL;sizeL--) {
					int molIndex = arr[li];
					if (!processed[molIndex]) {
						processed[molIndex] = 1;
						fingerPrintIntType *dataFP = data_getFingerPrint(data, molIndex);
						double intersection = data_getIntersection(queryFP, dataFP, fpLen);

						if (intersection >= minimumIntersection) {
							double sim = 0;
							int u = leftLength+L_q-intersection;
							sim = intersection/u;
							minHeapInt_addIfKeyHigher(solutionHeap, sim, molIndex);

							double r1 = minHeapInt_peekKey(solutionHeap);
							if (r1 > r) {
								r=r1;
								double l_thresh = r * L_q - adjuster;
								if (l_thresh > leftLength) {
									//no need to explore further
									leftLength=startLength-1;
									break;
								}
								minimumIntersection = (r/(r+1)*(leftLength+L_q));
							}
						}
					}

					//update index
					li--;
				}

				for (leftLength--;leftLength >= startLength; leftLength--) {
					int val = mList->_sizeIndexStart[leftLength];
					if (val >= 0) {
						int leftIndex = mList->_sizeIndexEnd[leftLength];
						sizeL = leftIndex-val+1;
						break;
					}
				}
			}

			if (rightLength<= endLength) {
				r = minHeapInt_peekKey(solutionHeap);
				double adjuster = 0.5;
				double l_thresh = (double)L_q/r+adjuster; //added 0.5 to avoid floating point
				if (rightLength > l_thresh) {
					//no need explre further
					rightLength = endLength+1;
					break;
				}
				more2Explore=1;
				double minimumIntersection = (r/(r+1)*(rightLength+L_q));
				minimumIntersection = ((double)floor(G_PRECISION * minimumIntersection))/G_PRECISION;

				//explore right
				for (;sizeR;sizeR--) {
					int molIndex = arr[ri];
					if (!processed[molIndex]) {
						processed[molIndex] = 1;

						fingerPrintIntType *dataFP = data_getFingerPrint(data, molIndex);
						double intersection = data_getIntersection(queryFP, dataFP, fpLen);
						if (intersection >= minimumIntersection) {
							double sim = 0;
							int u = rightLength+L_q-intersection;
							sim = intersection/u;
							minHeapInt_addIfKeyHigher(solutionHeap, sim, molIndex);

							double r1 = minHeapInt_peekKey(solutionHeap);
							if (r1 > r) {
								r=r1;
								l_thresh = (double)L_q/r+adjuster;
								if (rightLength > l_thresh) {
									//no need explre further
									rightLength = endLength+1;
									break;
								}
								minimumIntersection = (r/(r+1)*(rightLength+L_q));
							}
						}
					}

					//update index
					ri++;
				}

				for (rightLength++;rightLength <= endLength; rightLength++) {
					int val = mList->_sizeIndexStart[rightLength];

					if (val >= 0) {
						sizeR = mList->_sizeIndexEnd[rightLength]-val+1;
						break;
					}
				}
			}
		}while(more2Explore);
	}
}

void _runTopKHelper_explore (int len, _index *index, fingerPrintIntType *queryFP,
		minHeapInttype *solutionHeap, int *candidates, int *featureIndex,
		int p, int L_q, double minimumIntersection, int fpLen, double *pruned) {
	//from all feature list we will only explore the molecules of length len
	for (int i=p; i<L_q; i++) {
		u_long feature = featureIndex[i];

		_moleculeList *mList = index->moleculeList[feature];

		if (mList->_sizeIndexStart[len] < 0)
			continue;

		int startIndex 	= mList->_sizeIndexStart[len];
		int endIndex 	= mList->_sizeIndexEnd  [len];
		u_long *arr = mList->arr->data;
		for (int i = startIndex; i <= endIndex; i++) {
			u_long fpIndex = arr[i];
			candidates[fpIndex]++;
		}
	}

	//process the candidates
	DataSmall *data = index->data;
	int startIndex = data->_sizeIndexStart[len];
	if (startIndex >= 0) {
		int endIndex = data->_sizeIndexEnd[len];

		double m = minimumIntersection-p;
		for (int fpIndex = startIndex; fpIndex <= endIndex; fpIndex++) {
			int c = candidates[fpIndex];
			if (c) {

#ifdef G_GET_PRUNING_RATE
				*pruned = *pruned + 1;
#endif

				if (c >= m) {
					fingerPrintIntType *dataFP = data_getFingerPrint (data, fpIndex);
					int intersection= data_getIntersection(queryFP, dataFP, fpLen);
					int u 			= len+L_q-intersection;
					double sim 		= (double)intersection/u;

					minHeapInt_addIfKeyHigher(solutionHeap, sim, fpIndex);
				}
			}
		}
	}
}

minHeapInttype* _runTopKSmallHelper_1(void *source, fingerPrintIntType *queryFP , int k, double *pruned, double *unPruned) {
	_index *index 	= (_index *)source;
	DataSmall *data = index->data;
	int dSize 		= data->size;
	int fpLen 		= data->fingerPrintLen;
	int maxNumFeatures = data->maxNumFeatures;
	int minNumFeatures = data->minNumFeatures;

	int *featureIndex = (int *)helper_malloc((index->numFeatures+1)*sizeof(int),__FILE__, __LINE__);
	int num1Bits = data_getIndexOf1Bits (featureIndex, queryFP, fpLen);

	_sortFeatures (featureIndex, num1Bits, index);

	int *candidates  = (int *) helper_calloc(dSize, sizeof(int),__FILE__, __LINE__);

	int L_q = num1Bits;

 	double defaultKey = 1.0/(L_q+data->maxNumFeatures-1);
	minHeapInttype *solutionHeap = newMinHeapInttype_withDefault(k, defaultKey);

	int moveLeft=1, moveRight=1;
	int leftLen	= L_q;
	int rightLen= leftLen+1;

#ifdef G_GET_PRUNING_RATE
	int *allC  = (int *) helper_calloc(data->size, sizeof(int),	__FILE__, __LINE__);
	for (int i=0; i<num1Bits; i++) {
		int feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList *mList= index->moleculeList[feature];
		u_long *molid = mList->arr->data;
		u_long msize  = mList->arr->size;
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

	while (moveLeft || moveRight){
		if (moveLeft) {
			double r = minHeapInt_peekKey_or_default(solutionHeap);
			int L_l = ceil (r*L_q);
			if (L_l < minNumFeatures)
				L_l = minNumFeatures;

			if (leftLen < L_l) {
				moveLeft = 0;
			} else {
				int p = floor (r*L_q);
				double m = (r/(r+1)*(leftLen+L_q));
				m = ((double)floor(G_PRECISION * m))/G_PRECISION;
				_runTopKHelper_explore (leftLen, index, queryFP, solutionHeap,
						candidates, featureIndex, p, L_q, m, fpLen, pruned);
				leftLen--;
			}
		}

		if (moveRight) {
			double r = minHeapInt_peekKey_or_default(solutionHeap);
			int L_u  = floor (L_q/r);
			if (L_u > maxNumFeatures)
				L_u = maxNumFeatures;

			if (rightLen > L_u) {
				moveRight = 0;
			} else {
				int p = floor (r*L_q);
				double m = (r/(r+1)*(rightLen+L_q));
				m = ((double)floor(G_PRECISION * m))/G_PRECISION;
				_runTopKHelper_explore (rightLen, index, queryFP, solutionHeap,
						candidates, featureIndex, p, L_q, m, fpLen, pruned);
				rightLen++;
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
				helper_err(__FILE__ , __LINE__, "Solution not matching");
			}
		}

		minHeapInt_free(solutionHeapLinear);
	}

	free (candidates);
	free (featureIndex);
	return (solutionHeap);
}

void _addFPIndex2moleculeList (_moleculeList *mList, int fpIndex, int num1Bits) {
	arrayList_flexible_put(&(mList->arr), fpIndex);

	int index = mList->arr->size-1;

	if (mList->minNum1Bits > num1Bits)
		mList->minNum1Bits = num1Bits;
	if (mList->maxNum1Bits < num1Bits)
		mList->maxNum1Bits = num1Bits;

	//adjust the start and end index of this size
	if (mList->_sizeIndexStart[num1Bits] == -1) {
		mList->_sizeIndexStart[num1Bits]=index;
	}
	mList->_sizeIndexEnd[num1Bits]=index;
}

workerFunctions_type_small invertedIndex_getWorkerFunction () {
	workerFunctions_type_small worker;

	worker.dataDistribution = _dataDistribution_small;
	worker.init_index 	= _init_index;
	worker.free_index	= _free_index;
	worker.runRange		= _runRangeSmallHelper_1;
	worker.runTopK		= _runTopKSmallHelper_1;
	return worker;
}
