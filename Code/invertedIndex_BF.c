/*
 * invertedIndex_BF.c
 *
 *  Created on: 22-Dec-2018
 *      Author: jithin
 */

#include "invertedIndex_BF.h"

#define _ARRAY_BUCKET_SIZE_II_BF 1000
//#define _STAT_PRUNING_BF

/******************* Typedefs ****************************/
typedef struct {
	arrayListtype_flexible *arr;
	int minNumFeatures;				// the size of the smallest mol
	int maxNumFeatures;				// the size of the largest mol
	u_long *_sizeIndexStart_1;		// stores the starting index of molecules of given size
	u_long *_sizeIndexEnd_1;		// stores the ending index(inclusive) of molecules of given size
	char *_sizeIndexValid;			//set to one if elements of this size exist
}_moleculeList_BF;

typedef struct {
	int numFeatures; 				// the number of features in the index
	DataBinary *data;
	_moleculeList_BF **moleculeList;	//the array of index entries
}_index_BF;

/******************* Helpers *****************************/
void _addFPIndex2moleculeList_BF 	(_moleculeList_BF *mList, int fpIndex, int numFeatures);

 /******************* Initializers ************************/
_moleculeList_BF* _new_moleculeList_BF (int maxNumFeatures, int dSize) {
	_moleculeList_BF *mList = (_moleculeList_BF *) helper_calloc(1, sizeof(_moleculeList_BF),
			__FILE__, __LINE__);

	mList->_sizeIndexStart_1= (u_long *) helper_malloc(sizeof(u_long) * (maxNumFeatures + 1),
			__FILE__, __LINE__);
	mList->_sizeIndexEnd_1	= (u_long *) helper_malloc(sizeof(u_long) * (maxNumFeatures + 1),
			__FILE__, __LINE__);
	mList->_sizeIndexValid	= (char *) helper_calloc((maxNumFeatures + 1), sizeof(char),
			__FILE__, __LINE__);

	mList->arr = new_arrayListtype_flexible(_ARRAY_BUCKET_SIZE_II_BF, dSize);
	mList->minNumFeatures = INT_MAX;
	mList->maxNumFeatures = 0;

	return mList;
}

void _free_moleculeList_BF (_moleculeList_BF *mList) {

	free_arrayListtype_flexible(mList->arr);
	free (mList->_sizeIndexStart_1);
	free (mList->_sizeIndexEnd_1);
	free (mList->_sizeIndexValid);
	free (mList);
}

/************************ index ************************/
void _addFingerPrint2index_BF (arrayListtype *fp, _index_BF *index, int fpIndex) {
	int numFeatures = fp->size;

	arrayList_Iterator_type *iterator = arrayList_getIterator(fp);

	int n = numFeatures;

	while (n--) {
		int feature = arrayList_IteratorGetNext(iterator);

		if (feature >= index->numFeatures)
			helper_err(__FILE__, __LINE__, "the feature value is more than the possible");

		_moleculeList_BF *mList = index->moleculeList[feature];
		_addFPIndex2moleculeList_BF (mList, fpIndex, numFeatures);
	}
	arrayList_IteratorFree(iterator);
}

void* _init_index_BF (DataBinary *data) {
	_index_BF *index = (_index_BF *) helper_calloc (1,sizeof(_index_BF),
			__FILE__, __LINE__);
	//allocate memory
	{
		int numFeatures = data->numFeatures;
		index->numFeatures 	= numFeatures;
		index->data			= data;

		_moleculeList_BF **mList = (_moleculeList_BF **) helper_malloc ((size_t)numFeatures * sizeof(_moleculeList_BF **),
				__FILE__, __LINE__);

		for (int i=0; i<numFeatures; i++) {
			mList[i] = _new_moleculeList_BF (data->maxNumFeatures, data->size);
		}

		index->moleculeList = mList;
	}

	// load data
	int size 		= data->size;

	for (int i=0; i<size; i++) {
		arrayListtype *fingerPrint = dataBinary_getFingerPrint(data, i);
		_addFingerPrint2index_BF(fingerPrint, index, i);
	}
	return index;
}

void _free_index_BF (void *source) {
	_index_BF *index = (_index_BF *) source;
	int numFeatures = index->numFeatures;
	for (int i=0; i<numFeatures; i++) {
		_free_moleculeList_BF(index->moleculeList[i]);
	}

	free (index->moleculeList);
	free (index);
}

void _dataDistribution_BF(void *source, char*fName) {
	FILE *fp = (FILE *) fopen (fName, "w");

	if (!fp)
		return;

	_index_BF *index = (_index_BF *)source;
	DataBinary *data = index->data;

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

	//feature list distribution
	u_long dSize = data->size;
	u_long lenLimit = dSize/10;

	arrayListtype *hevayHitters = new_arrayListtype(1000);
	int numFeatures = index->numFeatures;
	for (int f=0;f<numFeatures;f++) {
		_moleculeList_BF *mlist = index->moleculeList[f];

		u_long size = mlist->arr->size;
		fprintf (fp,"%d, %lu\n",f,size);

		if (size >= lenLimit)
			arrayList_put(hevayHitters, f);
	}

	fclose (fp);

	//print information about heavy hitters

	arrayList_Iterator_type *iterator = arrayList_getIterator(hevayHitters);

	double totalHH = 0;
	int maxHH = 0;
	for (int molId=0; molId < dSize; molId++) {
		arrayListtype *fp =  dataBinary_getFingerPrint (data, molId);

		int n = fp->size;

		int numHH=0;
		arrayList_Iterator_type *fiterator = arrayList_getIterator(fp);
		for (int i=0; i<n; i++){
			int feature = arrayList_IteratorGetNext(fiterator);

			int found = 0;
			int notenoughelements = arrayList_IteratorFind(iterator, feature, &found);
			if (notenoughelements)
				break;

			if (found) {
				numHH++;
			}
		}
		arrayList_IteratorFree(fiterator);
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

/**
 * sort features in descending order size of mList
 */
void _sortFeatures_BF (u_long *featureIndex, arrayListtype *queryFP, _index_BF *index) {
	int L_q = queryFP->size;
	int size[L_q];

	arrayList_Iterator_type *iterator = arrayList_getIterator(queryFP);
	for (int i=0; i<L_q; i++) {
		u_long feature = arrayList_IteratorGetNext(iterator);
		featureIndex[i] = feature;
		_moleculeList_BF *mList= index->moleculeList[feature];
		size[i] = mList->arr->size;
	}
	arrayList_IteratorFree(iterator);

	//sort
	for (int i=0; i<L_q; i++) {
		for (int k=i+1; k<L_q; k++) {
			if (size[i] < size[k]) {
				//swap len
				{
					int temp= size[i];
					size[i]	= size[k];
					size[k]	= temp;
				}

				//swap the mlist
				{
					u_long temp = featureIndex[i];
					featureIndex[i] = featureIndex[k];
					featureIndex[k] = temp;
				}
			}
		}
	}
}

int _getIndexRangeOfCandidates_BF(int k, int L_q, double similarityThreshold,
		u_long *startIndex,  u_long *endIndex, _moleculeList_BF *mList, int maxNumFeatures) {
	int discard = 1;
	//what is the max similarity we can get
	double numerator = k + 1;
	double denominator = L_q - (k + 1);
	// not perfect but good enough
	int minL_m = mList->minNumFeatures;
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
			if (B_u > maxNumFeatures)
				B_u = maxNumFeatures;
			int B_l = floor(similarityThreshold * L_q);

			if (B_u < B_l)
				break;

			//find molecules associated with this feature of size < B_u and >= B_l
			int i;
			for (i = B_l; i < B_u; i++) {
				if (mList->_sizeIndexValid[i]) {
					*startIndex = mList->_sizeIndexStart_1[i];
					break;
				}
			}
			if (i >= B_u)
				break;

			for (i = B_u - 1; i >= B_l; i--) {
				if (mList->_sizeIndexValid[i]) {
					*endIndex = mList->_sizeIndexEnd_1[i];
					break;
				}
			}

			discard = 0;
		} while (0);
	}
	return discard;
}

void _updatesCandidatesRange_BF(int *candidates, double similarityThreshold,
		int L_q, int k, _moleculeList_BF *mList, int maxNumFeatures) {
	u_long start=-1;
	u_long end=-1;


	int discard = _getIndexRangeOfCandidates_BF(k, L_q, similarityThreshold,
			&start, &end, mList, maxNumFeatures);

	if (!discard) {
		//update the candidates
		arrayListtype_flexible *arr = mList->arr;
		arrayList_flexible_iterateAndWork_3(arr, start, end, candidates);
	}
}

int _get_L_l_BF (int L_q, double r, int minNumFeatures) {
	int L_l = ceil (r*L_q);
	if (L_l < minNumFeatures)
		return minNumFeatures;
	else
		return L_l;
}

int _get_L_u_BF (int L_q, double r, int maxNumFeatures) {
	int L_u = floor((double)L_q/r);
	if (L_u > maxNumFeatures)
		return maxNumFeatures;
	else
		return L_u;
}

minHeapInttype* _runTopKHelper_BF (void *source, arrayListtype *queryFP , int k,
		double *pruned, double *unPruned) {
	_index_BF *index = (_index_BF *) source;
	DataBinary *data = index->data;

	int L_q = queryFP->size;
	u_long featureIndex[L_q];

	_sortFeatures_BF (featureIndex, queryFP, index);

	int *processed  = (int *) helper_calloc(data->size, sizeof(int),	__FILE__, __LINE__);

	double defaultKey = 1.0/(L_q+index->data->maxNumFeatures-1);
	minHeapInttype *solutionHeap = newMinHeapInttype_withDefault(k, defaultKey);

	//populate the initial K elements
	{
		for (int maxIntersection=L_q; maxIntersection>0; maxIntersection--) {
			double r = minHeapInt_peekKey_or_default(solutionHeap);
			int feature = featureIndex[maxIntersection-1];

			_moleculeList_BF *mList = index->moleculeList[feature];

			int diff = mList->minNumFeatures - maxIntersection;
			if (diff < 0)
				diff = 0;
			double denominator = (L_q + diff);
			double maxSim = ((double) maxIntersection)/denominator;

			if (maxSim < r) {
				//we can discard all other list
				break;
			} else {
				//we need to process this

				//process for L_q down to L_l and
				//process for L_q+1 upto to L_u

				//process left side
				int minNumFeatures = mList->minNumFeatures;
				int L_l = _get_L_l_BF(L_q, r, minNumFeatures);
				for (int l = L_q; l >= L_l; l--) {
					if (mList->_sizeIndexValid[l]) {
						int startIndex = mList->_sizeIndexStart_1[l];
						double minimumIntersection = (r/(r+1)*(l+L_q));
						minimumIntersection = ((double)floor(G_PRECISION * minimumIntersection))/G_PRECISION;
						int endIndex= mList->_sizeIndexEnd_1[l];
						u_long *fpId 	= mList->arr->data;

						for (int f = endIndex; f>= startIndex; f--) {
							int fpIndex = fpId[f];
							if (!processed[fpIndex]) {
								processed[fpIndex] = 1;

								arrayListtype *dataFP = dataBinary_getFingerPrint (data,   fpIndex);
								int intersection      = dataBinary_getIntersection(queryFP, dataFP);

								if (intersection >= minimumIntersection) {
									int u = l+L_q-intersection;
									double sim = (double)intersection/u;

									minHeapInt_addIfKeyHigher(solutionHeap, sim, fpIndex);
								}
							}
						}

						r = minHeapInt_peekKey_or_default(solutionHeap);
						L_l = _get_L_l_BF(L_q, r, minNumFeatures);
					}
				}

				//process right
				int maxNumFeatures = mList->maxNumFeatures;
				int L_u = _get_L_u_BF(L_q, r, maxNumFeatures);
				for (int l = L_q+1; l <= L_u; l++) {
					if (mList->_sizeIndexValid[l]) {
						int startIndex = mList->_sizeIndexStart_1[l];
						double minimumIntersection = (r/(r+1)*(l+L_q));
						minimumIntersection = ((double)floor(G_PRECISION * minimumIntersection))/G_PRECISION;
						int endIndex= mList->_sizeIndexEnd_1[l];
						u_long *fpID 	= mList->arr->data;

						for (int f = startIndex; f<= endIndex; f++) {
							int fpIndex = fpID[f];
							if (!processed[fpIndex]) {
								processed[fpIndex] = 1;

								arrayListtype *dataFP = dataBinary_getFingerPrint (data,   fpIndex);
								int intersection      = dataBinary_getIntersection(queryFP, dataFP);

								if (intersection >= minimumIntersection) {
									int u = l+L_q-intersection;
									double sim = (double)intersection/u;

									minHeapInt_addIfKeyHigher(solutionHeap, sim, fpIndex);
								}
							}
						}

						r = minHeapInt_peekKey_or_default(solutionHeap);
						L_u = _get_L_u_BF(L_q, r, maxNumFeatures);
					}
				}
			}
		}
	}
#ifdef _STAT_PRUNING_BF
	int numProcessed=0;
	for (int i=0; i<data->size; i++)
		numProcessed += processed[i];

	printf ("%d,\n",numProcessed);
#endif

	free (processed);

	//check the solution
	if (CHECK_SOLUTION_G) {
		static int c;
		printf ("%d comparing the 2 solutions\n",c++);
		if (c==40) {
			c=c+1-1;
		}
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
				helper_err(__FILE__ , __LINE__, "Solution not matching");
			}
		}

		minHeapInt_free(solutionHeapLinear);
	}

	return solutionHeap;
}

void _runTopKHelper_BF_explore (int len, _index_BF *index, arrayListtype *queryFP,
		minHeapInttype *solutionHeap, int *candidates, u_long *featureIndex,
		int p, int L_q, double minimumIntersection, double *pruned) {
	//from all feature list we will only explore the molecules of length len

	//creating a special query mol
	arrayListtype *specialQuery_fp;
	if (p){
		specialQuery_fp = new_arrayListtype(p);
		minHeapPlain_Type *heap = new_minHeapPlain_Type(p);
		for (int i=0; i<p; i++) {
			u_long feature = featureIndex[i];
			minHeapPlain_add(heap, feature, 0);
		}
		for (int i=0; i<p; i++) {
			u_long feature;
			minHeapPlain_pop(heap, &feature);
			arrayList_put(specialQuery_fp, feature);
		}
		minHeapPlain_free(heap);
	}

	for (int i=p; i<L_q; i++) {
		u_long feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList_BF *mList	= index->moleculeList[feature];

		if (!mList->_sizeIndexValid[len])
			continue;

		int startIndex 	= mList->_sizeIndexStart_1[len];
		int endIndex 	= mList->_sizeIndexEnd_1  [len];
		u_long *arr = mList->arr->data;
		for (int i = startIndex; i <= endIndex; i++) {
			u_long fpIndex = arr[i];
			candidates[fpIndex]++;
		}
	}

	//process the candidates
	DataBinary *data = index->data;
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
					arrayListtype *dataFP = dataBinary_getFingerPrint (data, fpIndex);
					int intersection=c;
					if (p)
						intersection += dataBinary_getIntersection_fast(specialQuery_fp, dataFP);
					int u 			= len+L_q-intersection;
					double sim 		= (double)intersection/u;

					minHeapInt_addIfKeyHigher(solutionHeap, sim, fpIndex);
				}
			}
		}
	}
}

minHeapInttype* _runTopKHelper_BF_1 (void *source, arrayListtype *queryFP , int k,
		double *pruned, double *unPruned) {
	G_LOCATION;
	_index_BF *index 	= (_index_BF *) source;
	DataBinary *data 	= index->data;
	int maxNumFeatures 	= data->maxNumFeatures;

	int L_q = queryFP->size;
	u_long featureIndex[L_q];

	_sortFeatures_BF (featureIndex, queryFP, index);

	int *candidates  = (int *) helper_calloc(data->size, sizeof(int),	__FILE__, __LINE__);

	double defaultKey = 1.0/(L_q+index->data->maxNumFeatures-1);
	minHeapInttype *solutionHeap = newMinHeapInttype_withDefault(k, defaultKey);

#ifdef G_GET_PRUNING_RATE
	int *allC  = (int *) helper_calloc(data->size, sizeof(int),	__FILE__, __LINE__);
	for (int i=0; i<L_q; i++) {
		u_long feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList_BF *mList= index->moleculeList[feature];

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


	int moveLeft=1, moveRight=1;
	int leftLen	= L_q;
	int rightLen= leftLen+1;

	while (moveLeft || moveRight){
		if (moveLeft) {
			double r = minHeapInt_peekKey_or_default(solutionHeap);
			int L_l = ceil (r*L_q);

			if (leftLen < L_l) {
				moveLeft = 0;
			} else {
				int p = floor (r*L_q);
				double m = (r/(r+1)*(leftLen+L_q));
				m = ((double)floor(G_PRECISION * m))/G_PRECISION;
				_runTopKHelper_BF_explore (leftLen, index, queryFP, solutionHeap,
						candidates, featureIndex, p, L_q, m, pruned);
				leftLen--;
			}
		}

		if (moveRight) {
			double r = minHeapInt_peekKey_or_default(solutionHeap);
			int L_u  = _get_L_u_BF(L_q, r, maxNumFeatures);

			if (rightLen > L_u) {
				moveRight = 0;
			} else {
				int p = floor (r*L_q);
				double m = (r/(r+1)*(rightLen+L_q));
				m = ((double)floor(G_PRECISION * m))/G_PRECISION;
				_runTopKHelper_BF_explore (rightLen, index, queryFP, solutionHeap,
						candidates, featureIndex, p, L_q, m, pruned);
				rightLen++;
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
				helper_err(__FILE__ , __LINE__, "Solution not matching");
			}
		}

		minHeapInt_free(solutionHeapLinear);
	}

	return solutionHeap;
}

arrayListtype* _runRangeHelper_BF(void * source, arrayListtype *queryFP , double r,
		double *pruned, double *unPruned) {
	_index_BF *index = (_index_BF *)source;

	DataBinary *data = index->data;
	int dSize = data->size;
	int maxNumFeatures = index->data->maxNumFeatures;

	int L_q = queryFP->size;
	u_long featureIndex[L_q];

	_sortFeatures_BF (featureIndex, queryFP, index);

	int *candidates = (int *) helper_calloc(dSize, sizeof(int), __FILE__, __LINE__);
#ifdef _STAT_PRUNING_BF
	int *uprunedCandidates = (int *) helper_calloc(dSize, sizeof(int), __FILE__, __LINE__);
#endif

	int k = 0; //number of features pruned


	for (int i=0; i<L_q; i++) {
		u_long feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList_BF *mList= index->moleculeList[feature];

		//add all molecules to the candidate list
		_updatesCandidatesRange_BF(candidates, r, L_q, k, mList, maxNumFeatures);

#ifdef _STAT_PRUNING_BF
		arrayList_flexible_iterateAndWork_3(arr, 0, arr->size-1, uprunedCandidates);
#endif
		k++;
	}

	arrayListtype *solution = new_arrayListtype(_ARRAY_BUCKET_SIZE_II_BF);

	//process the candidates
#ifdef _STAT_PRUNING_BF
	int numCandidates	=0;
	int numCandidatesU	=0;
#endif

	int L_l = ceil (r*L_q);
	int L_u = floor((double)L_q/r);
	if (L_u > maxNumFeatures)
		L_u = maxNumFeatures;

	for (int l=L_l; l<=L_u; l++) {
		int startIndex = data->_sizeIndexStart[l];
		if (startIndex >= 0) {
			int endIndex = data->_sizeIndexEnd[l];
			double minimumIntersection = (r/(r+1)*(l+L_q));
			minimumIntersection = ((double)floor(G_PRECISION * minimumIntersection))/G_PRECISION;
			for (int fpIndex = startIndex; fpIndex <= endIndex; fpIndex++) {
				if (candidates[fpIndex]) {
#ifdef _STAT_PRUNING_BF
					numCandidates++;
#endif
					arrayListtype *dataFP = dataBinary_getFingerPrint (data,   fpIndex);
					int intersection      = dataBinary_getIntersection(queryFP, dataFP);

					if (intersection >= minimumIntersection) {
						arrayList_put(solution, fpIndex);
					}
				}
			}
		}
	}

#ifdef _STAT_PRUNING_BF
	for (int i=0; i<dSize; i++)
		if (uprunedCandidates[i])
			numCandidatesU++;
	printf ("%d,%d\n",numCandidates, numCandidatesU);
#endif
	//check the solution
	if (CHECK_SOLUTION_G) {
		printf ("comparing the 2 solutions : ");
		arrayListtype *solutionLinear = new_arrayListtype(_ARRAY_BUCKET_SIZE_II_BF);
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

	return (solution);
}

arrayListtype* _runRangeHelper_BF_1(void * source, arrayListtype *queryFP , double r,
		double *pruned, double *unPruned) {
	_index_BF *index = (_index_BF *)source;

	DataBinary *data = index->data;
	int dSize = data->size;
	int maxNumFeatures = index->data->maxNumFeatures;

	int L_q = queryFP->size;
	u_long featureIndex[L_q];

	_sortFeatures_BF (featureIndex, queryFP, index);

	arrayListtype *solution = new_arrayListtype(_ARRAY_BUCKET_SIZE_II_BF);

#ifdef G_GET_PRUNING_RATE
	int *allC  = (int *) helper_calloc(data->size, sizeof(int),	__FILE__, __LINE__);
	for (int i=0; i<L_q; i++) {
		u_long feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList_BF *mList= index->moleculeList[feature];

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

	int *candidates = (int *) helper_calloc(dSize, sizeof(int), __FILE__, __LINE__);

	int k = floor (r*L_q); //number of features to be pruned
	int L_l = ceil (r*L_q);
	int L_u = floor((double)L_q/r);
	if (L_u > maxNumFeatures)
		L_u = maxNumFeatures;

	double minIntersection[L_u];
	for (int l=L_l; l<=L_u; l++) {
		double m = (r/(r+1)*(l+L_q));
		m = ((double)floor(G_PRECISION * m))/G_PRECISION;

		minIntersection[l] = m;
	}

	int m = ceil (minIntersection[L_l]) + k;
	if (m > L_q)
		m = L_q;

	//creating a special query mol
	arrayListtype *specialQuery_fp=NULL;
	if (k) {
		specialQuery_fp = new_arrayListtype(k);
		minHeapPlain_Type *heap = new_minHeapPlain_Type(k);
		for (int i=0; i<k; i++) {
			u_long feature = featureIndex[i];
			minHeapPlain_add(heap, feature, 0);
		}
		for (int i=0; i<k; i++) {
			u_long feature;
			minHeapPlain_pop(heap, &feature);
			arrayList_put(specialQuery_fp, feature);
		}
		minHeapPlain_free(heap);
	}

	for (int i=k; i<m; i++) {
		u_long feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList_BF *mList	= index->moleculeList[feature];

		int startIndex=-1;
		for (int l=L_l; l<=L_u; l++) {
			if (mList->_sizeIndexValid[l]) {
				startIndex = mList->_sizeIndexStart_1[l];
				break;
			}
		}
		if (startIndex >= 0) {
			int endIndex=-1;
			for (int l=L_u; l>=L_l; l--) {
				if (mList->_sizeIndexValid[l]) {
					endIndex = mList->_sizeIndexEnd_1[l];
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

	//now we may compare the number of intersections
	for (int i=m; i<L_q; i++) {
		u_long feature = featureIndex[i];

		//get the arraylist of molecules
		_moleculeList_BF *mList	= index->moleculeList[feature];

		for (int l=L_l; l<=L_u; l++) {
			u_long *arr = mList->arr->data;
			if (mList->_sizeIndexValid[l]) {
				int startIndex 	=  mList->_sizeIndexStart_1[l];
				int endIndex 	=  mList->_sizeIndexEnd_1  [l];

				//double minimumIntersection = minIntersection[l]-1;

				for (int i = startIndex; i <= endIndex; i++) {
					u_long fpIndex = arr[i];
					candidates[fpIndex]++;
				}
			}
		}
	}


	//process the candidates which were not processed earlier
	for (int l=L_l; l<=L_u; l++) {
		int startIndex = data->_sizeIndexStart[l];
		if (startIndex >= 0) {
			int endIndex = data->_sizeIndexEnd[l];
			double minimumIntersection = minIntersection[l];
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
						arrayListtype *dataFP = dataBinary_getFingerPrint (data,   fpIndex);
						int intersection = c;
						if (k) {
							intersection += dataBinary_getIntersection_fast(specialQuery_fp, dataFP);
						}
						//int intersection = dataBinary_getIntersection_fast(queryFP, dataFP);

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
		static int x;
		fprintf (stderr, "%d : comparing the 2 solutions : \n",x++);
		if (x == 49) {
			x++;
			x--;
		}
		arrayListtype *solutionLinear = new_arrayListtype(_ARRAY_BUCKET_SIZE_II_BF);
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
			//sort the sloution
			minHeapPlain_Type *minHeap = new_minHeapPlain_Type(sSize+2);
			arrayList_Iterator_type *iteratorI = arrayList_getIterator(solution);
			while (iteratorI->nextAvailable) {
				u_long val = arrayList_IteratorGetNext(iteratorI);
				minHeapPlain_add(minHeap, val, 0);
			}
			arrayList_IteratorFree(iteratorI);

			arrayList_Iterator_type *iteratorL = arrayList_getIterator(solutionLinear);

			for (int i=0; i<sSize; i++) {
				int a = arrayList_IteratorGetNext(iteratorL);
				u_long b;
				int count;
				if (minHeapPlain_popMultiple(minHeap, &b, &count)) {
					helper_err(__FILE__, __LINE__, "not enough elements");
				} if (count > 1) {
					printf ("molid : %lu, occurs %d times\n", b, count);
				}

				if (a != b) {
					printf ("%d,%lu; ",a,b);
					helper_err(__FILE__, __LINE__, "Solutions do not match");
				}
			}
			arrayList_IteratorFree(iteratorL);
			minHeapPlain_free(minHeap);
		}

		free_arrayListtype(solutionLinear);
	}

	free (candidates);

	return (solution);
}

void _addFPIndex2moleculeList_BF (_moleculeList_BF *mList, int fpIndex, int numFeatures) {
	arrayList_flexible_put(&mList->arr, fpIndex);

	int index = mList->arr->size-1;

	if (mList->minNumFeatures > numFeatures)
		mList->minNumFeatures = numFeatures;
	if (mList->maxNumFeatures < numFeatures)
		mList->maxNumFeatures = numFeatures;

	//adjust the start and end index of this size
	if (! mList->_sizeIndexValid[numFeatures]) {
		mList->_sizeIndexValid[numFeatures] = 1;
		mList->_sizeIndexStart_1[numFeatures]=index;
	}
	mList->_sizeIndexEnd_1[numFeatures]=index;
}

workerFunctions_type invertedIndex_BF_getWorkerFunction () {
	workerFunctions_type worker;

	worker.dataDistribution = _dataDistribution_BF;
	worker.init_index 	= _init_index_BF;
	worker.free_index	= _free_index_BF;
	worker.runRange		= _runRangeHelper_BF_1;
	worker.runTopK		= _runTopKHelper_BF_1;
	return worker;
}
