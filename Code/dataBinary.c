/*
 * dataBinary.c
 *
 *  Created on: 16-Dec-2018
 *      Author: jithin
 */

#include "dataBinary.h"


void dataBinary_addMolecule (DataBinary *data, minHeapPlain_Type *minHeap) {
	arrayListtype *arr = new_arrayListtype(100);

	if (minHeap->size > data->maxNumFeatures)
		data->maxNumFeatures = minHeap->size;

	while (minHeap->size) {
		u_long val;
		int count;

		minHeapPlain_popMultiple(minHeap, &val, &count);

		if (count > 1) {
			helper_err(__FILE__, __LINE__, "invalid feature values");
		}

		arrayList_put(arr, val);
	}

	data->_fpArray[data->size] = arr;
	data->molIDsORG[data->size] = data->size+1;
	data->size++;
}

DataBinary* new_dataBinary (int maxSize, int numFeatures) {	//allocate memory
	DataBinary *data = (DataBinary *)    helper_calloc(1,       sizeof(DataBinary),      __FILE__, __LINE__);
	data->_fpArray   = (arrayListtype **)helper_calloc(maxSize, sizeof(arrayListtype *), __FILE__, __LINE__);
	data->molIDsORG  = (arrayListtype  *)helper_calloc(maxSize, sizeof(u_long), __FILE__, __LINE__);

	data->_maxSize  	= maxSize;
	data->numFeatures 	= numFeatures;
	return data;
}

arrayListtype* dataBinary_getFingerPrint(DataBinary *data, int index) {
	return data->_fpArray[index];
}

void _dataBinary_swapFingerprints (DataBinary *data, int a, int b) {
	arrayListtype *fp_a = dataBinary_getFingerPrint(data, a);
	arrayListtype *fp_b = dataBinary_getFingerPrint(data, b);

	data->_fpArray[a] = fp_b;
	data->_fpArray[b] = fp_a;

	//swap the molIDs too
	u_long temp = data->molIDsORG[a];
	data->molIDsORG[a] = data->molIDsORG[b];
	data->molIDsORG[b] = temp;
}

/**
 * Find the similarity
 */
double dataBinary_tanimotoSimilarity (arrayListtype *queryFP, arrayListtype *dataFP) {
	int n=dataBinary_getIntersection(queryFP, dataFP);
	int d=queryFP->size + dataFP->size - n;

	double sim;
	if (d==0)
		sim = 0;
	else
		sim = (double)n/d;

	return sim;
}

double dataBinary_tanimotoSimilarity_fast (arrayListtype *queryFP, arrayListtype *dataFP) {
	int n=dataBinary_getIntersection_fast(queryFP, dataFP);
	int d=queryFP->size + dataFP->size - n;

	double sim;
	if (d==0)
		sim = 0;
	else
		sim = (double)n/d;

	return sim;
}

/**
 * Find the intersection
 */

int dataBinary_getIntersection (arrayListtype *queryFP, arrayListtype *dataFP) {
	int n=-1;

	arrayList_Iterator_type *iterator1 = arrayList_getIterator(queryFP);
	arrayList_Iterator_type *iterator2 = arrayList_getIterator(dataFP);

	int f1 = 0;
	int f2 = 0;

	int done = 0;
	while (!done) {
		if (f1==f2) {
			n++;
			if (iterator1->nextAvailable) {
				f1 = arrayList_IteratorGetNext(iterator1);
				if (iterator2->nextAvailable) {
					f2 = arrayList_IteratorGetNext(iterator2);
				} else {
					//done = 1;
					break;
				}
			} else {
				//done = 1;
				break;
			}
		}

		while (f1 < f2) {
			if (iterator1->nextAvailable) {
				f1 = arrayList_IteratorGetNext(iterator1);
			} else {
				done = 1;
				break;
			}
		}

		if (done)
			break;

		while (f2 < f1) {
			if (iterator2->nextAvailable) {
				f2 = arrayList_IteratorGetNext(iterator2);
			} else {
				done = 1;
				break;
			}
		}
	}

	arrayList_IteratorFree(iterator1);
	arrayList_IteratorFree(iterator2);

	return n;
}

int dataBinary_getIntersection_fast (arrayListtype *queryFP, arrayListtype *dataFP) {
	int n=-1;

	arrayList_Iterator_type *iterator1 = arrayList_getIterator(queryFP);
	arrayList_Iterator_type *iterator2 = arrayList_getIterator(dataFP);

	int f1 = 0;
	int f2 = 0;

	int done = 0;
	while (!done) {
		if (f1==f2) {
			n++;
			if (iterator1->nextAvailable) {
				f1 = arrayList_IteratorGetNext(iterator1);
				if (iterator2->nextAvailable) {
					f2 = arrayList_IteratorGetNext(iterator2);
				} else {
					//done = 1;
					break;
				}
			} else {
				//done = 1;
				break;
			}
		}

		if (f1 < f2) {
			arrayList_IteratorMoveTo_Binary(iterator1, f2);
			//arrayList_IteratorMoveTo(iterator1, f2);
			if (iterator1->nextAvailable) {
				f1 = arrayList_IteratorGetNext(iterator1);
			} else {
				done = 1;
				break;
			}
		}

		if (done)
			break;

		if (f2 < f1) {
			arrayList_IteratorMoveTo_Binary(iterator2, f1);
			//arrayList_IteratorMoveTo(iterator2, f1);
			if (iterator2->nextAvailable) {
				f2 = arrayList_IteratorGetNext(iterator2);
			} else {
				done = 1;
				break;
			}
		}
	}

	arrayList_IteratorFree(iterator1);
	arrayList_IteratorFree(iterator2);

	return n;
}

/**
 * Find the union
 */
int dataBinary_getUnion (arrayListtype *queryFP, arrayListtype *dataFP) {
	int n=queryFP->size + dataFP->size - dataBinary_getIntersection(queryFP, dataFP);
	return n;
}

/**
 * sort the data based on the len
 */
void dataBinary_sort (DataBinary *data) {
	int maxNumFeatures = data->maxNumFeatures;
	//allocate memory
	data->_sizeIndexStart= (int *) helper_malloc((size_t)(maxNumFeatures + 1) * sizeof(int),  __FILE__, __LINE__);
	data->_sizeIndexEnd	 = (int *) helper_malloc((size_t)(maxNumFeatures + 1) * sizeof(int),  __FILE__, __LINE__);

	for (int i=0; i<=maxNumFeatures; i++) {
		data->_sizeIndexStart[i]=-1;
	}

	//maintain the data sorted by len
	int *lens = (int *) helper_malloc((size_t)(data->size) * sizeof(int),  __FILE__, __LINE__);

	//populate the sizes
	int sizeData 	= data->size;

	int maxLen = 0;
	for (int i=0; i<sizeData; i++) {
		arrayListtype *fp 	= dataBinary_getFingerPrint(data, i);
		int numFeatures		= fp->size;
		lens[i] = numFeatures;

		if (maxLen < numFeatures)
			maxLen = numFeatures;
	}

	//sort
	int currentIndex = 0;
	for (int l =1; l<maxLen; l++) {
		for (int k = currentIndex; k<sizeData; k++) {
			if (lens[k]==l) {
				//swap lengths
				lens[k] = lens[currentIndex];
				lens[currentIndex] = l;

				//swap the fingerprints
				_dataBinary_swapFingerprints (data, k, currentIndex);

				currentIndex++;
			}
		}
	}

	int maxSize=0;

	for (int i=0; i<sizeData; i++) {
		int l = lens[i];

		data->_sizeIndexStart[l]=i;

		while ((l == lens[i]) && (i<sizeData))
			i++;
		i--;
		data->_sizeIndexEnd[l]=i;

		int size = data->_sizeIndexEnd[l] - data->_sizeIndexStart[l] + 1;
		if (size > maxSize)
			maxSize = size;
	}
	data->_isSorted = 1;
	data->maxNumOfMoleculesOfAnySize = maxSize;

	free (lens);
}

void dataBinary_printFeatures(arrayListtype *dataFP) {
	arrayList_Iterator_type *iter = arrayList_getIterator(dataFP);

	while (iter->nextAvailable) {
		int feature = arrayList_IteratorGetNext(iter);
		printf("%d, ",feature);
	}
	printf("\n");

	arrayList_IteratorFree(iter);
}
