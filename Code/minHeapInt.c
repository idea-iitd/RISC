/*
 * maxHeapInt.c
 *
 *  Created on: 02-Dec-2018
 *      Author: jithin
 */

#include "minHeapInt.h"

minHeapInttype* newMinHeapInttype (int k) {
	return newMinHeapInttype_withDefault(k,0);
}

minHeapInttype* newMinHeapInttype_withDefault(int k, double defaultKey) {
	minHeapInttype *heap = (minHeapInttype *) helper_malloc(sizeof(minHeapInttype),
			__FILE__, __LINE__);

	heap->k 	= k;
	heap->size 	= 0;
	heap->heap 	= (_intheap_element *) helper_calloc(k+1, sizeof(_intheap_element),
			__FILE__, __LINE__);

	/*
	_intheap_element *arr = heap->heap;

	for (int i=0; i<=k; i++) {
		//arr[i].key = DBL_MAX;
		arr[i].key = 0;
	}
	*/
	heap->defaultKey = defaultKey;

	return heap;
}

/**
 * return 0 is added successfully
 */
int minHeapInt_addIfKeyHigher (minHeapInttype *heap, double key, int val) {
	int notAdded = 1;

	if (key > 0) {
		int s = heap->size;
		if (s < heap->k) {
			notAdded = 0;
			_intheap_element element;
			element.key = key;
			element.val = val;

			//insert the element
			heap->size++;

			_intheap_element *arr = heap->heap;

			arr[s] = element;

			//move it to the correct position
			while (s>0) {
				int parent = (s-1)/2;

				if (parent < 0)
					break;

				if (arr[parent].key > arr[s].key) {
					//swap
					_intheap_element temp;

					temp 	= arr[s];
					arr[s]	= arr[parent];
					arr[parent] = temp;

					//move up the tree
					s=parent;
				} else {
					break;
				}
			}
		} else {
			_intheap_element *arr = heap->heap;

			int s = 0;
			if (arr[s].key >= key) {
				//no need to add
			} else {
				notAdded = 0;

				arr[s].key = key;
				arr[s].val = val;

				int k = heap->k;

				while (1) {
					//move the node to the correct position
					int l = s*2 + 1; //left child
					int r = l + 1; // right child

					//which of these children to be considered

					int a; //smaller key

					if (l >= k) {
						//we are at the end of heap
						break;
					} else if (r >=k) {
						//left is fine but not right
						a = l;
					} else if (arr[l].key > arr[r].key) {
						a = r;
					} else {
						a = l;
					}

					if (arr[s].key < arr[a].key) {
						//nothing to be done we are at the correct position
						break;
					} else {
						//swap
						_intheap_element temp;

						temp 	= arr[s];
						arr[s]	= arr[a];
						arr[a] 	= temp;

						s = a;
					}
				}
			}
		}
	}
	return notAdded;
}


void minHeapInt_getValues (minHeapInttype *heap, int *val) {
	int k = heap->k;

	_intheap_element *arr = heap->heap;

	for (int i=0; i<k; i++) {
		val[i] = arr[i].val;
	}
}

void minHeapInt_getKeys (minHeapInttype *heap, double *val) {
	int k = heap->k;

	_intheap_element *arr = heap->heap;

	for (int i=0; i<k; i++) {
		val[i] = arr[i].key;
	}
}

void _minHeapInt_print (_intheap_element *arr, int i, int k) {
	int left 	= i*2 +1;
	int right 	= left+1;
	printf ("%d, %f\n",arr[i].val, arr[i].key);

	double l,r;

	if (left < k && right < k) {
		l = arr[left].key;
		r = arr[right].key;

		if (l<r) {
			_minHeapInt_print (arr, left, k);
			_minHeapInt_print (arr, right, k);
		} else {
			_minHeapInt_print (arr, right, k);
			_minHeapInt_print (arr, left, k);
		}
	} else if (left < k) {
		_minHeapInt_print (arr, left, k);
	} else if (right < k) {
		_minHeapInt_print (arr, right, k);
	}
}

void minHeapInt_print (minHeapInttype *heap) {
	int k = heap->k;

	_intheap_element *arr = heap->heap;

	printf ("val, key\n");

	_minHeapInt_print(arr, 0, k);
	//for (int i=0; i<k; i++) {
	//	printf ("%d, %f\n",arr[i].val, arr[i].key);
	//}
}

double minHeapInt_peekKey (minHeapInttype *heap) {
	return heap->heap[0].key;
}

double minHeapInt_peekKey_or_default(minHeapInttype *heap) {
	if (heap->size < heap->k)
		return heap->defaultKey;
	else
		return heap->heap[0].key;
}
void minHeapInt_free (minHeapInttype *heap) {
	free (heap->heap);
	free (heap);
}
