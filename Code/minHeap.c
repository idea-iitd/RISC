/*
 * maxHeap.c
 *
 *  Created on: 02-Dec-2018
 *      Author: jithin
 */

#include "minHeap.h"

minHeaptype* newMinHeaptype (int k, void (*print)(void *, double, int printHeader)) {
	minHeaptype *heap = (minHeaptype *) helper_malloc(sizeof(minHeaptype),
			__FILE__, __LINE__);

	heap->k 	= k;
	heap->size 	= 0;
	heap->heap 	= (_minheap_element *) helper_calloc(k+1, sizeof(_minheap_element),
			__FILE__, __LINE__);

	_minheap_element *arr = heap->heap;

	for (int i=0; i<=k; i++) {
		//arr[i].key = DBL_MAX;
		arr[i].key = 0;
	}

	heap->print = print;
	return heap;
}

/**
 * return 0 is added successfully
 */
int minHeap_addIfKeyHigher (minHeaptype *heap, double key, void *val) {
	int notAdded = 1;

	int s = heap->size;
	if (s < heap->k) {
		notAdded = 0;

		//insert the element
		heap->size++;

		_minheap_element *arr = heap->heap;

		arr[s].element 	= val;
		arr[s].key		= key;

		//move it to the correct position
		while (s>0) {
			int parent = (s-1)/2;

			if (parent < 0)
				break;

			if (arr[parent].key > arr[s].key) {
				//swap
				_minheap_element temp;

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
		_minheap_element *arr = heap->heap;

		int s = 0;
		if (arr[s].key >= key) {
			//no need to add
		} else {
			notAdded = 0;

			arr[s].key 		= key;
			arr[s].element 	= val;

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
					_minheap_element temp;

					temp 	= arr[s];
					arr[s]	= arr[a];
					arr[a] 	= temp;

					s = a;
				}
			}
		}
	}

	return notAdded;
}

void minHeap_print (minHeaptype *heap) {
	int size = heap->size;
	int printHeader = 1;
	heap->print (NULL, -1,printHeader);
	printHeader = 0;

	for (int i=0; i<size; i++) {
		void *element;
		double key = minHeap_pop(heap, &element);
		printf ("%3d) ",size-i);
		heap->print (element, key,printHeader);
		free (element);
	}
}

double minHeap_peekKey (minHeaptype *heap) {
	return heap->heap[0].key;
}

void minHeap_free (minHeaptype *heap) {
	for (int i=0; i<heap->size; i++) {
		free (heap->heap->element);
	}
	free (heap->heap);
	free (heap);
}

/**
 * return the key
 */
double minHeap_pop (minHeaptype *heap, void **element) {
	_minheap_element *arr = heap->heap;
	int len = heap->size;
	double key;

	if (len <=0) {
		key = -1;
		*element = NULL;
	} else {
		key 	= arr[0].key;
		*element= arr[0].element;

		//re heapify
		len--;
		heap->size = len;
		if (len > 0) {
			arr[0].key 		= arr[len].key;
			arr[0].element 	= arr[len].element;

			int s = 0;

			int k = heap->size;

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
					_minheap_element temp;

					temp 	= arr[s];
					arr[s]	= arr[a];
					arr[a] 	= temp;

					s = a;
				}
			}
		}
	}

	return key;
}
