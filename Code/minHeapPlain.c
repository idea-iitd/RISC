/*
 * minHeapPlain.c
 *
 *  Created on: 10-Dec-2018
 *      Author: jithin
 */
#include "minHeapPlain.h"

minHeapPlain_Type* new_minHeapPlain_Type (u_long maxSize) {
	minHeapPlain_Type *heap = (minHeapPlain_Type *) helper_calloc (1, sizeof(minHeapPlain_Type),
			__FILE__, __LINE__);
	heap->maxSize = maxSize;
	heap->arr	  = (u_long *) helper_malloc(maxSize * sizeof(u_long), __FILE__, __LINE__);
	return heap;
}

/**
 * if the minheap is full,
 * create a new one and add unique elements.
 */
void minHeapPlain_add_WO_fail_unique (minHeapPlain_Type **heap_org, u_long val) {
	minHeapPlain_Type *heap = *heap_org;
	int notAdded = minHeapPlain_add(heap, val, 1);

	if (notAdded) {
		//if the heap is full, pop and populate unique element to a new heap
		u_long newSize = heap->size * 2;
		if (newSize < heap->size) {
			//we have overflowed
			newSize = -1; //this will set to the largest value
		}
		minHeapPlain_Type *tempHeap = new_minHeapPlain_Type(newSize);

		while (heap->size > 0) {
			u_long value;
			int count;
			minHeapPlain_popMultiple(heap, &value, &count);
			minHeapPlain_add		(tempHeap, value, 0);
		}
		minHeapPlain_add (tempHeap, val, 0);
		minHeapPlain_free(heap);
		*heap_org = tempHeap;
	}
}

/**
 * returb 0 if not added
 */
int minHeapPlain_add (minHeapPlain_Type *heap, u_long val, int softFail) {
	int notAdded = 0;

	if (heap->size >= heap->maxSize) {
		if(!softFail)
			helper_err(__FILE__, __LINE__, "not enough memory in heap");
		notAdded = 1;
	} else {
		u_long *arr= heap->arr;
		u_long s	= heap->size++;

		arr[s]	= val;
		//move the element to the correct possition
		while (s != 0) {
			u_long parent = (s-1)/2;

			if (arr[parent] > arr[s]) {
				//swap
				u_long temp 	= arr[s];
				arr[s] 		= arr[parent];
				arr[parent] = temp;
			} else {
				break;
			}

			s = parent;
		}
	}

	return notAdded;
}

int minHeapPlain_peek (minHeapPlain_Type *heap, u_long * val) {
	int notEnoughElement;

	if (heap->size) {
		*val 	= heap->arr[0];
		notEnoughElement = 0;
	} else {
		notEnoughElement = 1;
	}

	return notEnoughElement;
}

int minHeapPlain_pop (minHeapPlain_Type *heap, u_long * val) {
	int notEnoughElement;

	if (heap->size) {
		notEnoughElement = 0;
		heap->size--;
		u_long size = heap->size;

		u_long *arr= heap->arr;
		*val 	= arr[0];
		u_long p = 0;

		//heapify
		arr[p] = arr[size];

		while (p < size) {
			u_long lIndex = p*2 + 1;
			u_long rIndex = lIndex+1;

			u_long c; //candidate

			if (lIndex >= size) {
				//we are done
				break;
			} else if (rIndex>= size) {
				c = lIndex;
			} else {
				u_long l = arr[lIndex];
				u_long r = arr[rIndex];

				if (l < r) {
					c = lIndex;
				} else {
					c = rIndex;
				}
			}

			if (arr[p] > arr[c]) {
				//swap
				u_long temp= arr[p];
				arr[p] 	= arr[c];
				arr[c] 	= temp;

				p = c;
			} else {
				//we are done
				break;
			}
		}
	} else {
		notEnoughElement = 1;
	}

	return notEnoughElement;
}

int minHeapPlain_popMultiple (minHeapPlain_Type *heap, u_long *value, int *count) {
	int notEnoughElement;

	u_long val;
	int c = 1;
	notEnoughElement = minHeapPlain_pop(heap, &val);

	if (notEnoughElement) {
		//nothing to do
	} else {
		u_long valCurrent;

		while (1) {
			int failed = minHeapPlain_peek(heap, &valCurrent);
			if (!failed) {
				if (valCurrent == val) {
					minHeapPlain_pop(heap, &valCurrent);
					c++;
				} else {
					break;
				}
			} else {
				break;
			}
		}

		*value = val;
		*count = c;
	}

	return notEnoughElement;
}

void minHeapPlain_free (minHeapPlain_Type *heap) {
	free (heap->arr);
	free (heap);
}
