/*
 * maxHeapInt.h
 *
 *  Created on: 02-Dec-2018
 *      Author: jithin
 */

#ifndef MINHEAP_H_
#define MINHEAP_H_

#include "global.h"
#include "helper.h"
#include <float.h>

typedef struct {
	void *element;
	double key;
}_minheap_element;


typedef struct {
	_minheap_element *heap;
	int k; //the size of heap is fixed at k
	int size;
	void (*print)(void *, double, int printHeader);
}minHeaptype;

minHeaptype* newMinHeaptype 	(int k, void (*print)(void *, double, int printHeader));
int  minHeap_addIfKeyHigher 	(minHeaptype *heap, double key, void *element);
void minHeap_print 				(minHeaptype *heap);
double minHeap_pop 				(minHeaptype *heap, void **element);
double minHeap_peekKey 			(minHeaptype *heap);
void   minHeap_free 			(minHeaptype *heap);
#endif /* MINHEAPINT_H_ */
