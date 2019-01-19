/*
 * maxHeapInt.h
 *
 *  Created on: 02-Dec-2018
 *      Author: jithin
 */

#ifndef MINHEAPINT_H_
#define MINHEAPINT_H_

#include "global.h"
#include "helper.h"
#include <float.h>

typedef struct {
	int val;
	double key;
}_intheap_element;


typedef struct {
	_intheap_element *heap;
	int k; //the size of heap is fixed at k
	int size;
	double defaultKey; //the default value returned if size < k
}minHeapInttype;

minHeapInttype* newMinHeapInttype 			 (int k);
minHeapInttype* newMinHeapInttype_withDefault(int k, double defaultKey);

int  minHeapInt_addIfKeyHigher 		(minHeapInttype *heap, double key, int val);

void minHeapInt_getValues 			(minHeapInttype *heap, int *val);
void minHeapInt_getKeys 			(minHeapInttype *heap, double *val);
void minHeapInt_print 				(minHeapInttype *heap);
void minHeapInt_free 				(minHeapInttype *heap);

double minHeapInt_peekKey 			(minHeapInttype *heap);
double minHeapInt_peekKey_or_default(minHeapInttype *heap);

#endif /* MINHEAPINT_H_ */
