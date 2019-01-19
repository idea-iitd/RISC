/*
 * minHeapPlain.h
 *
 *  Created on: 10-Dec-2018
 *      Author: jithin
 */

#ifndef MINHEAPPLAIN_H_
#define MINHEAPPLAIN_H_

#include "helper.h"

typedef struct {
	u_long maxSize;
	u_long size;
	u_long *arr;
}minHeapPlain_Type;

minHeapPlain_Type* new_minHeapPlain_Type (u_long maxSize);

int minHeapPlain_add 		(minHeapPlain_Type *heap, u_long  val, int softFail);
int minHeapPlain_pop 		(minHeapPlain_Type *heap, u_long *val);
int minHeapPlain_peek 		(minHeapPlain_Type *heap, u_long *val);
int minHeapPlain_popMultiple(minHeapPlain_Type *heap, u_long *val, int *count);
void minHeapPlain_free 		(minHeapPlain_Type *heap);

void minHeapPlain_add_WO_fail_unique (minHeapPlain_Type **heap_org, u_long val);

#endif /* MINHEAPPLAIN_H_ */
