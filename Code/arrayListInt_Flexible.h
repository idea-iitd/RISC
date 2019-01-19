/*
 * arrayListInt_Flexible.h
 *
 *  Created on: 24-Dec-2018
 *      Author: jithin
 */

#ifndef ARRAYLISTINT_FLEXIBLE_H_
#define ARRAYLISTINT_FLEXIBLE_H_

#include "global.h"
#include "helper.h"

typedef struct {
	u_long size; //number of elements
	u_long maxSize;
	u_long upperLimit;
	u_long *data;
}arrayListtype_flexible;

typedef struct {
	int nextIndex; 					//current index
	u_long nextAvailable; 			//number of more elements available
	u_long *data;
}arrayList_flexible_Iterator_type;

arrayListtype_flexible* new_arrayListtype_flexible (u_long len, u_long upperLimit);

u_long  arrayList_flexible_get 			(arrayListtype_flexible *arr,  u_long index);
void arrayList_flexible_put 			(arrayListtype_flexible **arr, u_long val);
void arrayList_flexible_iterateAndWork_3(arrayListtype_flexible *arr,  u_long startIndex, u_long endIndex, int *dest);
void free_arrayListtype_flexible		(arrayListtype_flexible *arr);

//iterator
void 	arrayList_flexible_IteratorMoveTo		(arrayList_flexible_Iterator_type *iterator, u_long val);
void 	arrayList_flexible_IteratorMoveTo_Binary(arrayList_flexible_Iterator_type *iterator, u_long val);
int  	arrayList_flexible_IteratorFind			(arrayList_flexible_Iterator_type *iterator, u_long val, int *found);
u_long 	arrayList_flexible_IteratorFind_Binary	(arrayList_flexible_Iterator_type *iterator, u_long val, int *found);

u_long  arrayList_flexible_IteratorGetNext 		(arrayList_flexible_Iterator_type *iterator);
void    arrayList_flexible_IteratorFree 		(arrayList_flexible_Iterator_type *iterator);

arrayList_flexible_Iterator_type* arrayList_flexible_getIterator (arrayListtype_flexible *arr);

#endif /* ARRAYLISTINT_FLEXIBLE_H_ */
