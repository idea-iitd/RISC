/*
 * arrayListInt.h
 *
 *  Created on: 01-Dec-2018
 *      Author: jithin
 */

#ifndef ARRAYLISTINT_H_
#define ARRAYLISTINT_H_

#include "global.h"
#include "helper.h"

struct _array_list_bucket {
	u_long *data;
	int size;
	struct _array_list_bucket *next;
	struct _array_list_bucket *prev;
};

typedef struct {
	u_long size; //number of elements
	unsigned int bucketMaxLen;
	struct _array_list_bucket *head;
	struct _array_list_bucket *tail;
}arrayListtype;

typedef struct {
	int nextIndex; 					//current index
	unsigned int _numBucketsBefore; // number of buckets before this
	u_long nextAvailable; 				//number of more elements available
	struct _array_list_bucket *bucket;
	int bucketMaxLen;
}arrayList_Iterator_type;


arrayListtype* new_arrayListtype (unsigned int bucketMaxLen);
u_long arrayList_get (arrayListtype *arr, u_long index);
void arrayList_put (arrayListtype *arr, u_long val);
void arrayList_iterateAndWork 	(arrayListtype *arr, u_long startIndex, u_long endIndex, void (*worker) (u_long));
void arrayList_iterateAndWork_2 (arrayListtype *arr, u_long startIndex, u_long endIndex, void (*worker) (int*, int), int *dest);
void arrayList_iterateAndWork_3 (arrayListtype *arr, u_long startIndex, u_long endIndex, int *dest);

void free_arrayListtype(arrayListtype *arr);
void arrayList_IteratorMoveTo		(arrayList_Iterator_type *iterator, u_long val);
void arrayList_IteratorMoveTo_Binary(arrayList_Iterator_type *iterator, u_long val);
int  arrayList_IteratorFind			(arrayList_Iterator_type *iterator, u_long val, int *found);
u_long arrayList_IteratorFind_Binary(arrayList_Iterator_type *iterator, u_long val, int *found);

//iterator
arrayList_Iterator_type* 	arrayList_getIterator 	(arrayListtype *arr);
arrayList_Iterator_type* 	arrayList_getIterator_2 (arrayListtype *arr, int startindex);
struct _array_list_bucket* 	arrayList_getBucket 	(arrayListtype *arr, int startIndex);
u_long  arrayList_IteratorGetNext 	(arrayList_Iterator_type *iterator);
void arrayList_IteratorFree 	(arrayList_Iterator_type *iterator);
void arrayList_resetIterator 	(arrayList_Iterator_type *iterator, arrayListtype *arr);

#endif /* ARRAYLISTINT_H_ */
