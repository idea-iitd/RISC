/*
 * arrayListInt.c
 *
 *  Created on: 01-Dec-2018
 *      Author: jithin
 */

#include "arrayListInt.h"

struct _array_list_bucket* _new_array_list_bucket (unsigned int len) {
	struct _array_list_bucket *bucket = (struct _array_list_bucket *) helper_calloc
			(1, sizeof(struct _array_list_bucket), __FILE__, __LINE__);

	bucket->data = (u_long *) helper_calloc (len, sizeof(u_long), __FILE__, __LINE__);

	return bucket;
}

arrayListtype* new_arrayListtype (unsigned int bucketMaxLen) {
	arrayListtype *list = (arrayListtype *) helper_calloc(1, sizeof(arrayListtype),
			__FILE__, __LINE__);

	list->head 			= _new_array_list_bucket(bucketMaxLen);
	list->tail 			= list->head;
	list->bucketMaxLen 	= bucketMaxLen;

	return list;
}

u_long arrayList_get (arrayListtype *arr, u_long index) {
	if (index > arr->size)
		helper_err(__FILE__, __LINE__, "index out of bound");

	//get to the correct bucket
	int skip = index/arr->bucketMaxLen;

	struct _array_list_bucket *bucket = arr->head;

	while (skip--) {
		bucket = bucket->next;
	}

	int i = index % arr->bucketMaxLen;

	return bucket->data[i];
}

void arrayList_put (arrayListtype *arr, u_long val) {
	struct _array_list_bucket *bucket = arr->tail;

	if (bucket->size >= arr->bucketMaxLen) {
		struct _array_list_bucket *b = _new_array_list_bucket(arr->bucketMaxLen);
		b->prev 	= bucket;
		bucket->next= b;
		bucket 		= b;
		arr->tail 	= bucket;
	}

	bucket->data[bucket->size++]=val;
	arr->size++;
}

void arrayList_iterateAndWork (arrayListtype *arr, u_long startIndex, u_long endIndex, void (*worker) (u_long)) {
	if (endIndex > arr->size || startIndex < 0)
		helper_err(__FILE__, __LINE__, "index out of bound");

	//get to the correct bucket
	int skip = startIndex/arr->bucketMaxLen;

	struct _array_list_bucket *bucket = arr->head;

	while (skip--) {
		bucket = bucket->next;
	}

	//iterate
	int numElements = endIndex-startIndex+1;
	int index = startIndex % arr->bucketMaxLen;

	while (numElements--) {
		if (index >= arr->bucketMaxLen) {
			index = 0;
			bucket = bucket->next;
		}

		u_long val = bucket->data[index++];

		(*worker) (val);
	}
}

void arrayList_iterateAndWork_2 (arrayListtype *arr, u_long startIndex, u_long endIndex,
			void (*worker) (int*, int), int *dest) {
	if (endIndex > arr->size || startIndex < 0) {
		char c[100];
		sprintf(c, "index out of bound startIndex:%lu, endIndex:%lu, arraySize:%lu",startIndex, endIndex, arr->size);
		helper_err(__FILE__, __LINE__, c);
	}

	//get to the correct bucket
	int skip = startIndex/arr->bucketMaxLen;

	struct _array_list_bucket *bucket = arr->head;

	while (skip--) {
		bucket = bucket->next;
	}

	//iterate
	int numElements = endIndex-startIndex+1;
	int index = startIndex % arr->bucketMaxLen;

	while (numElements--) {
		if (index >= arr->bucketMaxLen) {
			index = 0;
			bucket = bucket->next;
		}

		u_long val = bucket->data[index++];

		(*worker) (dest, val);
	}
}

void arrayList_iterateAndWork_3 (arrayListtype *arr, u_long startIndex, u_long endIndex, int *dest) {
	if (endIndex > arr->size || startIndex < 0) {
		char c[100];
		sprintf(c, "index out of bound startIndex:%lu, endIndex:%lu, arraySize:%lu",startIndex, endIndex, arr->size);
		helper_err(__FILE__, __LINE__, c);
	}

	//get to the correct bucket
	int skip = startIndex/arr->bucketMaxLen;

	struct _array_list_bucket *bucket = arr->head;

	while (skip--) {
		bucket = bucket->next;
	}

	//iterate
	int numElements = endIndex-startIndex+1;
	int index = startIndex % arr->bucketMaxLen;

	while (numElements--) {
		if (index >= arr->bucketMaxLen) {
			index = 0;
			bucket = bucket->next;
		}

		u_long val = bucket->data[index++];

		dest[val]=1;
	}
}

void free_arrayListtype(arrayListtype *arr) {
	struct _array_list_bucket *bucket=arr->head;

	while (bucket) {
		struct _array_list_bucket *temp = bucket->next;
		free (bucket->data);
		free (bucket);
		bucket = temp;
	}

	free(arr);
}


void arrayList_resetIterator (arrayList_Iterator_type *iter, arrayListtype *arr) {
	iter->nextIndex		=0;
	iter->_numBucketsBefore	=0;
	iter->bucketMaxLen 	= arr->bucketMaxLen;
	iter->bucket 		= arr->head;
	iter->nextAvailable = arr->size;
}

/**
 * if NULL is passed the iterator will be initalized without any value
 */
arrayList_Iterator_type* arrayList_getIterator (arrayListtype *arr) {
	arrayList_Iterator_type *iter = (arrayList_Iterator_type *) helper_calloc(1, sizeof(arrayList_Iterator_type),
			__FILE__, __LINE__);
	iter->nextIndex	=0;
	if (arr) {
		iter->bucketMaxLen 	= arr->bucketMaxLen;
		iter->bucket 		= arr->head;
		iter->nextAvailable = arr->size;
	}
	return iter;
}

arrayList_Iterator_type* arrayList_getIterator_2 (arrayListtype *arr, int startIndex) {
	arrayList_Iterator_type *iter = (arrayList_Iterator_type *) helper_malloc(sizeof(arrayList_Iterator_type),
			__FILE__, __LINE__);

	iter->nextAvailable = arr->size - startIndex;

	struct _array_list_bucket *bucket = arr->head;

	//got to correct bucket
	for  (;startIndex >= arr->bucketMaxLen; startIndex -= arr->bucketMaxLen)
		bucket = bucket->next;

	iter->nextIndex		= startIndex;
	iter->bucket 		= bucket;
	iter->bucketMaxLen 	= arr->bucketMaxLen;
	return iter;
}

struct _array_list_bucket* arrayList_getBucket (arrayListtype *arr, int startIndex) {
	struct _array_list_bucket *bucket = arr->head;

	//got to correct bucket
	for  (;startIndex >= arr->bucketMaxLen; startIndex -= arr->bucketMaxLen)
		bucket = bucket->next;
	return bucket;
}

/**
 * The function does not do error check,
 * you should check for next element availability before invoking this
 */
u_long arrayList_IteratorGetNext (arrayList_Iterator_type *iterator) {
	iterator->nextAvailable--;
	if (iterator->nextIndex >= iterator->bucketMaxLen){
		iterator->nextIndex = 0;
		iterator->bucket = iterator->bucket->next;
	}
	return iterator->bucket->data[iterator->nextIndex++];
}

void arrayList_IteratorFree (arrayList_Iterator_type *iterator) {
	free (iterator);
}

/**
 * This assumes the array list is sorted
 */
void arrayList_IteratorMoveTo(arrayList_Iterator_type *iterator, u_long val) {
	while (iterator->nextAvailable) {
		if (iterator->nextIndex >= iterator->bucketMaxLen){
			iterator->nextIndex = 0;
			iterator->bucket = iterator->bucket->next;
		}

		u_long v = iterator->bucket->data[iterator->nextIndex];

		if (v>=val)
			break;
		else {
			iterator->nextAvailable--;
			iterator->nextIndex++;
		}
	}
}

/**
 * This assumes the array list is sorted
 */
void arrayList_IteratorMoveTo_Binary(arrayList_Iterator_type *iterator, u_long val) {
	int bucketMaxLen = iterator->bucketMaxLen;

	int lastIndex; //in the bucket
	//goto the right bucket
	while (iterator->nextAvailable > 0) {
		int dip;

		if ((iterator->nextIndex + iterator->nextAvailable) > bucketMaxLen) {
			lastIndex 	= bucketMaxLen - 1;
			dip 		= bucketMaxLen - iterator->nextIndex;
		} else {
			lastIndex 	= iterator->nextIndex + iterator->nextAvailable - 1;
			dip 		= iterator->nextAvailable;
		}

		//is the last element lesser than greater than val
		if (iterator->bucket->data[lastIndex] >= val) {
			//this is the bucket to be explored
			break;
		}

		iterator->nextAvailable -= dip;
		iterator->bucket		 = iterator->bucket->next;
		iterator->nextIndex		 = 0;
		iterator->_numBucketsBefore++;
	}

	//search the bucket
	if (iterator->nextAvailable && iterator->bucket){
		//we do binary search between 0 and lastindex
		u_long *arr = iterator->bucket->data;
		int l=iterator->nextIndex, r=lastIndex;
		int c;
		int index=l;
		while (l <= r) {
			c=(l+r)/2;

			u_long v = arr[c];
			if (v == val) {
				//we have found the element
				index = c;
				break;
			} else if (v < val) {
				//explore second half
				index = c;
				l = c + 1;
			} else {
				//explore first half
				r = c - 1;
			}
		}

		//adjust the iterator
		iterator->nextAvailable -= (index - iterator->nextIndex);
		iterator->nextIndex 	 = index;

		while (iterator->nextAvailable) {
			if (iterator->nextIndex >= iterator->bucketMaxLen){
				iterator->nextIndex = 0;
				iterator->bucket = iterator->bucket->next;
			}

			u_long v = iterator->bucket->data[iterator->nextIndex];

			if (v>=val)
				break;
			else {
				iterator->nextAvailable--;
				iterator->nextIndex++;
			}
		}
	}
}

int arrayList_IteratorFind(arrayList_Iterator_type *iterator, u_long val, int *found) {
	int notEnoughElements=1;
	*found = 0;

	int bucketMaxLen = iterator->bucketMaxLen;

	//goto the right bucket
	while (iterator->nextAvailable > 0) {
		int lastIndex; //in the bucket
		int dip;

		if ((iterator->nextIndex + iterator->nextAvailable) > bucketMaxLen) {
			lastIndex 	= bucketMaxLen - 1;
			dip 		= bucketMaxLen - iterator->nextIndex;
		} else {
			lastIndex = iterator->nextIndex + iterator->nextAvailable - 1;
			dip 		= iterator->nextAvailable;
		}
		//is the las element lesser than greater than val
		if (iterator->bucket->data[lastIndex] >= val) {
			//this is the bucket to be explored
			break;
		}

		iterator->nextAvailable -= dip;
		iterator->bucket		 = iterator->bucket->next;
		iterator->nextIndex 	 = 0;
		iterator->_numBucketsBefore++;
	}

	while (iterator->nextAvailable) {
		if (iterator->nextIndex >= iterator->bucketMaxLen){
			iterator->nextIndex = 0;
			iterator->bucket = iterator->bucket->next;
		}

		u_long v = iterator->bucket->data[iterator->nextIndex];

		if (v>val) {
			notEnoughElements = 0;
			break;
		} else if (v==val) {
			notEnoughElements = 0;
			*found = 1;
			break;
		} else {
			iterator->nextIndex++;
			iterator->nextAvailable--;
		}
	}

	return notEnoughElements;
}

/*
 * return the index of the molecule
 */
u_long arrayList_IteratorFind_Binary(arrayList_Iterator_type *iterator, u_long val, int *found) {
	*found = 0;
	int ret = -1;

	int bucketMaxLen = iterator->bucketMaxLen;
	int lastIndex; //in the bucket

	//goto the right bucket
	while (iterator->nextAvailable > 0) {
		int dip;

		if ((iterator->nextIndex + iterator->nextAvailable) > bucketMaxLen) {
			lastIndex 	= bucketMaxLen - 1;
			dip 		= bucketMaxLen - iterator->nextIndex;
		} else {
			lastIndex = iterator->nextIndex + iterator->nextAvailable - 1;
			dip 		= iterator->nextAvailable;
		}
		//is the las element lesser than greater than val
		if (iterator->bucket->data[lastIndex] >= val) {
			//this is the bucket to be explored
			break;
		}

		iterator->nextAvailable -= dip;
		iterator->bucket		 = iterator->bucket->next;
		iterator->nextIndex 	 = 0;
		iterator->_numBucketsBefore++;
	}

	//search the bucket
	if (iterator->nextAvailable && iterator->bucket){
		//we do binary search between 0 and lastindex
		u_long *arr = iterator->bucket->data;
		int l=0, r=lastIndex;
		int c;
		while (l <= r) {
			c=(l+r)/2;

			u_long v = arr[c];
			if (v == val) {
				//we have found the element
				ret = iterator->_numBucketsBefore * bucketMaxLen + c;
				*found = 1;
				break;
			} else if (v < val) {
				//explore second half
				l = c + 1;
			} else {
				//explore first half
				r = c - 1;
			}
		}
	}

	return ret;
}
