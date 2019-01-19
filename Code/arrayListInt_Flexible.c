/*
 * arrayListInt_Flexible.c
 *
 *  Created on: 24-Dec-2018
 *      Author: jithin
 */

#include "arrayListInt_Flexible.h"


arrayListtype_flexible* new_arrayListtype_flexible (u_long len, u_long upperLimit) {
	arrayListtype_flexible *list = (arrayListtype_flexible *) helper_calloc(1, sizeof(arrayListtype_flexible),
				__FILE__, __LINE__);
	list->data			= (u_long *) helper_malloc (len * sizeof(u_long), __FILE__, __LINE__);
	list->maxSize 		= len;
	list->upperLimit 	= upperLimit;
	return list;
}

u_long  arrayList_flexible_get (arrayListtype_flexible *arr,  u_long index) {
	return arr->data[index];
}

void arrayList_flexible_put (arrayListtype_flexible **list_para, u_long val) {
	arrayListtype_flexible *list = *list_para;

	if (list->size < list->maxSize) {
		//we are fine
		list->data[list->size++] = val;
	} else {
		//we need to allcate more memory and copy the current contents
		u_long newSize = list->maxSize * 2;
		if (list->upperLimit > 0) {
			if (newSize > list->upperLimit)
				newSize = list->upperLimit;
		}
		arrayListtype_flexible *listNew = new_arrayListtype_flexible(newSize, list->upperLimit);
		listNew->size 		= list->size;
		listNew->upperLimit = list->upperLimit;
		memcpy(listNew->data, list->data, sizeof(u_long) * list->size);

		listNew->data[listNew->size++] = val;
		*list_para = listNew;
		free_arrayListtype_flexible(list);
	}
}

void arrayList_flexible_iterateAndWork_3(arrayListtype_flexible *list,  u_long startIndex, u_long endIndex, int *dest) {
	u_long *arr = list->data;

	for (u_long i=startIndex; i<= endIndex; i++) {
		u_long val = arr[i];
		dest[val] = 1;
	}
}

void free_arrayListtype_flexible (arrayListtype_flexible *list) {
	free (list->data);
	free (list);
}

void arrayList_flexible_IteratorFree (arrayList_flexible_Iterator_type *iterator) {
	free (iterator);
}


/**
 * This assumes the array list is sorted
 */
void arrayList_flexible_IteratorMoveTo(arrayList_flexible_Iterator_type *iterator, u_long val) {
	while (iterator->nextAvailable) {
		u_long v = iterator->data[iterator->nextIndex];

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
void arrayList_flexible_IteratorMoveTo_Binary(arrayList_flexible_Iterator_type *iterator, u_long val) {
	if (iterator->nextAvailable){
		//we do binary search between 0 and lastindex
		u_long *arr = iterator->data;
		int l=iterator->nextIndex, r=iterator->nextIndex + iterator->nextAvailable - 1;
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
			u_long v = iterator->data[iterator->nextIndex];

			if (v>=val)
				break;
			else {
				iterator->nextAvailable--;
				iterator->nextIndex++;
			}
		}
	}
}


/**
 * The function does not do error check,
 * you should check for next element availability before invoking this
 */
u_long arrayList_flexible_IteratorGetNext (arrayList_flexible_Iterator_type *iterator) {
	iterator->nextAvailable--;
	return iterator->data[iterator->nextIndex++];
}

/**
 * if NULL is passed the iterator will be initalized without any value
 */
arrayList_flexible_Iterator_type* arrayList_flexible_getIterator (arrayListtype_flexible *arr) {
	arrayList_flexible_Iterator_type *iter = (arrayList_flexible_Iterator_type *) helper_calloc
			(1, sizeof(arrayList_flexible_Iterator_type),__FILE__, __LINE__);
	iter->nextIndex	=0;
	if (arr) {
		iter->data 			= arr->data;
		iter->nextAvailable = arr->size;
	}
	return iter;
}

int arrayList_flexible_IteratorFind(arrayList_flexible_Iterator_type *iterator, u_long val, int *found) {
	int notEnoughElements=1;
	*found = 0;

	while (iterator->nextAvailable) {
		u_long v = iterator->data[iterator->nextIndex];

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
u_long arrayList_flexible_IteratorFind_Binary(arrayList_flexible_Iterator_type *iterator, u_long val, int *found) {
	*found = 0;
	int ret = -1;

	if (iterator->nextAvailable){
		//we do binary search between first index and lastindex
		u_long *arr = iterator->data;
		int l=iterator->nextIndex, r=iterator->nextAvailable + iterator->nextIndex - 1;
		int c=l;
		while (l <= r) {
			c=(l+r)/2;

			u_long v = arr[c];
			if (v == val) {
				//we have found the element
				ret = c;
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

		if (l > r) {
			//we have not found the element, we shall move the
			//iterator back so that subsequent calls do not fail
			while (c > iterator->nextIndex) {
				u_long v = arr[c];
				if (v > val)
					c--;
				else {
					iterator->nextAvailable -= (c - iterator->nextIndex);
					iterator->nextIndex = c;
					break;
				}
			}
		}
	}

	return ret;
}

