/*
 * dataBinary.h
 *
 *  Created on: 16-Dec-2018
 *      Author: jithin
 */

#ifndef DATABINARY_H_
#define DATABINARY_H_

#include "global.h"
#include "options.h"
#include "arrayListInt.h"
#include "minHeapPlain.h"
#include "minHeapInt.h"

typedef struct {
	int _maxSize; 					// the amount of elements for which memory is allocated
	int size; 						// the #elements stored
	arrayListtype **_fpArray;		//NOT TO BE USED DIRECTLY, array of fingerprints,
	u_long *molIDsORG;				// the original molIDS
	int *_sizeIndexStart;			// stores the starting index of molecules of given size
	int *_sizeIndexEnd;				// stores the ending index(inclusive) of molecules of given size
	int maxNumOfMoleculesOfAnySize;
	int _isSorted;
	int numFeatures;				//in the entire dataset
	int maxNumFeatures;				//in any molecule
}DataBinary;


typedef struct {
	void  (*dataDistribution)(void *index, char *fname);			//an experimental function created to analyse results
	void* (*init_index)	(DataBinary *data);	//init index
	void  (*free_index)	(void *index);
	minHeapInttype* (*runTopK) 	(void *index, arrayListtype *queryFP, int k, double *pruned, double *unPruned);
	arrayListtype*  (*runRange)	(void *index, arrayListtype *queryFP, double r, double *pruned, double *unPruned);
}workerFunctions_type;


void dataBinary_addMolecule (DataBinary *data, minHeapPlain_Type *minHeap);
void dataBinary_sort 		(DataBinary *data);

void 	dataBinary_printFeatures	 (arrayListtype *dataFP);
double 	dataBinary_tanimotoSimilarity(arrayListtype *queryFP, arrayListtype *dataFP);
double dataBinary_tanimotoSimilarity_fast (arrayListtype *queryFP, arrayListtype *dataFP);
int 	dataBinary_getIntersection 	 (arrayListtype *queryFP, arrayListtype *dataFP);
int 	dataBinary_getIntersection_fast (arrayListtype *queryFP, arrayListtype *dataFP);
int 	dataBinary_getUnion 		 (arrayListtype *queryFP, arrayListtype *dataFP);


arrayListtype* 	dataBinary_getFingerPrint(DataBinary *data, int index);
DataBinary* 	new_dataBinary 			 (int maxSize, int numFeatures);

#endif /* DATABINARY_H_ */
