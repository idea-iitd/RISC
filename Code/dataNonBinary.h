/*
 * dataNonBinary.h
 *
 *  Created on: 30-Dec-2018
 *      Author: jithin
 */

#ifndef DATANONBINARY_H_
#define DATANONBINARY_H_

#include "global.h"
#include "options.h"
#include "arrayListInt.h"
#include "arrayListInt_Flexible.h"
#include "minHeapPlain.h"
#include "minHeapInt.h"

#include <limits.h>

typedef struct {
	u_long cummulativeFeatureVal;
	arrayListtype_flexible *features;
	arrayListtype_flexible *values;
}MolNonBinary;


typedef struct {
	u_long _maxSize; 				// the amount of elements for which memory is allocated
	u_long size; 					// the #elements stored
	MolNonBinary **_fpArray;		//NOT TO BE USED DIRECTLY, array of fingerprints,
	char *(*ids);					//array of fingerprint id
	u_long numFeatures;				//in the entire dataset
	u_long maxNumFeatures;			//in any molecule
	u_long minLen;
	u_long maxLen;
	int *_sizeIndexStart;
	int *_sizeIndexEnd;
	int *_cvMin; 					//Cumulative value max and min within each len
	int *_cvMax;
	int maxCV;
	int _isSorted;
	int maxNumOfMoleculesOfAnySize;
}DataNonBinary;


typedef struct {
	void  (*dataDistribution)(void *index, char *fname);			//an experimental function created to analyse results
	void* (*init_index)	(DataNonBinary *data);	//init index
	void  (*free_index)	(void *index);
	minHeapInttype* (*runTopK) 	(void *index, MolNonBinary *queryFP, int k, double *pruned, double *unPruned);
	arrayListtype*  (*runRange)	(void *index, MolNonBinary *queryFP, double r, double *pruned, double *unPruned);
}workerFunctions_NB_type;

void dataNonBinary_addMolecule  (DataNonBinary *data, char *ch, arrayList_Iterator_type *iterator, arrayListtype *featureIDs);
void dataNonBinary_sort 		(DataNonBinary *data);

void dataNonBinary_foldAndWrite2File  (char *buffer, FILE *fp, u_long fold, int molID);

double 	dataNonBinary_tanimotoSimilarity (MolNonBinary *queryFP, MolNonBinary *dataFP);
/*
int 	dataNonBinary_getIntersection 	 (arrayListtype *queryFP, arrayListtype *dataFP);
int 	dataNonBinary_getUnion 		     (arrayListtype *queryFP, arrayListtype *dataFP);
*/

MolNonBinary* 	dataNonBinary_getFingerPrint (DataNonBinary *data, u_long index);
DataNonBinary* 	new_dataNonBinary 			 (u_long maxSize, u_long numFeatures);

#endif /* DATANONBINARY_H_ */
