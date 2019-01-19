/*
 * data.h
 *
 *  Created on: 24-Nov-2018
 *      Author: jithin
 */

#ifndef DATA_H_
#define DATA_H_

#include "global.h"
#include "helper.h"
#include "arrayListInt.h"
#include "minHeapInt.h"

typedef unsigned long fingerPrintIntType;
#define  NUMBITS_fingerPrintIntType (sizeof(fingerPrintIntType)*8)

typedef struct {
	int maxSize; 					// the amount of elements for which memory is allocated
	int size; 						// the #elements stored
	int fingerPrintLen;				//number of elements of type fingerPrintIntType in each fingerprint
	fingerPrintIntType *_fpArray;	//NOT TO BE USED DIRECTLY, array of fingerprints,
	char *(*ids);					//array of fingerprint id
	int *_sizeIndexStart;			// stores the starting index of molecules of given size
	int *_sizeIndexEnd;				// stores the ending index(inclusive) of molecules of given size
	int maxNumOfMoleculesOfAnySize;
	int _isSorted;
	int maxNumFeatures;				//in any molecule, valid only after sorting
	int minNumFeatures;				//in any molecule, valid only after sorting
}DataSmall;

typedef struct {
	void  (*dataDistribution)(void *index, char *fName);			//an experimental function created to analyse results
	void* (*init_index)	(DataSmall *data);	//init index
	void  (*free_index)	(void *index);
	minHeapInttype* (*runTopK) 	(void *index, fingerPrintIntType *queryFP, int k, double *pruned, double *unPruned);
	arrayListtype*  (*runRange)	(void *index, fingerPrintIntType *queryFP, double r, double *pruned, double *unPruned);
}workerFunctions_type_small;

// Functions
void data_init 	  		(DataSmall *data, int size, int fingerPrintSize);
void data_addFingerprint(DataSmall *data, char *fingerPrintS, char *id);
void data_sort 			(DataSmall *data);
double 	data_tanimotoSimilarity (fingerPrintIntType *queryFP, fingerPrintIntType *dataFP, int fpLen);
int 	data_getIntersection 	(fingerPrintIntType *queryFP, fingerPrintIntType *dataFP, int fpLen);
int 	data_getUnion 			(fingerPrintIntType *queryFP, fingerPrintIntType *dataFP, int fpLen);
int  	data_getIndexOf1Bits 	(int *FeatureIndex, fingerPrintIntType *fp, int len);
int 	data_getNum1Bits 		(fingerPrintIntType *fp, int fpLen);

//used only for testing
void data_printHex 		(DataSmall *data);
int  data_printFeatures (fingerPrintIntType *fp, int fpLen);
//inline functions
fingerPrintIntType* data_getFingerPrint (DataSmall *data, u_long index);

#endif /* DATA_H_ */
