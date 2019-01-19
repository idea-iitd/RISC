/*
 * experiments.h
 *
 *  Created on: 25-Nov-2018
 *      Author: jithin
 */

#ifndef EXPERIMENTS_H_
#define EXPERIMENTS_H_

#include "global.h"
#include "data.h"
#include "minHeap.h"
#include "minHeapInt.h"
#include "arrayListInt.h"
#include "binaryFingerPrint.h"
#include "nonBinaryFingerPrint.h"
#include "dataBinary.h"
#include "dataNonBinary.h"
#include "invertedIndexSmall.h"
#include "AOR.h"
#include "divideSkip.h"
#include "chemFP.h"

extern int 	kStart_G;
extern int 	kEnd_G;
extern int 	k_G[];
extern int 	rSize_G;
extern double r_G[];

extern int num_exp_G;

extern int CHECK_SOLUTION_G;

extern int releasePublic;

void experiments_work 			(DataBinary *data, 		DataBinary *dataQueries, 	char *rnamePartial);
void experiments_work_small 	(DataSmall  *data, 		DataSmall  *dataQueries, 	char *rnamePartial);
void experiments_nonBinary_work (DataNonBinary *data, 	DataNonBinary *dataQueries, char *rnamePartial);
void notImplemented (void *source, char *fname);

#endif /* EXPERIMENTS_H_ */
