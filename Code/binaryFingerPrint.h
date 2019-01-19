/*
 * binaryFingerPrint.h
 *
 *  Created on: 16-Dec-2018
 *      Author: jithin
 */

#ifndef BINARYFINGERPRINT_H_
#define BINARYFINGERPRINT_H_

#include "global.h"
#include "options.h"
#include "experiments.h"
#include "dataBinary.h"
#include "minHeapPlain.h"
#include "invertedIndex_BF.h"
#include "divideSkip_BF.h"
#include "AOR_BF.h"

void binaryFingerPrint_run (enum options_datasetTypes d);
void binaryFingerPrint_run_public ();

void runLinear_BF (DataBinary *data, DataBinary *dataQueries, char *rnamePartial);

workerFunctions_type linear_BF_getWorkerFunction ();

#endif /* BINARYFINGERPRINT_H_ */
