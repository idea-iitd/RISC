/*
 * nonBinaryFingerPrint.h
 *
 *  Created on: 30-Dec-2018
 *      Author: jithin
 */

#ifndef NONBINARYFINGERPRINT_H_
#define NONBINARYFINGERPRINT_H_

#include "global.h"
#include "options.h"
#include "arrayListInt.h"
#include "minHeapPlain.h"
#include "dataNonBinary.h"
#include "experiments.h"

#include <dirent.h>
#include <stdlib.h>

void nonBinaryFingerPrint_run_public ();
void nonBinaryFingerPrint_run (enum options_datasetTypes d);
workerFunctions_NB_type linear_NBF_getWorkerFunction ();
workerFunctions_NB_type invertedIndex_NBF_getWorkerFunction ();

#endif /* NONBINARYFINGERPRINT_H_ */
