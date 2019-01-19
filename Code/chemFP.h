/*
 * chenFP.h
 *
 *  Created on: 23-Nov-2018
 *      Author: jithin
 */

#ifndef CHEMFP_H_
#define CHEMFP_H_

#include "global.h"
#include "data.h"
#include "helper.h"
#include "experiments.h"
#include "options.h"
#include "invertedIndexSmall.h"
#include "divideSkip.h"
#include "AOR.h"

void chemFP_run_public ();
void chemFP_run(char *fileID, int numBits);
workerFunctions_type_small linear_getWorkerFunction ();

#endif /* CHEMFP_H_ */
