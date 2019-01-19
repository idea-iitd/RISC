/*
 * global.h
 *
 *  Created on: 23-Nov-2018
 *      Author: jithin
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define G_DATA_DIR   "../../Data/"
#define G_RESULT_DIR "../../Data/Result/"
#define G_PRECISION 1000
#define _G_AOR 40
#define G_LOCATION //do{fprintf(stderr,"Entered File:%s, Function:%s, line:%d\n", __FILE__, __FUNCTION__, __LINE__);}while(0)

typedef unsigned long u_long;

//#define G_GET_PRUNING_RATE

#endif /* GLOBAL_H_ */
