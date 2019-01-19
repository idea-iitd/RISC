/*
 * helper.h
 *
 *  Created on: 20-Nov-2018
 *      Author: jithin
 */

#ifndef HELPER_H_
#define HELPER_H_

#include "global.h"

void helper_perror 	(char *fName, int line, char *m);
void helper_err 	(char *fName, int line, char *m);
void helper_errMem 	(char *fName, int line);
void helper_scanfErr(int expectedValue, int ourVal, char *fName, int line);
void helper_errFP   (FILE *fp, char *fName, char *fName_code, int line);

void *helper_malloc (size_t __size, char *fname, int lineNum);
void *helper_calloc (size_t __nmemb, size_t __size, char *fname, int lineNum);

void helper_printProgress (int n);

#endif /* HELPER_H_ */
