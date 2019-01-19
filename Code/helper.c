/*
 * helper.c
 *
 *  Created on: 23-Nov-2018
 *      Author: jithin
 */
#include "helper.h"

void helper_printProgress (int n) {
	if (n % 1000 == 0) {
		int m = n/10+1;

		while (m--) {
			fprintf (stderr,"\b");
		}
		fprintf (stderr,"%d",n);
	}
}

/**
 * print error filename, line and message
 */
void helper_err (char *fName, int line, char *m) {
	fprintf(stderr, "%s (line %d) %s \n",fName,line, m);
	exit(-1);
}

/**
 * print error filename, line and message
 */
void helper_errMem (char *fName, int line) {
	fprintf(stderr, "%s (line %d) : Not enough memory allocated\n",fName,line);
	exit(-1);
}

/**
 * print error filename, line and message
 */
void helper_perror (char *fName, int line, char *m) {
	char message[500];
	sprintf(message, "%s (line %d) %s ",fName,line, m);
	perror(message);
	exit(-1);
}

void helper_scanfErr (int expectedValue, int ourVal, char *fName, int line) {
	if (ourVal != expectedValue) {
		helper_perror(fName, line, "scanf failed");
	}
}

/**
 * Check if file pointer is NULL
 */
void helper_errFP (FILE *fp, char *fName, char *fName_code, int line) {
	if (fp == NULL) {
		helper_perror(fName_code, line, fName);
	}
}

void *helper_malloc (size_t __size, char *fname, int lineNum) {
	void *ret = malloc(__size);

	if (ret == NULL) {
		char message[500];
		sprintf(message, "%s (line %d) malloc asking for %zu bytes", fname, lineNum, __size);
		perror(message);
		exit (-1);
	}

	return ret;
}

void *helper_calloc (size_t __nmemb, size_t __size, char *fname, int lineNum){
	void *ret = calloc(__nmemb, __size);

	if (ret == NULL) {
		char message[500];
		sprintf(message, "%s (line %d) calloc asking for %zu elments of size %zu bytes",
				fname, lineNum, __nmemb, __size);
		perror(message);
		exit (-1);
	}

	return ret;
}
