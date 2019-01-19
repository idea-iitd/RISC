/*
 * data.c
 *
 *  Created on: 25-Nov-2018
 *      Author: jithin
 */
#include "data.h"

/**
 * Allocate memory
 * size : maximum number of molucules,
 * fingerPrintSize : the number of number of elements of type fingerPrintIntType in each fingerprint
 */
void data_init (DataSmall *data, int size, int fingerPrintSize) {
	data->maxSize		= size;
	data->size			= 0;
	data->_fpArray 		= (fingerPrintIntType *) helper_malloc (size*fingerPrintSize* sizeof (fingerPrintIntType),
			__FILE__, __LINE__);
	data->ids			= (char **) helper_malloc (size * sizeof(char *),
			__FILE__, __LINE__);
	data->fingerPrintLen= fingerPrintSize;

	int m = fingerPrintSize * NUMBITS_fingerPrintIntType;
	data->_sizeIndexStart	= (int *) helper_malloc((size_t)(m + 1) * sizeof(int),
			__FILE__, __LINE__);
	data->_sizeIndexEnd		= (int *) helper_malloc((size_t)(m + 1) * sizeof(int),
			__FILE__, __LINE__);

	for (int i=0; i<=m; i++) {
		data->_sizeIndexStart[i]=-1;
	}

	data->_isSorted = 0;
}


int data_getNum1Bits (fingerPrintIntType *fp, int fpLen) {
	int n=0;
	for (int i=0; i< fpLen; i++) {
		n += __builtin_popcountl (fp[i]);
		//printf ("%lu,", fp[i]);
	}
	//printf (",%d;\n", n);
	return n;
}

void _swapFingerprints (DataSmall *data, int a, int b) {
	fingerPrintIntType *fp_a = data_getFingerPrint(data, a);
	fingerPrintIntType *fp_b = data_getFingerPrint(data, b);

	for (int i=0; i <data->fingerPrintLen; i++) {
		fingerPrintIntType temp = fp_a[i];
		fp_a[i] = fp_b[i];
		fp_b[i] = temp;
	}
}

/**
 * sort the data based on the len
 */
void data_sort (DataSmall *data) {
	//maintain the data sorted by len
	int *lens= (int *)helper_malloc((data->size)*sizeof(int),__FILE__, __LINE__);

	//populate the sizes
	int sizeData 	= data->size;
	int fpLen 		= data->fingerPrintLen;

	data->maxNumFeatures = 0;
	data->minNumFeatures = data->maxNumFeatures;

	for (int i=0; i<sizeData; i++) {
		fingerPrintIntType *fp = data_getFingerPrint(data, i);
		int numBits = data_getNum1Bits(fp, fpLen);
		lens[i] = numBits;

		if (numBits > data->maxNumFeatures)
			data->maxNumFeatures = numBits;
		if (numBits < data->minNumFeatures)
			data->minNumFeatures = numBits;
	}

	//data_printHex(data);

	//sort
	int maxLen =  data->fingerPrintLen * NUMBITS_fingerPrintIntType;
	int currentIndex = 0;
	for (int l =1; l<maxLen; l++) {
		for (int k = currentIndex; k<sizeData; k++) {
			if (lens[k]==l) {
				//swap lengths
				lens[k] = lens[currentIndex];
				lens[currentIndex] = l;

				//swap the fingerprints
				_swapFingerprints (data, k, currentIndex);

				currentIndex++;
			}
		}
	}

	int maxSize=0;

	for (int i=0; i<sizeData; i++) {
		int l = lens[i];

		data->_sizeIndexStart[l]=i;

		while ((l == lens[i]) && (i<sizeData))
			i++;
		i--;
		data->_sizeIndexEnd[l]=i;

		int size = data->_sizeIndexEnd[l] - data->_sizeIndexStart[l] + 1;
		if (size > maxSize)
			maxSize = size;
	}
	data->_isSorted = 1;
	data->maxNumOfMoleculesOfAnySize = maxSize;


	//check
	/*
	{
		data_printHex(data);

		for (int i=0; i<sizeData; i++) {
			fingerPrintIntType *fp = data_getFingerPrint(data, i);
			int numBits = data_getNum1Bits(fp, fpLen);

			if (numBits != lens[i])
				helper_err(__FILE__, __LINE__, "length not matching");
		}
		printf ("\n len \n");
		for (int i=0; i<sizeData; i++) {
			printf ("%d,",lens[i]);
		}
		printf ("\n len \n");
		for (int i=0; i<=maxLen; i++) {
			int s = data->_sizeIndexStart[i];
			if (s >=0) {
				int e = data->_sizeIndexEnd[i];
				printf ("%d:(%d,%d)\n",i,s,e);
			}
		}
	}
	 */
	free (lens);
}

/**
 * convert the string of hexadecimal to a array of fingerprint
 */
void _convertFingerPrintHexStringToInt(char *fingerPrintS, int sizeI, fingerPrintIntType *fp) {
	int sizeX = sizeof(fingerPrintIntType)*2; // the number of hex characters to be used in one int

	int i=0;

	for (char *s=fingerPrintS; ; s+=sizeX) {
		char temp;
		int flag=0;

		if (strlen(s)>sizeX) {
			temp = s[sizeX];
			s[sizeX] = '\0';
			flag = 1;
		}
		fp[i] = strtoul(s, NULL, 16);

		if (flag)
			s[sizeX]=temp;
		else
			break;

		i++;
	}
}


inline fingerPrintIntType* data_getFingerPrint (DataSmall *data, u_long i) {
	i = i * data->fingerPrintLen;
	return &(data->_fpArray[i]);
}

/**
 * add the fingerprint to the dataset
 */
void data_addFingerprint(DataSmall *data, char *fingerPrintS, char *id) {
	int index = data->size;
	if (index >= data->maxSize) {
		char m[100];
		sprintf(m, "Not enough memory allocated, maxSize = %d", data->maxSize);
		helper_err(__FILE__, __LINE__, m);
	}


	fingerPrintIntType *fp = data_getFingerPrint (data, index);
	_convertFingerPrintHexStringToInt(fingerPrintS, data->fingerPrintLen, fp);

	data->ids[index] = strdup(id);

	data->size++;
}

/**
 * Find the similarity
 */
double data_tanimotoSimilarity (fingerPrintIntType *queryFP, fingerPrintIntType *dataFP, int fpLen) {
	int n=0,d=0;

	for (int i=0; i< fpLen; i++) {
		n += __builtin_popcountl (queryFP[i] & dataFP[i]);
		d += __builtin_popcountl (queryFP[i] | dataFP[i]);
	}

	double sim;
	if (d==0)
		sim = 0;
	else
		sim = (double)n/d;

	return sim;
}

/**
 * Find the intersection
 */
int data_getIntersection (fingerPrintIntType *queryFP, fingerPrintIntType *dataFP, int fpLen) {
	int n=0;

	for (int i=0; i< fpLen; i++) {
		n += __builtin_popcountl (queryFP[i] & dataFP[i]);
	}

	return n;
}

/**
 * Find the intersection
 */
int data_getUnion (fingerPrintIntType *queryFP, fingerPrintIntType *dataFP, int fpLen) {
	int n=0;

	for (int i=0; i< fpLen; i++) {
		n += __builtin_popcountl (queryFP[i] | dataFP[i]);
	}

	return n;
}

/**
 * populate FeatureIndex with the features having non zero value
 * and return the number of such features
 */
int  data_getIndexOf1Bits (int *FeatureIndex, fingerPrintIntType *fp, int len) {
	int bitLen = NUMBITS_fingerPrintIntType;

	fingerPrintIntType  mask = 1;
	mask = mask  << (bitLen -1);
	int n=0;
	for (int i=0; i<len; i++) {
		int j = bitLen;
		for (fingerPrintIntType num = fp[i]; num; num = num << 1) {
			j--;
			if (num & mask) {
				FeatureIndex[n++] = bitLen * (len-i-1) + j;
			}
		}
	}

	FeatureIndex[n]=-1;
	return n;
}

int data_printFeatures (fingerPrintIntType *fp, int len) {
	int bitLen = NUMBITS_fingerPrintIntType;

	fingerPrintIntType  mask = 1;
	mask = mask  << (bitLen -1);
	int n=0;
	for (int i=0; i<len; i++) {
		int j = bitLen;
		for (fingerPrintIntType num = fp[i]; num; num = num << 1) {
			j--;
			if (num & mask) {
				int f = bitLen * (len-i-1) + j;
				printf ("%d,",f);
				n++;
			}
		}
	}
	printf ("\n");
	return n;
}

/**
 * print the number in hex
 */
void data_printHex (DataSmall *data) {
	int bitLen = NUMBITS_fingerPrintIntType;

	fingerPrintIntType  mask = 1;
	mask = mask  << (bitLen -1);
	for (int i=0; i<data->size; i++) {
		printf ("%d\t:",i);
		fingerPrintIntType *fp = data_getFingerPrint(data, i);

		for (int k=0; k<data->fingerPrintLen; k++) {
			fingerPrintIntType val = fp[k];

			for (int j=0; j<bitLen; j+=4) {

				char text[5]="XXXX";

				for (int u=0; u<4; u++) {
					if (val & mask)
						text[u]='1';
					else
						text[u]='0';

					val = val<<1;
				}

				char c='0';
				if      (strcmp(text, "0000")==0)
					c='0';
				else if (strcmp(text, "0001")==0)
					c='1';
				else if (strcmp(text, "0010")==0)
					c='2';
				else if (strcmp(text, "0011")==0)
					c='3';
				else if (strcmp(text, "0100")==0)
					c='4';
				else if (strcmp(text, "0101")==0)
					c='5';
				else if (strcmp(text, "0110")==0)
					c='6';
				else if (strcmp(text, "0111")==0)
					c='7';
				else if (strcmp(text, "1000")==0)
					c='8';
				else if (strcmp(text, "1001")==0)
					c='9';
				else if (strcmp(text, "1010")==0)
					c='a';
				else if (strcmp(text, "1011")==0)
					c='b';
				else if (strcmp(text, "1100")==0)
					c='c';
				else if (strcmp(text, "1101")==0)
					c='d';
				else if (strcmp(text, "1110")==0)
					c='e';
				else if (strcmp(text, "1111")==0)
					c='f';

				printf ("%c",c);
			}
		}
		printf (",\n");
	}
	printf ("\n");
}
