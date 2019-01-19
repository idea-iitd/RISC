/*
 * dataNonBinary.c
 *
 *  Created on: 30-Dec-2018
 *      Author: jithin
 */

#include "dataNonBinary.h"

MolNonBinary* _new_MolNonBinary () {
	MolNonBinary *mol = (MolNonBinary *) helper_calloc(1, sizeof(MolNonBinary), __FILE__, __LINE__);
	mol->features 	  = new_arrayListtype_flexible(10, 100000);
	mol->values 	  = new_arrayListtype_flexible(10, 100000);

	return mol;
}

u_long _get_newFeature_NB (u_long feature, arrayList_Iterator_type *iterator) {
	int found;
	u_long newFeature = arrayList_IteratorFind_Binary(iterator, feature, &found);
	if (!found)
		helper_err(__FILE__, __LINE__, "feature missing");

	return newFeature;
}

void _MolNonBinary_add (MolNonBinary *mol, u_long feature, u_long val) {
	arrayList_flexible_put(&mol->features, feature);
	arrayList_flexible_put(&mol->values, val);
	mol->cummulativeFeatureVal += val;
}

void _sort_MolNonBinary (MolNonBinary *mol) {
	int n = mol->features->size;

	u_long *f = mol->features->data;
	u_long *v = mol->values->data;

	for (int i=0; i<n;i++) {
		for (int j=i+1;j<n;j++) {
			if (f[i] > f[j]) {
				//swap
				u_long temp = f[i];
				f[i] = f[j];
				f[j] = temp;

				temp = v[i];
				v[i] = v[j];
				v[j] = temp;
			}
		}
	}
}

void dataNonBinary_addMolecule  (DataNonBinary *data, char *buffer, arrayList_Iterator_type *iterator,
		arrayListtype *featureIDs) {
	int numFeatures = 0;

	if (buffer[0] != '#') {
		//its not a valid fingerprint
		return;
	}

	arrayList_resetIterator(iterator, featureIDs);
	//load features to the min heap,
	char *str 	= strtok (buffer," ");
	int l = strlen(str);
	char *molName = (char *) helper_calloc(l+2, sizeof(char), __FILE__, __LINE__);

	strcpy(molName, str);

	str = strtok ((char *)NULL, " "); //first value is molID
	int val=-1;

	MolNonBinary *mol = _new_MolNonBinary ();

	do {
		u_long feature;
		val = -1;
		sscanf(str, "%lu:%d",&feature,&val);

		if (val > 0) {
			u_long newFeature = _get_newFeature_NB (feature, iterator);

			_MolNonBinary_add (mol, newFeature, val);
			numFeatures++;
		}

		str = strtok ((char *)NULL, " ");
	} while (str);

	_sort_MolNonBinary (mol);

	if (data->maxNumFeatures < numFeatures)
		data->maxNumFeatures = numFeatures;

	if (data->size >= data->_maxSize) {
		//memory is not sufficient
		helper_err(__FILE__, __LINE__, "memory is not sufficient");
	}

	data->ids     [data->size] 	= molName;
	data->_fpArray[data->size++]= mol;

	{
		u_long len = mol->features->size;
		if (data->minLen > len)
			data->minLen = len;
		if (data->maxLen < len)
			data->maxLen = len;
	}
}


void dataNonBinary_foldAndWrite2File  (char *buffer, FILE *fp, u_long fold, int id) {
	if (buffer[0] != '#') {
		//its not a valid fingerprint
		return;
	}

	//load features to the min heap,
	char *str 	= strtok (buffer," ");
	str = strtok ((char *)NULL, " "); //first value is molID
	int val=-1;


	arrayListtype_flexible *featuresList= new_arrayListtype_flexible(1000, 100000);
	arrayListtype_flexible *valuesList 	= new_arrayListtype_flexible(1000, 100000);

	do {
		u_long feature;
		val = -1;
		sscanf(str, "%lu:%d",&feature,&val);

		if (val > 0) {
			//folding
			feature = feature % fold;

			arrayList_flexible_put(&featuresList, feature);
			arrayList_flexible_put(&valuesList, val);
		}

		str = strtok ((char *)NULL, " ");
	} while (str);

	int n = featuresList->size;

	if (n > 0) {
		u_long *features= featuresList->data;
		u_long *values 	= valuesList->data;

		//sorting the features
		for (int i=0;i<n; i++){
			for (int j=i+1;j<n;j++){
				if (features[i] > features[j]) {
					u_long temp;

					temp 		= features[i];
					features[i] = features[j];
					features[j] = temp;

					temp 		= values[i];
					values[i] 	= values[j];
					values[j] 	= temp;
				}
			}
		}
		//write to file considering the duplicate values
		fprintf (fp, "#%d",id);
		for (int i=0;i<n; i++){
			u_long feature 	= features[i];
			u_long value 	= values[i];
			for (int j=i+1;j<n;j++){
				if (features[j]==feature) {
					//take the max value
					if (values[j] > value)
						value = values[j];
					i++;
				} else {
					break;
				}
			}
			fprintf (fp, " %lu:%lu",feature, value);
		}
		fprintf (fp, "\n");
	}
}


MolNonBinary* dataNonBinary_getFingerPrint (DataNonBinary *data, u_long index) {
	return data->_fpArray[index];
}

DataNonBinary* new_dataNonBinary (u_long maxSize, u_long numFeatures) {
	DataNonBinary *data = (DataNonBinary *)helper_calloc(1,       sizeof(DataNonBinary), __FILE__, __LINE__);
	data->_fpArray      = (MolNonBinary **)helper_calloc(maxSize, sizeof(MolNonBinary*), __FILE__, __LINE__);
	data->ids           = (char         **)helper_calloc(maxSize, sizeof(char*),         __FILE__, __LINE__);

	data->_maxSize  	= maxSize;
	data->numFeatures 	= numFeatures;
	data->minLen		= INT_MAX;
	return data;
}

double 	dataNonBinary_tanimotoSimilarity (MolNonBinary *queryFP, MolNonBinary *dataFP) {
	int index_1=0, l_1=queryFP->features->size;
	int index_2=0, l_2=dataFP->features->size;

	u_long *features_1 = queryFP->features->data, *values_1 = queryFP->values->data;
	u_long *features_2 = dataFP->features->data,  *values_2 = dataFP->values->data;

	double n=0, d=0;

	while ((index_1 < l_1) && (index_2 < l_2)) {
		u_long f1 = features_1[index_1];
		u_long f2 = features_2[index_2];

		if (f1==f2) {
			u_long v_1 = values_1[index_1];
			u_long v_2 = values_2[index_2];

			if (v_1 < v_2) {
				n += v_1;
				d += v_2;
			} else {
				n += v_2;
				d += v_1;
			}

			index_1++;
			index_2++;
		} else if (f1<f2) {
			d += values_1[index_1];
			index_1++;
		} else if (f1>f2) {
			d += values_2[index_2];
			index_2++;
		}
	}

	while (index_1 < l_1) {
		d += values_1[index_1];
		index_1++;
	}

	while (index_2 < l_2) {
		d += values_2[index_2];
		index_2++;
	}

	return n/d;
}

void _dataNonBinary_swapFingerprints (DataNonBinary *data, u_long a, u_long b) {
	MolNonBinary *fp_a = dataNonBinary_getFingerPrint(data, a);
	MolNonBinary *fp_b = dataNonBinary_getFingerPrint(data, b);

	data->_fpArray[a] = fp_b;
	data->_fpArray[b] = fp_a;

	char *tmp = data->ids[a];
	data->ids[a]= data->ids[b];
	data->ids[b] = tmp;
}

void dataNonBinary_print(MolNonBinary *fp) {
	u_long len = fp->features->size;
	for (int i=0; i<len; i++) {
		printf ("%lu:%lu, ",fp->features->data[i],fp->values->data[i]);
	}
	printf ("\n");
}

/**
 * sort the data based on the len and with in len sort basedon the cummulative val
 */
void dataNonBinary_sort (DataNonBinary *data) {
	//populate the sizes
	int sData 	= data->size;
	u_long maxLen = data->maxLen;
	u_long minLen = data->minLen;

	do {
		//allocate memory
		data->_sizeIndexStart= (int *) helper_malloc((size_t)(maxLen+1) * sizeof(int),  __FILE__, __LINE__);
		data->_sizeIndexEnd	 = (int *) helper_malloc((size_t)(maxLen+1) * sizeof(int),  __FILE__, __LINE__);

		for (int i=0; i<=maxLen; i++) {
			data->_sizeIndexStart[i]=-1;
		}

		//maintain the data sorted by len
		int *lens = (int *) helper_malloc((size_t)(sData) * sizeof(int),  __FILE__, __LINE__);

		for (int i=0; i<sData; i++) {
			MolNonBinary *fp 	= dataNonBinary_getFingerPrint(data, i);

			lens[i] = fp->features->size;
		}

		//sort
		int currentIndex = 0;
		for (int l =minLen; l<=maxLen; l++) {
			for (int k = currentIndex; k<sData; k++) {
				if (lens[k]==l) {
					//swap lengths
					lens[k] = lens[currentIndex];
					lens[currentIndex] = l;

					//swap the fingerprints
					_dataNonBinary_swapFingerprints (data, k, currentIndex);

					currentIndex++;
				}
			}
		}

		int maxSize=0;
		for (int i=0; i<sData; i++) {
			int l = lens[i];

			data->_sizeIndexStart[l]=i;

			while ((l == lens[i]) && (i<sData))
				i++;
			i--;
			data->_sizeIndexEnd[l]=i;

			int size = data->_sizeIndexEnd[l] - data->_sizeIndexStart[l] + 1;
			if (size > maxSize)
				maxSize = size;
		}
		data->maxNumOfMoleculesOfAnySize = maxSize;

		free (lens);
	} while (0);

	//sort based on cumulative value
	do {
		int *cv = (int *) helper_malloc((size_t)(sData) * sizeof(int),  __FILE__, __LINE__);

		int minCV = INT_MAX;
		int maxCV = 0;

		for (int i=0; i<sData; i++) {
			MolNonBinary *fp 	= dataNonBinary_getFingerPrint(data, i);

			int v = fp->cummulativeFeatureVal;
			cv[i] = v;

			if (v<minCV)
				minCV = v;
			if (v>maxCV)
				maxCV = v;
		}

		//sort
		for (int len=minLen; len<maxLen; len++) {
			int startINdex = data->_sizeIndexStart[len];
			if (startINdex >= 0) {
				int endINdex = data->_sizeIndexEnd[len];

				//sort in this len range
				{
					int currentIndex = startINdex;
					for (int v=minCV; v<=maxCV; v++) {
						for (int k = currentIndex; k<=endINdex; k++) {
							if (cv[k]==v) {
								//swap cv
								cv[k] = cv[currentIndex];
								cv[currentIndex] = v;

								//swap the fingerprints
								_dataNonBinary_swapFingerprints (data, k, currentIndex);

								currentIndex++;
							}
						}
					}
				}
			}
		}

		//allocate memory
		data->_cvMin = (int *) helper_malloc((size_t)(maxLen+1) * sizeof(int),  __FILE__, __LINE__);
		data->_cvMax = (int *) helper_malloc((size_t)(maxLen+1) * sizeof(int),  __FILE__, __LINE__);

		for (int len=minLen; len<maxLen; len++) {
			int startINdex = data->_sizeIndexStart[len];
			if (startINdex >= 0) {
				int endINdex = data->_sizeIndexEnd[len];

				data->_cvMin[len]=cv[startINdex];
				data->_cvMax[len]=cv[endINdex];

				if (data->maxCV < cv[endINdex])
					data->maxCV = cv[endINdex];
			}
		}

		free (cv);

	} while (0);

	data->_isSorted = 1;
}
