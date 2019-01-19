/*
 * options.c
 *
 *  Created on: 20-Nov-2018
 *      Author: jithin
 */
#include "global.h"
#include "options.h"
#include "helper.h"
/**
 * display the user messages
 */
char* _getMessage_data(enum options_datasetTypes d) {
	switch (d) {
	case chemFP_1024 :
		return ("chemFP_1024");
	case chemFP_166 :
		return ("chemFP_166");
	case chemFP_881 :
		return ("chemFP_811");
	case chemFP_2048 :
		return ("chemFP_2048");
	case chemFP_ZINC_400K:
		return ("chemFP_ZINC_400K");
	case chemFP_DUD:
		return ("chemFP_DUD");
	case _EXIT_data :
		return ("EXIT");
	case DUD:
		return ("DUD");
	case pubChem:
		return ("pubChem");
	case ZINC:
		return ("ZINC");
	case pubChem_RDKit:
		return ("pubChem_RDKit");
	case ZINC_RDKit:
		return ("ZINC_RDKit");
	case pubChem3DE3fp_14:
		return ("pubChem3D_E3fp_14");
	case pubChem3DLOOB:
		return ("pubChem3D_LOOB");
	case _START_data:
	case _END_data:
		return ("");
	}
	return ("");
}

/**
 * display the user messages
 */
char* _getMessage_data_Public(enum options_datasetTypes_Public d) {
	switch (d) {
	case _EXIT_data_Public :
		return ("EXIT");
	case _START_data_Public:
	case _END_data_Public:
		return ("");
	case chemFP_Public:
		return ("FPS");
	case BINARY_Public:
		return ("Binary Data");
	case NON_BINARY_Public:
		return ("Non Binary");
	}
	return ("");
}

char* options_getDataSetName (enum options_datasetTypes d) {
	return _getMessage_data(d);
}

/**
 * display the user messages
 */
char* _getMessage_exp(enum options_ExperimentsTypes d) {
	switch (d) {
	case TopK :
		return ("TopK");
	case Range :
		return ("Threshold");
	case _EXIT_exp :
		return ("EXIT");
	case _START_exp:
	case _END_exp:
		return ("");
	}
	return ("");
}

char* options_getExperminetnName(enum options_ExperimentsTypes d) {
	return _getMessage_exp(d);
}

/**
 * display the user messages
 */
char* _getMessage_index(enum options_IndexTypes d) {
	switch (d) {
	case Linear:
		return ("LINEAR");
	case RISC:
		return ("RISC");
	case _EXIT_index :
		return ("EXIT");
	case AOR:
		return ("AOR");
	case DivideSkip:
		return ("DivideSkip");
	case _START_index:
	case _END_index:
		return ("");
	}
	return ("");
}

char* options_getIndexingSchemeName(enum options_IndexTypes d) {
	return _getMessage_index(d);
}

int _getOption() {
	printf("Your option : ");

	int o;
	int r = scanf("%d", &o);
	helper_scanfErr(1, r, __FILE__, __LINE__);

	return o;
}

/**
 * get the user input
 */
enum options_datasetTypes options_getDatasetTypes () {
	 printf ("Which Data-set do you want to be loaded? \n");
	for (int i=_START_data+1; i<_END_data; i++){
		printf ("%d : %s\n",i, _getMessage_data(i));
	}

	int o = _getOption();

	printf ("%s\n",_getMessage_data(o));

	return (enum options_datasetTypes)o;
}

/**
 * get the user input
 */
enum options_datasetTypes_Public options_getDatasetTypes_Public () {
	 printf ("Which Data-set do you want to be loaded? \n");
	for (int i=_START_data_Public+1; i<_END_data_Public; i++){
		printf ("%d : %s\n",i, _getMessage_data_Public(i));
	}

	int o = _getOption();

	printf ("%s\n",_getMessage_data_Public(o));

	return (enum options_datasetTypes_Public)o;
}


/**
 * get user input
 */
enum options_ExperimentsTypes options_getExperimentsTypes () {
	 printf ("Which experiment do you want to be run? \n");
	for (int i=_START_exp+1; i<_END_exp; i++){
		printf ("%d : %s\n",i, _getMessage_exp(i));
	}

	int o = _getOption();
	printf ("%s\n",_getMessage_exp(o));
	return (enum options_ExperimentsTypes)o;
}

/**
 * get user input
 */
enum options_IndexTypes options_getIndexTypes (){
	printf ("Which Indexing scheme do you want to be use? \n");
	for (int i=_START_index+1; i<_END_index; i++){
		printf ("%d : %s\n",i, _getMessage_index(i));
	}

	int o = _getOption();
	printf ("%s\n",_getMessage_index(o));
	return (enum options_ExperimentsTypes)o;
}
