/*
 * options.h
 *
 *  Created on: 20-Nov-2018
 *      Author: jithin
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

enum options_datasetTypes_Public {
	_START_data_Public=-1,
	_EXIT_data_Public,
	chemFP_Public,
	BINARY_Public,
	NON_BINARY_Public,
	_END_data_Public
};

enum options_datasetTypes_Public options_getDatasetTypes_Public ();
char* options_getDataSetName_Public (enum options_datasetTypes_Public d);


enum options_datasetTypes {
	_START_data=-1,
	_EXIT_data,
	chemFP_166,
	chemFP_881,
	chemFP_1024,
	chemFP_2048,
	chemFP_DUD,
	chemFP_ZINC_400K,
	DUD,
	pubChem,
	ZINC,
	ZINC_RDKit,
	pubChem_RDKit,
	pubChem3DLOOB,
	pubChem3DE3fp_14,
	_END_data
};

enum options_datasetTypes options_getDatasetTypes ();
char* options_getDataSetName (enum options_datasetTypes d);

/**
 *
 * Experiments
 *
 */
enum options_ExperimentsTypes {
	_START_exp=-1,
	_EXIT_exp,
	Range,
	TopK,
	_END_exp
};

enum options_ExperimentsTypes options_getExperimentsTypes ();
char* options_getExperminetnName(enum options_ExperimentsTypes d);

/**
 * Indexing scheme
 */
enum options_IndexTypes {
	_START_index=-1,
	_EXIT_index,
	Linear,
	RISC,
	DivideSkip,
	AOR,
	_END_index
};
enum options_IndexTypes options_getIndexTypes ();
char* options_getIndexingSchemeName(enum options_IndexTypes d);

#endif /* OPTIONS_H_ */
