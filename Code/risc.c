/*
 * risc.c
 *
 *  Created on: 18-Nov-2018
 *      Author: jithin
 */

#include "global.h"
#include "options.h"
#include "chemFP.h"
#include "binaryFingerPrint.h"
#include "nonBinaryFingerPrint.h"

int main () {
	 printf ("Welcome to RISC\n");

	 if (releasePublic) {
		 enum options_datasetTypes_Public d = options_getDatasetTypes_Public ();

		 switch (d) {
		 case _START_data_Public:
		 case _EXIT_data_Public:
		 case _END_data_Public:
			 break;
		 case chemFP_Public:
			 chemFP_run_public ();
			 break;
		 case BINARY_Public:
			 binaryFingerPrint_run_public();
			 break;
		 case NON_BINARY_Public:
			 nonBinaryFingerPrint_run_public();
			 break;
		 }
	 } else {
		 enum options_datasetTypes d = options_getDatasetTypes ();

		 switch (d) {
		 case chemFP_1024:
			 chemFP_run ("1024", 1024);
			 break;
		 case chemFP_166:
			 chemFP_run ("0166",166);
			 break;
		 case chemFP_881:
			 chemFP_run ("0881",881);
			 break;
		 case chemFP_2048:
			 chemFP_run ("2048", 2048);
			 break;
		 case chemFP_ZINC_400K:
			 chemFP_run ("zinc_400K", 34138);
			 break;
		 case chemFP_DUD:
			 chemFP_run ("DUD", 32198);
			 break;
		 case DUD:
		 case pubChem:
		 case pubChem_RDKit:
		 case ZINC:
		 case ZINC_RDKit:
			 binaryFingerPrint_run (d);
			 break;
		 case pubChem3DE3fp_14:
		 case pubChem3DLOOB:
			 nonBinaryFingerPrint_run (d);
			 break;
		 case _END_data:
		 case _EXIT_data:
		 case _START_data:
			 break;
		 }
	 }

	 printf ("Bye\n");
	 return 0;
 }

