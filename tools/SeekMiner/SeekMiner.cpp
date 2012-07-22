/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "cmdline.h"


int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuffer	= 1024;
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	istream*			pistm;
	vector<string>		vecstrGenes;
	char				acBuffer[ c_iBuffer ];
	ushort				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1;
	}

	if(!sArgs.input_arg || !sArgs.quant_arg ||
		!sArgs.dset_arg || !sArgs.search_dset_arg ||
		!sArgs.query_arg || !sArgs.dir_platform_arg ||
		!sArgs.dir_in_arg || !sArgs.dir_prep_in_arg){
		fprintf(stderr, "Arguments missing!\n");
		return 1;
	}

	/* Functional network */
	/*CDatabase FN(false);
	FN.Open(sArgs.func_db_arg, vecstrGenes, 1, sArgs.func_n_arg);
	vector<float> FNquant;
	CSeekTools::ReadQuantFile(sArgs.func_quant_arg, FNquant);
	vector<CSeekPlatform> FNvp;
	map<string, ushort> FNmapstriPlatform;
	vector<string> FNvecstrPlatforms;
	CSeekTools::ReadPlatforms(sArgs.func_prep_arg, FNvp, FNvecstrPlatforms,
		FNmapstriPlatform);
	vector<string> FNvecstrDatasets;
	FNvecstrDatasets.push_back("global_functional_network");
	ap<string, string> FNmapstrstrDP;
	FNmapstrstrDP["global_functional_network"] = FNvecstrPlatforms[0];
	vector<CSeekDataset*> FNvc;
	CSeekTools::LoadDatabase(FN, sArgs.func_prep_arg,
		FNvecstrDatasets, FNmapstrstrDP, FNmapstriPlatform, FNvp, FNvc);
	CSeekTools::ReadDatabaselets(FN, vecstrAllQuery, FNvc);*/
	//=====================================================================

	bool useNibble = false;
	if(sArgs.is_nibble_flag==1) useNibble = true;

	CSeekCentral csk;
	csk.Initialize(sArgs.input_arg, sArgs.quant_arg, sArgs.dset_arg,
		sArgs.search_dset_arg, sArgs.query_arg, sArgs.dir_platform_arg,
		sArgs.dir_in_arg, sArgs.dir_prep_in_arg, useNibble, sArgs.num_db_arg,
		!!sArgs.norm_subavg_flag, !!sArgs.norm_platsubavg_flag,
		!!sArgs.norm_platstdev_flag, false);

	/* Random Number Generator Initializations */
	const gsl_rng_type *T;
	gsl_rng *rnd;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rnd = gsl_rng_alloc(T);
	float RATE = 0.95;
	ushort FOLD = 5;
	enum PartitionMode PART_M = CUSTOM_PARTITION;

	csk.CVSearch(rnd, PART_M, FOLD, RATE);
	csk.Destruct();

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }
