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

	bool useNibble = false;
	if(sArgs.is_nibble_flag==1) useNibble = true;


	CSeekCentral *func = new CSeekCentral();
	if(!func->Initialize(sArgs.input_arg, sArgs.func_quant_arg,
		sArgs.func_dset_arg, sArgs.func_dset_arg, sArgs.query_arg,
		sArgs.func_prep_arg, sArgs.func_db_arg, sArgs.func_prep_arg,
		useNibble, sArgs.func_n_arg, sArgs.buffer_arg,
		true, true, true, true)){
		return -1;
	}

	func->EqualWeightSearch();
	const vector< vector<AResultFloat> > &vfunc = func->GetAllResult();
	const vector<CSeekQuery> &vq = func->GetAllQuery();

	vector<vector<string> > newQuery;
	newQuery.resize(vfunc.size());
	ushort i,j;
	ushort TOP = 10;
	for(i=0; i<vfunc.size(); i++){
		newQuery[i] = vector<string>();
		const vector<ushort> &queryGenes = vq[i].GetQuery();
		for(j=0; j<queryGenes.size(); j++){
			newQuery[i].push_back(func->GetGene(queryGenes[j]));
		}
		for(j=0; j<TOP; j++){
			newQuery[i].push_back(func->GetGene(vfunc[i][j].i));
		}
	}

	func->Destruct();
	delete func;

	CSeekTools::Write2DArrayText("/tmp/expanded_query.txt", newQuery);

	CSeekCentral *csk = new CSeekCentral();

	if(!csk->Initialize(sArgs.input_arg, sArgs.quant_arg, sArgs.dset_arg,
		sArgs.search_dset_arg, "/tmp/expanded_query.txt",
		sArgs.dir_platform_arg,
		sArgs.dir_in_arg, sArgs.dir_prep_in_arg, useNibble, sArgs.num_db_arg,
		sArgs.buffer_arg,
		!!sArgs.norm_subavg_flag, !!sArgs.norm_platsubavg_flag,
		!!sArgs.norm_platstdev_flag, false))
			return -1;

	/* Random Number Generator Initializations */
	const gsl_rng_type *T;
	gsl_rng *rnd;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rnd = gsl_rng_alloc(T);
	float RATE = 0.95;
	ushort FOLD = 5;
	enum PartitionMode PART_M = CUSTOM_PARTITION;

	csk->CVSearch(rnd, PART_M, FOLD, RATE);
	const vector<vector<float> > &csk_weight = csk->GetAllWeight();
	vector<vector<float> > csk_weight_copy;
	csk_weight_copy.resize(csk_weight.size());
	for(i=0; i<csk_weight.size(); i++){
		csk_weight_copy[i] = vector<float>();
		for(j=0; j<csk_weight[i].size(); j++)
			csk_weight_copy[i].push_back(csk_weight[i][j]);
	}

	csk->Destruct();
	delete csk;


	CSeekCentral *csfinal = new CSeekCentral();
	if(!csfinal->Initialize(sArgs.input_arg, sArgs.quant_arg, sArgs.dset_arg,
		sArgs.search_dset_arg, sArgs.query_arg, sArgs.dir_platform_arg,
		sArgs.dir_in_arg, sArgs.dir_prep_in_arg, useNibble, sArgs.num_db_arg,
		sArgs.buffer_arg,
		!!sArgs.norm_subavg_flag, !!sArgs.norm_platsubavg_flag,
		!!sArgs.norm_platstdev_flag, false))
		return -1;

	csfinal->WeightSearch(csk_weight_copy);
	csfinal->Destruct();
	delete csfinal;

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }
