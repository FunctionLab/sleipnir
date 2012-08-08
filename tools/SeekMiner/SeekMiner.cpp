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

	// Random Number Generator Initializations
	const gsl_rng_type *T;
	gsl_rng *rnd;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rnd = gsl_rng_alloc(T);
	float RATE = 0.95;
	ushort FOLD = 5;
	enum PartitionMode PART_M = CUSTOM_PARTITION;
	ushort i,j;
	ushort TOP = 1000;

/*
	CSeekCentral *func = new CSeekCentral();
	if(!func->Initialize(sArgs.input_arg, sArgs.func_quant_arg,
		sArgs.func_dset_arg, sArgs.func_dset_arg, sArgs.query_arg,
		sArgs.func_prep_arg, sArgs.func_db_arg, sArgs.func_prep_arg,
		useNibble, sArgs.func_n_arg, sArgs.buffer_arg,
		"fn", false, false, false, false,
		sArgs.score_cutoff_arg, sArgs.per_q_required_arg)){
		return -1;
	}

	func->EqualWeightSearch();
	const vector< vector<AResultFloat> > &vfunc = func->GetAllResult();
	const vector<CSeekQuery> &vq = func->GetAllQuery();

	vector<vector<string> > origQuery;
	vector<vector<string> > newQuery;
	newQuery.resize(vfunc.size());
	origQuery.resize(vfunc.size());

	for(i=0; i<vfunc.size(); i++){
		newQuery[i] = vector<string>();
		origQuery[i] = vector<string>();
		const vector<ushort> &queryGenes = vq[i].GetQuery();
		for(j=0; j<queryGenes.size(); j++){
			origQuery[i].push_back(func->GetGene(queryGenes[j]));
			newQuery[i].push_back(func->GetGene(queryGenes[j]));
		}
		for(j=0; j<200; j++)
			newQuery[i].push_back(func->GetGene(vfunc[i][j].i));
	}

	func->Destruct();
	delete func;
	
	CSeekTools::Write2DArrayText("/tmp/ex_query.txt", newQuery);
*/
/*	
	CSeekCentral *csk = new CSeekCentral();

	if(!csk->Initialize(sArgs.input_arg, sArgs.quant_arg, sArgs.dset_arg,
		sArgs.search_dset_arg, sArgs.query_arg,
		sArgs.dir_platform_arg,
		sArgs.dir_in_arg, sArgs.dir_prep_in_arg, useNibble, sArgs.num_db_arg,
		sArgs.buffer_arg, "normal",
		!!sArgs.norm_subavg_flag, !!sArgs.norm_platsubavg_flag,
		!!sArgs.norm_platstdev_flag, false,
		sArgs.score_cutoff_arg, sArgs.per_q_required_arg))
			return -1;


	//csk->CVCustomSearch(newQuery, rnd, PART_M, FOLD, RATE);
	csk->CVSearch(rnd, PART_M, FOLD, RATE);
	const vector<vector<float> > &csk_weight = csk->GetAllWeight();

	vector<vector<float> > csk_weight_copy;
	csk_weight_copy.resize(csk_weight.size());
	for(i=0; i<csk_weight.size(); i++){
		csk_weight_copy[i] = vector<float>();
		for(j=0; j<csk_weight[i].size(); j++)
			csk_weight_copy[i].push_back(csk_weight[i][j]);
	}

	const vector< vector<AResultFloat> > &vcsk = csk->GetAllResult();
	vector< vector<string> > vcNew;
	vcNew.resize(vcsk.size());
	for(i=0; i<vcsk.size(); i++){
		vcNew[i] = vector<string>();
		for(j=0; j<TOP; j++){
			vcNew[i].push_back(csk->GetGene(vcsk[i][j].i));
		}
	}
	csk->Destruct();
	delete csk;
*/
/*
	vector< vector<string> > vcIntersect;
	vcIntersect.resize(vcNew.size());
	for(i=0; i<vcNew.size(); i++){
		vcIntersect[i] = vector<string>();
		vector<string> s1, s2;
		vector<string> intersect;
		intersect.resize(TOP);
		vector<string>::iterator it;

		for(j=0; j<origQuery[i].size(); j++)
			vcIntersect[i].push_back(origQuery[i][j]);

		//int G = max((int)1, (int)(origQuery[i].size()*0.3));
		//int G = max((int)1, (int)(20 - origQuery[i].size()));

		for(j=0; j<TOP; j++)
			s1.push_back(vcNew[i][j]);
		for(j=0; j<20; j++)
			s2.push_back(newQuery[i][j]);

		sort(s1.begin(), s1.end());
		sort(s2.begin(), s2.end());

		//fprintf(stderr, "G: %d\n", G);
		//for(j=0; j<TOP; j++){
		//	s1.push_back(vcNew[i][j]);
		//	s2.push_back(newQuery[i][j]);
		//	sort(s1.begin(), s1.end());
		//	sort(s2.begin(), s2.end());
		//	it = set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
		//		intersect.begin());
		//	//if((int)(it - intersect.begin()) > G) break;
		//}
		it = set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
			intersect.begin());

		int size2 = (int) (it - intersect.begin());
		for(j=0; j<size2; j++) vcIntersect[i].push_back(intersect[j]);
	}

	CSeekTools::Write2DArrayText("/tmp/ex_query.txt", vcIntersect);
*/

	//vector<vector<string> > newQ;
	//CSeekTools::ReadMultipleQueries("/tmp/ex_query2.txt", newQ);

	CSeekCentral *csfinal = new CSeekCentral();
	if(!csfinal->Initialize(sArgs.input_arg, sArgs.quant_arg, sArgs.dset_arg,
		sArgs.search_dset_arg, 
		//"/tmp/ex_query2.txt", 
		sArgs.query_arg,
		sArgs.dir_platform_arg,
		sArgs.dir_in_arg, sArgs.dir_prep_in_arg, useNibble, sArgs.num_db_arg,
		sArgs.buffer_arg, "results", sArgs.output_text_flag,  
		!!sArgs.norm_subavg_flag, !!sArgs.norm_platsubavg_flag,
		!!sArgs.norm_platstdev_flag, false,
		sArgs.score_cutoff_arg, sArgs.per_q_required_arg))
		return -1;

	//csfinal->WeightSearch(csk_weight_copy);
	//csfinal->CVCustomSearch(newQ, rnd, PART_M, FOLD, RATE);
	csfinal->CVSearch(rnd, PART_M, FOLD, RATE);
	csfinal->Destruct();
	delete csfinal;

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; 
}
