/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
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
	const int lineSize = 1024;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		//cmdline_parser_print_help( );
		fprintf(stderr, "Use -h to get help.\n");
		return 1;
	}

	string method = sArgs.weighting_method_arg;
	string cv = sArgs.CV_partition_arg;
	int cv_fold = sArgs.CV_fold_arg;
	float rbp_p = sArgs.CV_rbp_p_arg;
	string dist_measure = sArgs.dist_measure_arg;

	if(dist_measure=="pearson"){
		string sinfo_dir = sArgs.dir_sinfo_arg;
		if(sinfo_dir=="NA"){
			fprintf(stderr, "Pearson selected. Please give the sinfo directory (-u).\n");
			return 1;
		}
		if(!!sArgs.norm_subavg_flag){
			fprintf(stderr, "Warning: -m flag is ignored due to --dist_measure pearson.\n");
		}
		if(!!sArgs.norm_subavg_plat_flag){
			fprintf(stderr, "Warning: -M flag is ignored due to --dist_measure pearson.\n");
		}
		
	}else if(dist_measure=="z_score"){
		if(!!sArgs.norm_subavg_plat_flag && !(!!sArgs.norm_subavg_flag)){
			fprintf(stderr, "Please enable -m flag. --norm_subavg_plat requires -m flag.\n");
			return 1;
		}
	}

	if(sArgs.check_dset_size_flag==1){
		string dsize_file = sArgs.dset_size_file_arg;
		if(dsize_file=="NA"){
			fprintf(stderr, "Dataset size file is missing\n");
			return 1;
		}
	}

	if(!sArgs.input_arg || !sArgs.quant_arg ||
		!sArgs.dset_arg ||
		!sArgs.query_arg || !sArgs.dir_platform_arg ||
		!sArgs.dir_in_arg || !sArgs.dir_prep_in_arg){
		fprintf(stderr, "Arguments missing!\n");
		return 1;
	}

	bool useNibble = false;
	if(sArgs.is_nibble_flag==1){
		fprintf(stderr, "Nibble integration is not supported! Please use a non-nibble CDatabase.\n");
		useNibble = true;
		return 1;
	}

	bool bOutputWeightComponent = !!sArgs.output_w_comp_flag;
	bool bSimulateWeight = !!sArgs.simulate_w_flag;

	// Random Number Generator Initializations
	gsl_rng_env_setup();

	const gsl_rng_type *T;
	T = gsl_rng_default;

	gsl_rng *rnd = gsl_rng_alloc(T);

	const gsl_rng_type *T2;
	T2 = gsl_rng_default;

	gsl_rng *random_ranking_rnd = gsl_rng_alloc(T2);

	float RATE = rbp_p;
	utype FOLD = (utype) cv_fold;
	enum CSeekQuery::PartitionMode PART_M;
	if(cv=="LOI"){
		PART_M = CSeekQuery::LEAVE_ONE_IN;
	}else if(cv=="LOO"){
		PART_M = CSeekQuery::LEAVE_ONE_OUT;
	}else if(cv=="XFOLD"){
		PART_M = CSeekQuery::CUSTOM_PARTITION;
	}

	utype i,j;
	//utype TOP = 1000;
	//utype TOP = 0;

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
		const vector<utype> &queryGenes = vq[i].GetQuery();
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


	enum CSeekDataset::DistanceMeasure eDistMeasure;
	if(dist_measure=="pearson"){
		eDistMeasure = CSeekDataset::CORRELATION;
	}else{
		eDistMeasure = CSeekDataset::Z_SCORE;
	}

	/*fprintf(stderr, "input: %s quant: %s dset: %s, search_dset: %s\n", sArgs.input_arg, sArgs.quant_arg, sArgs.dset_arg, sArgs.search_dset_arg);
	fprintf(stderr, "query: %s dir_plat: %s dir_in: %s, dir_prep: %s\n", sArgs.query_arg, sArgs.dir_platform_arg, sArgs.dir_in_arg, sArgs.dir_prep_in_arg);
	fprintf(stderr, "dir_gvar: %s dir_sinfo: %s useNibble: %d, num_db: %s\n", sArgs.dir_gvar_arg, sArgs.dir_sinfo_arg, sArgs.is_nibble_flag, sArgs.num_db_arg);
	getchar();*/

	CSeekCentral *csfinal = new CSeekCentral();
	CSeekDBSetting *dbSetting = new CSeekDBSetting(sArgs.dir_gvar_arg,
		sArgs.dir_sinfo_arg, sArgs.dir_platform_arg, sArgs.dir_prep_in_arg,
		sArgs.dir_in_arg, sArgs.input_arg, sArgs.quant_arg, sArgs.dset_arg,
		sArgs.dset_size_file_arg, sArgs.num_db_arg);
	vector<CSeekDBSetting*> cc;
	cc.push_back(dbSetting);

	string add_db = sArgs.additional_db_arg;
	if(add_db!="NA"){
		ifstream ifsm;
		ifsm.open(add_db.c_str());
		if(!ifsm.is_open()){
			fprintf(stderr, "Error opening file %s\n", add_db.c_str());
			return false;
		}
		char acBuffer[lineSize];
		utype c_iBuffer = lineSize;
		vector<map<string,string> > parameters; //an array of CDatabase's
		while(!ifsm.eof()){
			ifsm.getline(acBuffer, c_iBuffer-1);
			if(acBuffer[0]==0) break;
			acBuffer[c_iBuffer-1]=0;
			string strB = acBuffer;
			if(strB=="START"){
				map<string,string> p;
				while(!ifsm.eof()){
					ifsm.getline(acBuffer, c_iBuffer-1);
					if(acBuffer[0]==0){
						fprintf(stderr, "Invalid line (empty)\n");
						return 1;
					}
					strB = acBuffer;
					if(strB=="END") break;
					vector<string> tok;
					CMeta::Tokenize(acBuffer, tok); //separator is tab
					p[tok[0]] = tok[1];
				}
				parameters.push_back(p);
			}
		}
		ifsm.close();
		if(parameters.size()==0){
			fprintf(stderr, "Error, extra_db setting file must begin with START and end with END lines\n");
			return 1;
		}

		//i=0;
		for(i=0; i<parameters.size(); i++){
		string sinfo_dir = "NA";
		string gvar_dir = "NA";
		string platform_dir = "NA";
		string prep_dir = "NA";
		string db_dir = "NA";
		string dset_map_file = "NA";
		string gene_map_file = "NA";
		string quant_file = "NA";
		string dset_size_file = "NA";
		int num_db = -1;

		if(eDistMeasure==CSeekDataset::CORRELATION){
			if(parameters[i].find("SINFO_DIR")==parameters[i].end() ||
				parameters[i].find("SINFO_DIR")->second=="NA"){
				fprintf(stderr, "Please specify an sinfo directory for the extra db\n");
				return false;
			}
			sinfo_dir = parameters[i].find("SINFO_DIR")->second;
		}
		if(parameters[i].find("GVAR_DIR")!=parameters[i].end())
			gvar_dir = parameters[i].find("GVAR_DIR")->second;

		if(sArgs.check_dset_size_flag==1){
			if(parameters[i].find("DSET_SIZE_FILE")==parameters[i].end() ||
				parameters[i].find("DSET_SIZE_FILE")->second=="NA"){
				fprintf(stderr, "Please specify the dataset size file for the extra db\n");
				return false;
			}
		}

		if(parameters[i].find("DSET_SIZE_FILE")!=parameters[i].end() &&
			parameters[i].find("DSET_SIZE_FILE")->second!="NA")
			dset_size_file = parameters[i].find("DSET_SIZE_FILE")->second;
		
		if(parameters[i].find("PREP_DIR")==parameters[i].end() ||
			parameters[i].find("PLATFORM_DIR")==parameters[i].end() ||
			parameters[i].find("DB_DIR")==parameters[i].end() ||
			parameters[i].find("DSET_MAP_FILE")==parameters[i].end() ||
			parameters[i].find("GENE_MAP_FILE")==parameters[i].end() ||
			parameters[i].find("QUANT_FILE")==parameters[i].end() ||
			parameters[i].find("NUMBER_OF_DB")==parameters[i].end()){
			fprintf(stderr, "Some arguments are missing. Please make sure the following are provided:\n");
			fprintf(stderr, "PREP_DIR, DB_DIR, DSET_MAP_FILE, GENE_MAP_FILE, QUANT_FILE, NUMBER_OF_DB\n");
		}

		platform_dir = parameters[i].find("PLATFORM_DIR")->second;
		db_dir = parameters[i].find("DB_DIR")->second;
		prep_dir = parameters[i].find("PREP_DIR")->second;
		dset_map_file = parameters[i].find("DSET_MAP_FILE")->second;
		gene_map_file = parameters[i].find("GENE_MAP_FILE")->second;
		quant_file = parameters[i].find("QUANT_FILE")->second;
		num_db = atoi(parameters[i].find("NUMBER_OF_DB")->second.c_str());

		CSeekDBSetting *dbSetting2 = new CSeekDBSetting(gvar_dir, sinfo_dir,
			platform_dir, prep_dir, db_dir, gene_map_file, quant_file, dset_map_file,
			dset_size_file, num_db);
		cc.push_back(dbSetting2);
		}
	}

	bool bVariance = false;
	if(method=="VAR"){
		bVariance = true;
	}

	if(sArgs.per_q_required_arg>=1.00001 || sArgs.per_q_required_arg<=-0.00001){
		fprintf(stderr, "Error, per_q_required needs to be <=1.0 and >=0.0\n");
		return -1;
	}

	if(sArgs.per_g_required_arg>=1.00001 || sArgs.per_g_required_arg<=-0.00001){
		fprintf(stderr, "Error, per_g_required needs to be <=1.0 and >=0.0\n");
		return -1;
	}

	if(!csfinal->Initialize(cc,
		sArgs.search_dset_arg, 
		//"/tmp/ex_query2.txt", 
		sArgs.query_arg,
		sArgs.output_dir_arg,
		sArgs.buffer_arg, !!sArgs.output_text_flag,
		bOutputWeightComponent, bSimulateWeight,
		eDistMeasure, bVariance,
		!!sArgs.norm_subavg_flag, !!sArgs.norm_subavg_plat_flag,
		false,
		!!sArgs.check_dset_size_flag,
		sArgs.score_cutoff_arg, 
		sArgs.per_q_required_arg, sArgs.per_g_required_arg,
		!!sArgs.square_z_flag,
		!!sArgs.random_flag, sArgs.num_random_arg, !!sArgs.neg_cor_flag,
		random_ranking_rnd, useNibble, 
		sArgs.num_threads_arg))
		return -1;

	if(method=="CV"){
		csfinal->CVSearch(rnd, PART_M, FOLD, RATE);
	}else if(method=="EQUAL"){
		csfinal->EqualWeightSearch();
	}else if(method=="ORDER_STAT"){
		csfinal->OrderStatistics();
	}else if(method=="USER"){
		string uw = sArgs.user_weight_list_arg;
		vector<string> uww;
		if(!CSeekTools::ReadListOneColumn(uw.c_str(), uww)){
			fprintf(stderr, "Error reading user weight list\n");
			return -1;
		}
		vector<vector<float> > fw;
		fw.resize(uww.size());
		for(i=0; i<uww.size(); i++){
			if(!CSeekTools::ReadArray(uww[i].c_str(), fw[i])){
				return -1;
			}
		}
		csfinal->WeightSearch(fw);
	}else if(method=="VAR"){
		for(i=0; i<cc.size(); i++){
			CSeekDBSetting* pc = cc[i];
			if(pc->GetValue("gvar")=="NULL"){
				fprintf(stderr, "Must specify gvar directory!\n");
				return -1;
			}
		}
		if(bSimulateWeight){
			fprintf(stderr, "simulate weight option is not supported for variance-based weighting\n");
			return -1;
		}
		csfinal->VarianceWeightSearch();
	}else if(method=="AVERAGE_Z"){
		csfinal->AverageWeightSearch();
	}else if(method=="CV_CUSTOM"){
		string uw = sArgs.user_gene_list_arg;
		vector<vector<string> > user_gene_list;
		if(!CSeekTools::ReadMultipleQueries(uw, user_gene_list, 2048)){
			fprintf(stderr, "Error reading user gene lists!\n");
			return -1;
		}	
		csfinal->CVCustomSearch(user_gene_list, rnd, PART_M, FOLD, RATE);
	}
	//csfinal->WeightSearch(csk_weight_copy);
	//csfinal->CVCustomSearch(newQ, rnd, PART_M, FOLD, RATE);
	//csfinal->EqualWeightSearch();
	//csfinal->CVSearch(rnd, PART_M, FOLD, RATE);
	//csfinal->OrderStatistics();
	fprintf(stderr, "Destructing...\n");
	csfinal->Destruct();
	fprintf(stderr, "Deleting...\n");
	delete csfinal;

	fprintf(stderr, "Deleting DBSetting...\n");
	//if(add_db!="NA"){
		for(i=0; i<cc.size(); i++){
			delete cc[i];
		}
	//}
	fprintf(stderr, "Finished deleting DBSetting...\n");
	cc.clear();
	gsl_rng_free(rnd);
	gsl_rng_free(random_ranking_rnd);

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; 
}
