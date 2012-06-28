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
	size_t				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	if( sArgs.input_arg ) {
		string strGeneInput = sArgs.input_arg;
		vector<string> vecstrGeneID;
		if(!CSeekTools::ReadListTwoColumns(strGeneInput, vecstrGeneID, vecstrGenes)){
			return false;
		}
	}

	bool useNibble = false;
	if(sArgs.is_nibble_flag==1){
		useNibble = true;
	}

	CDatabase DB(useNibble);

	if(sArgs.db_arg){
		string strDBInput = sArgs.db_arg;
		vector<string> vecstrDatasets, vecstrDP;
		if(!CSeekTools::ReadListTwoColumns(strDBInput, vecstrDatasets, vecstrDP)){
			return false;
		}
		map<string, string> mapstrstrDatasetPlatform;
		for(i=0; i<vecstrDatasets.size(); i++){
			mapstrstrDatasetPlatform[vecstrDatasets[i]] = vecstrDP[i];
		}

		string strQueryInput = sArgs.query_arg;
		vector<string> vecstrQuery;
		CSeekStrIntMap mapstriQuery;
		if(!CSeekTools::ReadListOneColumn(strQueryInput, vecstrQuery, mapstriQuery)){
			return false;
		}

		string strPlatformDirectory = sArgs.dir_platform_arg;
		vector<CSeekPlatform> vp;
		map<string, size_t> mapstriPlatform;
		vector<string> vecstrPlatforms;
		CSeekTools::ReadPlatforms(strPlatformDirectory, vp, vecstrPlatforms,
				mapstriPlatform);

		//printf("Done reading"); getchar();

		string strInputDirectory = sArgs.dir_in_arg;
		string strPrepInputDirectory = sArgs.dir_prep_in_arg;
		size_t iNumDBs = sArgs.num_db_arg;
		size_t iDatasets = vecstrDatasets.size();
		size_t iGenes = vecstrGenes.size();

		DB.Open(strInputDirectory, vecstrGenes, iDatasets, iNumDBs);

		vector<CSeekDataset*> vc;
		vector<char> cQuery;
		CSeekTools::LoadDatabase(DB, strPrepInputDirectory, cQuery, vecstrQuery,
			vecstrDatasets, mapstrstrDatasetPlatform, mapstriPlatform, vp, vc);

		/*
		DB.Open(strInputDirectory);
		size_t j,k;
		vc.clear();
		vc.resize(iDatasets);
		for(i=0; i<iDatasets; i++){
			vc[i] = new CSeekDataset();
			string strPrepInputDirectory = sArgs.dir_prep_in_arg;
			string strFileStem = vecstrDatasets[i];
			//string strFileStem = CMeta::Deextension(CMeta::Basename(vecstrDatasets[i].c_str()));
			string strAvgPath = strPrepInputDirectory + "/" + strFileStem + ".gavg";
			string strPresencePath = strPrepInputDirectory + "/" + strFileStem + ".gpres";
			vc[i]->ReadGeneAverage(strAvgPath);
			vc[i]->ReadGenePresence(strPresencePath);
		}

		CSeekTools::InitVector(cQuery, iGenes, (char) 0);

		for(i=0; i<vecstrQuery.size(); i++){
			k = DB.GetGene(vecstrQuery[i]);
			if(k==-1) continue;
			cQuery[k] = 1;
		}

		for(i=0; i<iDatasets; i++){
			vc[i]->InitializeQuery(cQuery);
		}

		vector<unsigned char> *Q = new vector<unsigned char>[vecstrQuery.size()];

		for(i=0; i<vecstrQuery.size(); i++){
			if(!DB.GetGene(vecstrQuery[i], Q[i])){
				cerr << "Gene does not exist" << endl;
			}
		}

		//printf("Before"); getchar();
		for(i=0; i<vecstrQuery.size(); i++){
			if(DB.GetGene(vecstrQuery[i])==-1){
				continue;
			}
			size_t m = DB.GetGene(vecstrQuery[i]);
			size_t l = 0;
			for(j=0; j<iDatasets; j++){
				CSeekIntIntMap *qu = vc[j]->GetQueryMap();
			    size_t query = qu->GetForward(m);
			    if(query==-1) continue;
			    for(k=0; k<iGenes; k++){
			    	unsigned char c = Q[i][k*iDatasets + j];
			    	vc[j]->SetQueryNoMapping(query, k, c);
			    }
			}
		}
		delete[] Q;
		*/
		size_t j;
		float RATE = 0.95;
		int FOLD = 5;
		enum PartitionMode PART_M = CUSTOM_PARTITION;

		const gsl_rng_type *T;
		gsl_rng *rnd;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rnd = gsl_rng_alloc(T);

		int d;
		omp_set_num_threads(8);
		int numThreads = omp_get_max_threads();

		for(i=0; i<1; i++){
			CSeekQuery query;
			query.InitializeQuery(cQuery);
			query.CreateCVPartitions(rnd, PART_M, FOLD);

			vector<float> master_rank;
			CSeekTools::InitVector(master_rank, iGenes, (float) 0);

			vector<float> sum_weight;
			CSeekTools::InitVector(sum_weight, iGenes, (float) 0);

			vector<int> counts;
			CSeekTools::InitVector(counts, iGenes, (int) 0);

			float **master_rank_threads = CSeekTools::Init2DArray(numThreads, iGenes, (float) 0);
			float **sum_weight_threads = CSeekTools::Init2DArray(numThreads, iGenes, (float) 0);
			int **counts_threads = CSeekTools::Init2DArray(numThreads, iGenes, (int) 0);

			printf("Entering search\n");

			#pragma omp parallel for \
			shared(vc) \
			private(d, j) \
			firstprivate(iDatasets) \
			schedule(static)

			for(d=0; d<iDatasets; d++){
				int tid = omp_get_thread_num();
				//printf("Dataset %d\n", d);
				CSeekIntIntMap *mapQ = vc[d]->GetQueryMap();
				CSeekIntIntMap *mapG = vc[d]->GetGeneMap();

				vector<int> this_q;
				for(j=0; j<mapQ->GetNumSet(); j++){
					this_q.push_back(mapQ->GetReverse(j));
				}

				if(mapQ->GetNumSet()==0){
					//printf("This dataset is skipped\n");
					continue;
				}

				//printf("Initializing\n");
				vc[d]->InitializeFloatMatrix();
				//printf("Weighting dataset\n");
				CSeekWeighter::CVWeighting(query, *vc[d]);
				float w = vc[d]->GetDatasetSumWeight();
				if(w==-1){
					//printf("Bad weight\n");
					vc[d]->FreeFloatMatrix();
					continue;
					//getchar();
				}
				vector<float> rank_normal;
				//printf("Doing linear combination\n");
				CSeekWeighter::LinearCombine(rank_normal, this_q, *vc[d]);
				/*for(j=0; j<1000; j++){
					size_t g = mapG->GetReverse(j);
					printf("Gene %d %.5f\n", g, rank_normal[g]);
				}*/
				vc[d]->FreeFloatMatrix();

				//printf("Adding contribution of dataset to master ranking: %.5f\n", w);
				for(j=0; j<mapG->GetNumSet(); j++){
					size_t g = mapG->GetReverse(j);
					master_rank_threads[tid][g] += rank_normal[g] * w;
					counts_threads[tid][g]++;
					sum_weight_threads[tid][g] += w;
				}
			}

			for(j=0; j<numThreads; j++){
				size_t k;
				for(k=0; k<iGenes; k++){
					master_rank[k] += master_rank_threads[j][k];
					counts[k] += counts_threads[j][k];
					sum_weight[k]+=sum_weight_threads[j][k];
				}
			}

			CSeekTools::Free2DArray(master_rank_threads);
			CSeekTools::Free2DArray(counts_threads);
			CSeekTools::Free2DArray(sum_weight_threads);


			printf("Aggregating genes\n");
			for(j=0; j<iGenes; j++){
				if(counts[j]<(int)(0.5*iDatasets)){
					master_rank[j] = -50.0;
				}else if(sum_weight[j]==0){
					master_rank[j] = -50.0;
				}else{
					master_rank[j] /= sum_weight[j];
				}
				//printf("Gene %d %.5f\n", j, master_rank[j]);
			}

			printf("Sorting genes\n");
			vector<AResult> a;
			a.clear();
			a.resize(iGenes);
			for(j=0; j<iGenes; j++){
				a[j].i = j;
				a[j].f = master_rank[j];
			}
			printf("Begin Sorting genes\n");
			sort(a.begin(), a.end());

			printf("Results:\n");
			size_t jj;
			size_t ii;
			for(ii=0, jj=0; jj<500; ii++){
				//if(cQuery[a[ii].i]==1) continue;
				printf("%d %.5f\n", a[ii].i, a[ii].f);
				jj++;
			}
		}

		/*for(i=0; i<iDatasets; i++){
			printf("Dataset %ld\n", i);
			CSeekMatrix<unsigned char> *cm = vc[i]->GetMatrix();
			for(j=0; j<cm->GetNumRow(); j++){
				printf("Row %ld\n", j);
				for(k=0; k<1000; k++){
					printf("%d ", cm->Get(j, k));
				}
				printf("\n");
			}
		}*/
		/*size_t j;
		for(i=0; i<vecstrQuery.size(); i++){
			printf("Query: %s\n", vecstrQuery[i].c_str());
			for(j=0; j<Q[i].size(); j++){
				printf("%d ", (int) Q[i][j]);
			}
			printf("\n");
			getchar();
		}*/

		//printf("Done"); getchar();

	}else{
		cerr << "Must give a db list." << endl;
		return 1;

	}

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }
