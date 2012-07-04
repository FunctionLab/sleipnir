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
	omp_set_num_threads(8);

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

		fprintf(stderr, "Start reading platform\n"); system("date +%s%N 1>&2");
		string strPlatformDirectory = sArgs.dir_platform_arg;
		vector<CSeekPlatform> vp;
		map<string, ushort> mapstriPlatform;
		vector<string> vecstrPlatforms;
		CSeekTools::ReadPlatforms(strPlatformDirectory, vp, vecstrPlatforms,
				mapstriPlatform);

		fprintf(stderr, "Done reading platform\n"); system("date +%s%N 1>&2");

		string strInputDirectory = sArgs.dir_in_arg;
		string strPrepInputDirectory = sArgs.dir_prep_in_arg;
		size_t iNumDBs = sArgs.num_db_arg;
		size_t iDatasets = vecstrDatasets.size();
		size_t iGenes = vecstrGenes.size();

		fprintf(stderr, "Start reading CDatabase header\n"); system("date +%s%N 1>&2");
		DB.Open(strInputDirectory, vecstrGenes, iDatasets, iNumDBs);
		fprintf(stderr, "Done reading CDatabase header\n"); system("date +%s%N 1>&2");

		//printf("Done opening"); getchar();

		vector<CSeekDataset*> vc;
		vector<char> cQuery;
		CSeekTools::LoadDatabase(DB, strPrepInputDirectory, cQuery, vecstrQuery,
			vecstrDatasets, mapstrstrDatasetPlatform, mapstriPlatform, vp, vc);

		ushort j;
		ushort d;
		float RATE = 0.95;
		ushort FOLD = 5;
		enum PartitionMode PART_M = CUSTOM_PARTITION;

		const gsl_rng_type *T;
		gsl_rng *rnd;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rnd = gsl_rng_alloc(T);

		ushort numThreads = omp_get_max_threads();

		for(i=0; i<1; i++){
			CSeekQuery query;
			query.InitializeQuery(cQuery);
			query.CreateCVPartitions(rnd, PART_M, FOLD);
			ushort iQuery = query.GetQuery().size();

			ushort ***rData = new ushort**[numThreads];
			for(j=0; j<numThreads; j++){
				rData[j] = CSeekTools::Init2DArray(iGenes, iQuery, (ushort)0);
			}

			vector<float> master_rank;
			CSeekTools::InitVector(master_rank, iGenes, (float) 0);

			vector<float> sum_weight;
			CSeekTools::InitVector(sum_weight, iGenes, (float) 0);

			vector<ushort> counts;
			CSeekTools::InitVector(counts, iGenes, (ushort) 0);

			float **master_rank_threads = CSeekTools::Init2DArray(numThreads, iGenes, (float) 0);
			float **sum_weight_threads = CSeekTools::Init2DArray(numThreads, iGenes, (float) 0);
			ushort **counts_threads = CSeekTools::Init2DArray(numThreads, iGenes, (ushort) 0);

			fprintf(stderr, "Entering search\n"); system("date +%s%N 1>&2");

			#pragma omp parallel for \
			shared(vc, rData) \
			private(d, j) \
			firstprivate(iDatasets) \
			schedule(dynamic)

			for(d=0; d<iDatasets; d++){
				ushort tid = omp_get_thread_num();
				//printf("Dataset %d\n", d);
				CSeekIntIntMap *mapQ = vc[d]->GetQueryMap();
				CSeekIntIntMap *mapG = vc[d]->GetGeneMap();

				vector<ushort> this_q;
				for(j=0; j<mapQ->GetNumSet(); j++){
					this_q.push_back(mapQ->GetReverse(j));
				}

				if(mapQ->GetNumSet()==0){
					//printf("This dataset is skipped\n");
					continue;
				}

				//printf("Initializing\n");
				vc[d]->InitializeDataMatrix(rData[tid], iGenes, iQuery);
				//printf("Weighting dataset\n");
				CSeekWeighter::CVWeighting(query, *vc[d]);
				float w = vc[d]->GetDatasetSumWeight();
				if(w==-1){
					//printf("Bad weight\n");
					//vc[d]->FreeDataMatrix();
					continue;
					//getchar();
				}
				vector<ushort> rank_normal;
				//printf("Doing linear combination\n");
				CSeekWeighter::LinearCombine(rank_normal, this_q, *vc[d]);
				/*for(j=0; j<1000; j++){
					size_t g = mapG->GetReverse(j);
					printf("Gene %d %d\n", g, rank_normal[g]);
				}*/
				//vc[d]->FreeDataMatrix();

				//printf("Adding contribution of dataset to master ranking: %.5f\n", w);

				//int BAD = 0;
				for(j=0; j<mapG->GetNumSet(); j++){
					ushort g = mapG->GetReverse(j);
					if(rank_normal[g]==0){
						//printf("BAD %d\n", g);
						//BAD++;
						continue;
					}
					master_rank_threads[tid][g] += (float) rank_normal[g] * w;
					sum_weight_threads[tid][g] += w;
					counts_threads[tid][g]++;
				}

				/*if(BAD==0){
					printf("Good!\n");
				}else{
					printf("BAD: %d %d\n", BAD, mapG->GetNumSet());
				}*/
			}

			for(j=0; j<numThreads; j++){
				ushort k;
				for(k=0; k<iGenes; k++){
					master_rank[k] += master_rank_threads[j][k];
					counts[k] += counts_threads[j][k];
					sum_weight[k]+=sum_weight_threads[j][k];
				}
			}

			CSeekTools::Free2DArray(master_rank_threads);
			CSeekTools::Free2DArray(counts_threads);
			CSeekTools::Free2DArray(sum_weight_threads);

			for(j=0; j<numThreads; j++){
				CSeekTools::Free2DArray(rData[j]);
			}
			delete[] rData;

			fprintf(stderr, "Aggregating genes\n");
			for(j=0; j<iGenes; j++){
				if(counts[j]<(int)(0.5*iDatasets)){
					master_rank[j] = -320;
				}else if(sum_weight[j]==0){
					master_rank[j] = -320;
				}else{
					master_rank[j] = (master_rank[j] / sum_weight[j] - 320) / 100.0;
				}
				//printf("Gene %d %.5f\n", j, master_rank[j]);
			}

			fprintf(stderr, "Sorting genes\n");
			vector<AResultFloat> a;
			a.clear();
			a.resize(iGenes);
			for(j=0; j<iGenes; j++){
				a[j].i = j;
				a[j].f = master_rank[j];
			}
			fprintf(stderr, "Begin Sorting genes\n");
			sort(a.begin(), a.end());

			fprintf(stderr, "Results:\n");
			ushort jj;
			ushort ii;
			for(ii=0, jj=0; jj<500; ii++){
				if(cQuery[a[ii].i]==1) continue;
				printf("%d %.5f\n", a[ii].i, a[ii].f);
				jj++;
			}

			fprintf(stderr, "Done search\n"); system("date +%s%N 1>&2");
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
