#include "stdafx.h"
#include "cmdline.h"

struct bresult{
	int i;
	int j;
	float f;
	float f2;
};

struct ascend{
	bool operator()(const bresult& lx, const bresult& rx) const{
		return lx.f < rx.f;
	}
};

struct descend{
	bool operator()(const bresult& lx, const bresult& rx) const{
		return lx.f > rx.f;
	}
};


bool calculate_correlation(
	vector<vector<vector<float> > > &mat,
	vector<CSeekIntIntMap*> &dm,
	utype g1, 
	utype g2, 
	float &r){

	int numDatasets = mat.size();
	int d, i;
	vector<float> v1, v2;

	vector<int> count_common;
	CSeekTools::InitVector(count_common, numDatasets, (int) 0);

	for(d=0; d<numDatasets; d++){
		if(CSeekTools::IsNaN(dm[d]->GetForward(g1))) continue;
		count_common[d]++;
	}

	for(d=0; d<numDatasets; d++){
		if(CSeekTools::IsNaN(dm[d]->GetForward(g2))) continue;
		count_common[d]++;
	}

	int total_common = 0;
	for(d=0; d<numDatasets; d++){
		if(count_common[d]==2) total_common++;
	}

	/*if(total_common<0.5*(float)numDatasets){
		r = -320;
		return true;
	}*/

	for(d=0; d<numDatasets; d++){
		if(count_common[d]==2){
			utype n1 = dm[d]->GetForward(g1);
			utype n2 = dm[d]->GetForward(g2);
			int numExp = mat[d][n1].size();
			for(i=0; i<numExp; i++){
				v1.push_back(mat[d][n1][i]);
				v2.push_back(mat[d][n2][i]);
				/*if(isinf(mat[d][n1][i])||isnan(mat[d][n1][i])){
					fprintf(stderr, "BAD: d %d g %d e %d v %.2f\n", d, n1, i, mat[d][n1][i]);
				}
				if(isinf(mat[d][n2][i])||isnan(mat[d][n2][i])){
					fprintf(stderr, "BAD: d %d g %d e %d v %.2f\n", d, n2, i, mat[d][n2][i]);
				}*/
			}
		}
	}

	//pearson correlation calculation
	double mean_x = 0;
	double mean_y = 0;
	int xx = 0;
	for(xx=0; xx<v1.size(); xx++){
		mean_x+=v1[xx];
		mean_y+=v2[xx];
	}
	mean_x /= (double) v1.size();
	mean_y /= (double) v2.size();
	
	double sum_xy = 0;
	double dev_x = 0;
	double dev_y = 0;
	for(xx=0; xx<v1.size(); xx++){
		sum_xy += (v1[xx] - mean_x) * (v2[xx] - mean_y);
		dev_x += (v1[xx] - mean_x) * (v1[xx] - mean_x);
		dev_y += (v2[xx] - mean_y) * (v2[xx] - mean_y);
	}
	dev_x = sqrt(dev_x);
	dev_y = sqrt(dev_y);
	r = (float) (sum_xy / dev_x / dev_y);
	if(isinf(r) || isnan(r)){
		r = -320;
	}
	//fprintf(stderr, "%.2f\t%d\n", r, total_common);

	return true;
			
}

int main(int iArgs, char **aszArgs){
	static const size_t c_iBuffer   = 1024;
	char acBuffer[c_iBuffer];
	
	gengetopt_args_info sArgs;
	ifstream ifsm;
	istream *pistm;
	vector<string> vecstrLine, vecstrGenes, vecstrDBs, vecstrQuery;
	utype i, qi, j, k, l;
	
	if(cmdline_parser(iArgs, aszArgs, &sArgs)){
		cmdline_parser_print_help();
		return 1;
	}

	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		pistm = &ifsm; 
	}else{
		pistm = &cin;
	}

	fprintf(stderr, "Reading gene list\n");	
	map<string, size_t> mapstriGenes;
	while( !pistm->eof( ) ) {
		pistm->getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) < 2 ) {
			cerr << "Ignoring line: " << acBuffer << endl;
			continue;
		}
		if( !( i = atoi( vecstrLine[ 0 ].c_str( ) ) ) ) {
			cerr << "Illegal gene ID: " << vecstrLine[ 0 ] << " for "
				<< vecstrLine[ 1 ] << endl;
			return 1;
		}
		i--;
		if( vecstrGenes.size( ) <= i )
			vecstrGenes.resize( i + 1 );
		vecstrGenes[ i ] = vecstrLine[ 1 ];
		mapstriGenes[vecstrGenes[i]] = i;
	}
			
	//char acBuffer[1024];

	fprintf(stderr, "Finished reading gene map\n");
	if(sArgs.pcl_flag==1){
		string pcl_dir = sArgs.pcl_dir_arg;
		string output_dir = sArgs.dir_out_arg;
		fprintf(stderr, "Arrived here\n");
		vector<string> pcl_list;
		vector< vector<string> > vecstrAllQuery;
		int numGenes = vecstrGenes.size();
	
		fprintf(stderr, "Reading query\n");	
		if(!CSeekTools::ReadMultipleQueries(sArgs.query_arg, vecstrAllQuery))
			return -1;
		
		fprintf(stderr, "Reading pcl list\n");	
		CSeekTools::ReadListOneColumn(sArgs.pcl_list_arg, pcl_list);
		vector< vector< vector<float> > > mat; //only needed for cor-calc
		vector<CSeekIntIntMap*> dm;
		dm.resize(pcl_list.size());
		mat.resize(pcl_list.size()); //only needed for cor-calc
		for(i=0; i<pcl_list.size(); i++){
			fprintf(stderr, "Reading %d: %s\n", i, pcl_list[i].c_str());
			dm[i] = new CSeekIntIntMap(vecstrGenes.size());
			
			string pclfile = pcl_dir + "/" + pcl_list[i] + ".bin";

			CPCL pcl;
			pcl.Open(pclfile.c_str());
			int totNumExperiments = pcl.GetExperiments();
			
			vector<utype> presentIndex;
			vector<string> presentGeneNames;
			for(j=0; j<vecstrGenes.size(); j++){
				utype g = pcl.GetGene(vecstrGenes[j]);
				if(CSeekTools::IsNaN(g)) continue; //gene does not exist in the dataset
				presentIndex.push_back(g);
				presentGeneNames.push_back(vecstrGenes[j]);
				dm[i]->Add(j);
			}
			
			//only needed for correlation calc
			mat[i].resize(presentIndex.size());
			for(j=0; j<presentIndex.size(); j++)
				mat[i][j].resize(totNumExperiments);
			
			for(j=0; j<presentIndex.size(); j++){
				float *val = pcl.Get(presentIndex[j]);
				for(k=0; k<pcl.GetExperiments(); k++){
					mat[i][j][k] = val[k];	
					if(isinf(val[k])||isnan(val[k])){
						fprintf(stderr, "loading error: %d of %d\n", k, pcl.GetExperiments());
					}
				}
			}

		}
	
		float present_cutoff = 0.5;
		vector<char> absent_g;
		CSeekTools::InitVector(absent_g, numGenes, (char) 0);

		for(i=0; i<numGenes; i++){
			int totAbsent = 0;
			for(j=0; j<pcl_list.size(); j++){
				CSeekIntIntMap *mapG = dm[j];
				if(CSeekTools::IsNaN(mapG->GetForward(i))){
					totAbsent++;
				}
			}
			if(totAbsent<0.5*(float)pcl_list.size()){
				absent_g[i] = 0;
			}else{
				absent_g[i] = 1;
			}
		}

		CSeekIntIntMap *geneMap	= new CSeekIntIntMap(numGenes);
		for(i=0; i<numGenes; i++){
			if(absent_g[i]==1) continue;
			geneMap->Add(i);
		}

		int numActualGenes = geneMap->GetNumSet();

		
		//Do pair per machine, Step 2============================
/*
		vector<utype> pairs;
		CSeekTools::ReadArray("/tmp/pairs_to_do", pairs);

		int numP = pairs.size()/2;
		int pi=0;
		const vector<utype> &allGenes = geneMap->GetAllReverse();
		vector<float> cor;
		cor.resize(numP);

		#pragma omp parallel for \
		shared(allGenes, cor, mat, dm, pairs) \
		private(pi) \
		firstprivate(numP) schedule(dynamic)
		for(pi=0; pi<numP; pi++){
			utype g1 = allGenes[pairs[pi*2]];
			utype g2 = allGenes[pairs[pi*2+1]];
			if(g1==g2){
				cor[pi] = -320;
				continue;
			}
			calculate_correlation(mat, dm, g1, g2, cor[pi]);
			if(pi%1000==0){
				fprintf(stderr, "  %d of %d\n", pi, numP);
			}
		}

		CSeekTools::WriteArray("/tmp/results_pairs", cor);		
*/		
		//=============================================================	
		
		
	
		vector< vector<float> > correlations;
		correlations.resize(numActualGenes);
		for(i=0; i<numActualGenes; i++){
			correlations[i].resize(numActualGenes);
		}

		//combining pairs, Step 3=================================
	/*
		int max_i = 49;
		char file[256];
		const vector<utype> &allGenes = geneMap->GetAllReverse();

		for(i=0; i<=max_i; i++){
			vector<utype> pairs;
			vector<float> cor;

			sprintf(file, "/memex/qzhu/p1/pairs_to_do.%d", i);
			CSeekTools::ReadArray(file, pairs);
			sprintf(file, "/memex/qzhu/p1/oct6/pairs_to_do.%d_results", i);
			CSeekTools::ReadArray(file, cor);

			int numP = pairs.size()/2;
			int pi=0;

			for(pi=0; pi<numP; pi++){
				correlations[pairs[pi*2]][pairs[pi*2+1]] = cor[pi];
				correlations[pairs[pi*2+1]][pairs[pi*2]] = cor[pi];
			}
		}

		vector<float> correlation1D;
		correlation1D.resize(numActualGenes*numActualGenes);
		int kk=0;
		for(i=0; i<numActualGenes; i++){
			for(j=0; j<numActualGenes; j++){
				correlation1D[kk] = correlations[i][j];
				kk++;
			}
		}
		CSeekTools::WriteArray("/memex/qzhu/p1/aggregated_dataset_correlation", correlation1D);

		fprintf(stderr, "Finished creating file\n");
		getchar();
*/
		//============================================================


		/*
		//generate pairs per machine, Step 1================================
		vector<utype> pairs;
		int numP = numActualGenes*(numActualGenes-1);
		pairs.resize(numP);
		int ki = 0;
		for(i=0; i<numActualGenes; i++){
			for(j=i+1; j<numActualGenes; j++){
				pairs[ki] = i; 
				ki++;
				pairs[ki] = j;
				ki++;
			}
		}

		if(ki!=numP){
			fprintf(stderr, "Cannot continue!\n");
			return -1;
		}

		int num_pairs_per_file = (numP / 2) / 49;
		
		ki = 0;
		vector<utype> pp;
		pp.resize(num_pairs_per_file*2);
		char destfile[256];
		int kj = 0;
		int ii = 0;

		fprintf(stderr, "Numpairs per file: %d\n", num_pairs_per_file);
		for(ii=0; ii<numP/2; ii++){
			if(ii%num_pairs_per_file==0 && ii>0){
				sprintf(destfile, "/memex/qzhu/p1/pairs_to_do.%d", ki);
				CSeekTools::WriteArray(destfile, pp);
				pp.clear();
				pp.resize(num_pairs_per_file*2);
				ki++;
				kj = 0;
			}
			pp[kj] = pairs[ii*2];
			pp[kj+1] = pairs[ii*2+1];
			kj+=2;
			if(ii==numP/2-1){
				sprintf(destfile, "/memex/qzhu/p1/pairs_to_do.%d", ki);
				pp.resize(kj);
				CSeekTools::WriteArray(destfile, pp);
			}
		}
		

		fprintf(stderr, "Finished\n");	
		getchar();
		//=====================================================================
		*/

		//Load aggregated dataset, Step 4=====================================

		vector<float> correlation1D;
		CSeekTools::ReadArray("/home/qzhu/Seek/aggregated_dataset_correlation", correlation1D);
		int pi = 0;
		for(i=0; i<numActualGenes; i++){
			for(j=0; j<numActualGenes; j++){
				correlations[i][j] = correlation1D[pi];
				pi++;
			}
		}
		correlation1D.clear();

		//======================================================================


		/*
		const vector<utype> &allGenes = geneMap->GetAllReverse();
		//start correlation calculations
		for(i=0; i<numActualGenes; i++){
			utype g1 = allGenes[i];
			fprintf(stderr, "Gene %d of %d: %d\n", i, numActualGenes, numActualGenes - (i+1));

			#pragma omp parallel for \
			shared(allGenes, correlations, mat, dm) \
			private(j) \
			firstprivate(numActualGenes, g1, i) schedule(dynamic)
			for(j=i+1; j<numActualGenes; j++){
				utype tid = omp_get_thread_num();
				utype g2 = allGenes[j];
				if(g1==g2){
					correlations[i][j] = -320; //Null value
					correlations[j][i] = -320;
					continue;
				}
				float r = -320;
				calculate_correlation(mat, dm, g1, g2, r);
				correlations[i][j] = r;
				correlations[j][i] = r;
				//if((j-(i+1))%500==0){
				//	fprintf(stderr, "   %d of %d from %d\n", j - (i+1), numActualGenes - (i+1), tid);
				//}
			}
			fprintf(stderr, "Done\n");
		}

		vector<float> correlation1D;
		correlation1D.resize(numActualGenes*numActualGenes);
		int kk=0;
		for(i=0; i<numActualGenes; i++){
			for(j=0; j<numActualGenes; j++){
				correlation1D[kk] = correlations[i][j];
				kk++;
			}
		}
		CSeekTools::WriteArray("/memex/qzhu/p1/aggregated_dataset_correlation", correlation1D);
		*/

		//Evaluation Last step==================================================
		
		const vector<utype> &allGenes = geneMap->GetAllReverse();		

		for(i=0; i<vecstrAllQuery.size(); i++){
			vector<utype> q;
			vector<char> is_query;

			fprintf(stderr, "Query %d begin\n", i); //getchar();	

			CSeekTools::InitVector(is_query, numGenes, (char) 0);
			for(j=0; j<vecstrAllQuery[i].size(); j++){
				q.push_back(mapstriGenes[vecstrAllQuery[i][j]]);
				is_query[mapstriGenes[vecstrAllQuery[i][j]]] = 1;
			}

			vector<float> gene_score;
			vector<char> gene_count;
			CSeekTools::InitVector(gene_score, numGenes, (float) CMeta::GetNaN());
			CSeekTools::InitVector(gene_count, numGenes, (char) 0);
			int qi=0;
			int tot_present_q = 0;
			for(qi=0; qi<q.size(); qi++){
				if(CSeekTools::IsNaN(geneMap->GetForward(q[qi]))) continue;
				tot_present_q++;
				utype qg = geneMap->GetForward(q[qi]);
				for(j=0; j<numActualGenes; j++){
					utype gi = allGenes[j];
					if(correlations[qg][j]==-320){
						continue;
					}
					if(CMeta::IsNaN(gene_score[gi])){
						gene_score[gi] = 0;
					}
					gene_score[gi] += correlations[qg][j];
					gene_count[gi]++;
					//fprintf(stderr, "%.2f\n", correlations[qg][j]);
				}
			}

			int total_non_zero = 0;
			for(j=0; j<numActualGenes; j++){
				utype gi = allGenes[j];
				if(gene_count[gi]!=tot_present_q){
					gene_score[gi] = -320;
					continue;
				}
				if(CMeta::IsNaN(gene_score[gi])){
					gene_score[gi] = -320;
					continue;
				}
				if(is_query[gi]==1){
					gene_score[gi] = -320;
					continue;
				}
				if(tot_present_q==0){
					gene_score[gi] = -320;
					continue;
				}
				gene_score[gi] /= (float) tot_present_q;
				total_non_zero++;
			}

			//fprintf(stderr, "Total non zero genes: %d\n", total_non_zero);

			for(j=0; j<numGenes; j++){
				if(CMeta::IsNaN(gene_score[j])){
					gene_score[j] = -320;
				}
			}

			sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), i);
			CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[i]);
			sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), i);
			CSeekTools::WriteArray(acBuffer, gene_score);
			
			//gs.clear();
			//b.clear();		

        }
	
		
	}
	
	
}
