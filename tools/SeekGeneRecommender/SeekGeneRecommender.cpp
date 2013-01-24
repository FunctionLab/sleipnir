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

bool rank_transform(vector< vector< vector<float> > > &mat, 
	vector<CSeekIntIntMap*> &dm, vector<char> &absent_g, 
	int numGenes){

	int d, dd, i, j, k;
	int numDatasets = mat.size();
	absent_g.clear();
	absent_g.resize(numGenes);
	
	for(i=0; i<numGenes; i++){
		int tot = 0;
		if(i%500==0 || i==numGenes-1) 
			fprintf(stderr, "Gene %d\n", i);
		for(d=0; d<numDatasets; d++){
			CSeekIntIntMap *mapG = dm[d];
			if(CSeekTools::IsNaN(mapG->GetForward(i))) continue;
			tot+=mat[d][0].size(); //mat[d][0].size() is the number of conditions
		}
		if(tot==0){
			absent_g[i] = 1;
			continue;
		}
		
		absent_g[i] = 0;
		vector<bresult> a;
		a.resize(tot);
		
		k=0; 
		for(d=0; d<numDatasets; d++){
			if(CSeekTools::IsNaN(dm[d]->GetForward(i))) continue;
			for(j=0; j<mat[d][0].size(); j++){
				a[k].i = d;
				a[k].j = j;
				a[k].f = mat[d][dm[d]->GetForward(i)][j];
				k++;
			}
		}
		
		sort(a.begin(), a.end(), ascend());
		float b1 = (float) (1.0 - (tot+1.0) / 2.0) / (float) ((tot+1.0)/2.0);
		float b2 = (float) (tot-1.0-(tot+1.0)/2.0) / (float)((tot+1.0)/2.0);
		
		for(j=0; j<tot; j++){
			int id1 = a[j].i;
			int id2 = a[j].j;
			int p = tot+1;
			mat[id1][dm[id1]->GetForward(i)][id2] = 
				((float)(j+1) - (float) p/2.0) / ((float)p/2.0);
		}
		a.clear();	
	}
	return true;	
}

bool weight_experiment(
	vector<vector<vector<float> > > &mat, vector<CSeekIntIntMap*> &dm, 
	int &numExperiments, float cut_off_percentage,
	vector<ushort> &q, vector<bresult> &a){

	int d, dd, j, qi;
	numExperiments = 0;
	int numDatasets = mat.size();
	for(d=0; d<numDatasets; d++){
		//for(j=0; j<mat[d][0].size(); j++){
			int q_size = 0;
			for(qi=0; qi<q.size(); qi++){
				int g = q[qi];
				if(CSeekTools::IsNaN(dm[d]->GetForward(g))) continue;
				q_size++;
			}
			int pq_size = (int)(0.5*q.size());
			if(q_size<pq_size || q_size<2){
				continue;
			}
			numExperiments+=mat[d][0].size();
		//}
	}
	
	a.clear();
	a.resize(numExperiments);
	int numThreads=8;
	omp_set_num_threads(numThreads);
	
	int k=0; 
	/*
	vector<int> array_i, array_j;
	vector<float> array_f, array_f2;
	*/

	vector< vector<int> > array_i, array_j;
	vector< vector<float> > array_f, array_f2;
	array_i.resize(numThreads);
	array_j.resize(numThreads);
	array_f.resize(numThreads);
	array_f2.resize(numThreads);
	
	for(j=0; j<numThreads; j++){
		array_i[j] = vector<int>();
		array_j[j] = vector<int>();
		array_f[j] = vector<float>();
		array_f2[j] = vector<float>();
	}
	
	#pragma omp parallel for \
	shared(numExperiments, mat, q, dm, a, array_i, array_j, array_f, array_f2) \
	private(d, j, qi) \
	firstprivate(numDatasets) schedule(dynamic)
	for(d=0; d<numDatasets; d++){
		ushort tid = omp_get_thread_num();
		for(j=0; j<mat[d][0].size(); j++){
			float mean_query = 0;
			float sample_variance = 0;
			int q_size = 0;
			float significance = 0;
			for(qi=0; qi<q.size(); qi++){
				int g= q[qi];
				if(CSeekTools::IsNaN(dm[d]->GetForward(g))) continue;
				q_size++;
				mean_query+=mat[d][dm[d]->GetForward(g)][j];
			}
			int pq_size = (int)(0.5*q.size());
			if(q_size<pq_size || q_size<2) continue;
			
			mean_query/=(float)q_size;
			for(qi=0; qi<q.size(); qi++){
				int g=q[qi];
				if(CSeekTools::IsNaN(dm[d]->GetForward(g))) continue;
				float diff = mat[d][dm[d]->GetForward(g)][j] - mean_query;
				sample_variance+=diff*diff;
			}
			sample_variance/=(float)q_size;
			significance = sqrt((float)q_size)*mean_query/sqrt(sample_variance + 1.0 /
				(3*numExperiments*numExperiments));
			
			array_i[tid].push_back(d);
			array_j[tid].push_back(j);
			array_f[tid].push_back(significance);
			array_f2[tid].push_back(mean_query);
			
			/*array_i.push_back(d);
			array_j.push_back(j);
			array_f.push_back(significance);
			array_f2.push_back(mean_query);*/
			/*a[k].i = d;
			a[k].j = j;
			a[k].f = sqrt((float)q_size)*mean_query/sqrt(sample_variance + 1.0 /
				(3*numExperiments*numExperiments));
			a[k].f2 = mean_query;
			k++;*/
		}
	}

	/*
	for(d=0; d<array_i.size(); d++){
		a[k].i = array_i[d];
		a[k].j = array_j[d];
		a[k].f = array_f[d];
		a[k].f2 = array_f2[d];
		k++;
	}
	*/
	
	
	for(j=0; j<numThreads; j++){
		for(d=0; d<array_i[j].size(); d++){
			a[k].i = array_i[j][d];
			a[k].j = array_j[j][d];
			a[k].f = array_f[j][d];
			a[k].f2 = array_f2[j][d];
			k++;
		}
	}
	
	//int nth=(int) ((float) numExperiments * cut_off_percentage);
	//std::nth_element(a.begin(), a.begin()+nth, a.end(), descend());
	sort(a.begin(), a.end(), descend());
	return true;
}

bool gene_scoring(
	int numGenes,
	vector<vector<vector<float> > > &mat, vector<CSeekIntIntMap*> &dm, 
	vector<bresult> &a, int numExperiments, vector<char> &absent_g, 
	float cut_off_percentage, vector<float> &gs){

	int i, ii, d, k;
	CSeekTools::InitVector(gs, numGenes, (float) -320);
	vector<int> tot;
	CSeekTools::InitVector(tot, numGenes, (int) 0);
	
	int numDatasets = mat.size();
	CSeekIntIntMap mapG(numGenes);
	
	for(i=0; i<numGenes; i++){
		if(absent_g[i]==1) continue;
		mapG.Add(i);
		for(d=0; d<numDatasets; d++){
			if(CSeekTools::IsNaN(dm[d]->GetForward(i))) continue;
			tot[i]+=mat[d][0].size(); //numConditions in a dataset
		}
	}
	
	int numThreads=8;
	omp_set_num_threads(numThreads);
	
	const vector<ushort> &allGenes = mapG.GetAllReverse();
	float cutoff = cut_off_percentage;
	int numExp = numExperiments;	

	if(numExp==0){
		return true;
	}

	#pragma omp parallel for \
	shared(allGenes, absent_g, a, tot, mat, dm, gs) \
	private(i, ii, k) \
	firstprivate(cutoff, numExp) schedule(dynamic)	
	//fprintf(stderr, "Number of experiments: %d\n", numExp);
	for(ii=0; ii<mapG.GetNumSet(); ii++){
		i = allGenes[ii];
		int num_valid_experiments = 0;
		float sc = 0;
		//fprintf(stderr, "Gene %d, %d of %d\n", i, ii, mapG.GetNumSet());
		for(k=0; k<tot[i]; k++){
			if(k+1 > (int) (cutoff*(float) numExp)){
				break;
			}
			int id1 = a[k].i;
			int id2 = a[k].j;
			
			if(CSeekTools::IsNaN(dm[id1]->GetForward(i))) continue;
			sc+=mat[id1][dm[id1]->GetForward(i)][id2] * a[k].f2;
			num_valid_experiments++;
		}
		//fprintf(stderr, "Valid experiments %d\n", num_valid_experiments);

		int valid_cutoff = (int) (cutoff / 2.0 * (float) numExp);
		if(num_valid_experiments<valid_cutoff){
			continue;
		}
		gs[i] = sc / (float)num_valid_experiments;
	}

    return true;
}

int main(int iArgs, char **aszArgs){
	static const size_t c_iBuffer   = 1024;
	char acBuffer[c_iBuffer];
	
	gengetopt_args_info sArgs;
	ifstream ifsm;
	istream *pistm;
	vector<string> vecstrLine, vecstrGenes, vecstrDBs, vecstrQuery;
	ushort i, qi, j, k, l;
	
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

	if(sArgs.pcl_flag==1){
		string pcl_dir = sArgs.pcl_dir_arg;
		string output_dir = sArgs.dir_out_arg;
		vector<string> pcl_list;
		vector< vector<string> > vecstrAllQuery;
		int numGenes = vecstrGenes.size();
		
		if(!CSeekTools::ReadMultipleQueries(sArgs.query_arg, vecstrAllQuery))
			return -1;
		
		CSeekTools::ReadListOneColumn(sArgs.pcl_list_arg, pcl_list);
		vector< vector< vector<float> > > mat;
		vector<CSeekIntIntMap*> dm;
		dm.resize(pcl_list.size());
		mat.resize(pcl_list.size());
		
		for(i=0; i<pcl_list.size(); i++){
			fprintf(stderr, "Reading %d: %s\n", i, pcl_list[i].c_str());
			dm[i] = new CSeekIntIntMap(vecstrGenes.size());
			
			string pclfile = pcl_dir + "/" + pcl_list[i] + ".bin";

			CPCL pcl;
			pcl.Open(pclfile.c_str());
			int totNumExperiments = pcl.GetExperiments() - 2;
			
			vector<ushort> presentIndex;
			vector<string> presentGeneNames;
			for(j=0; j<vecstrGenes.size(); j++){
				ushort g = pcl.GetGene(vecstrGenes[j]);
				if(CSeekTools::IsNaN(g)) continue; //gene does not exist in the dataset
				presentIndex.push_back(g);
				presentGeneNames.push_back(vecstrGenes[j]);
				dm[i]->Add(j);
			}
			
			mat[i].resize(presentIndex.size());
			for(j=0; j<presentIndex.size(); j++)
				mat[i][j].resize(totNumExperiments);
			
			for(j=0; j<presentIndex.size(); j++){
				float *val = pcl.Get(presentIndex[j]);
				for(k=2; k<pcl.GetExperiments(); k++)
					mat[i][j][k-2] = val[k];					
			}

		}
		
		fprintf(stderr, "Finished loading...\n"); //getchar();
		vector<char> absent_g;
		rank_transform(mat, dm, absent_g, vecstrGenes.size());
		fprintf(stderr, "Finished transforming...\n"); //getchar();
		
		for(i=0; i<vecstrAllQuery.size(); i++){
			vector<ushort> q;
			vector<char> is_query;

			fprintf(stderr, "Query %d begin\n", i); //getchar();	

			CSeekTools::InitVector(is_query, numGenes, (char) 0);
			for(j=0; j<vecstrAllQuery[i].size(); j++){
				q.push_back(mapstriGenes[vecstrAllQuery[i][j]]);
				is_query[mapstriGenes[vecstrAllQuery[i][j]]] = 1;
			}
			int numExperiments = 0;
			vector<bresult> b;
			vector<float> gs;
			fprintf(stderr, "Query %d finished initializing query\n", i); //getchar();
			weight_experiment(mat, dm, numExperiments, 0.05, q, b);

			fprintf(stderr, "Query %d finished weighting\n", i); //getchar();
			gene_scoring(vecstrGenes.size(), mat, dm, b, numExperiments, absent_g, 
				0.05, gs);
	
			fprintf(stderr, "Query %d finished scoring\n", i); //getchar();
			
			for(j=0; j<absent_g.size(); j++){
				if(absent_g[j]==1){
					gs[j] = -320;
				}
			}

			
			fprintf(stderr, "Query %d end\n", i); //getchar();

			sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), i);
			CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[i]);
			sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), i);
			CSeekTools::WriteArray(acBuffer, gs);

			/*
			fprintf(stdout, "Query %d\n", i);
			fprintf(stdout, "Weights:");
			//fprintf(stderr, "Weights:\n");
			for(j=0; j<b.size() && j<(int)(0.05*numExperiments); j++){
				if(j%10==0){
					fprintf(stdout, "\n");
				}	
				fprintf(stdout, "%.2f ", b[j].f);
				//fprintf(stderr, "%d %.2f\n", j, b[j].f);
			}
			fprintf(stdout, "\n");
			vector<AResultFloat> gsc;
			gsc.resize(numGenes);
			//fprintf(stderr, "numexp %d gsc %d gs %d\n", numExperiments, gsc.size(), gs.size());
			for(j=0; j<numGenes; j++){
				gsc[j].i = j;
				gsc[j].f = gs[j];
				//fprintf(stderr, "Gene %d %.2f\n", j, gsc[j].f);
			}
			sort(gsc.begin(), gsc.end());
        	fprintf(stdout, "Results:\n");
			int kk=0;
			for(k=0; k<numGenes && kk<500; k++){
				if(is_query[gsc[k].i]==1) continue;
				if(gsc[k].f==-320) continue;
				fprintf(stdout, "%s %.5f\n", vecstrGenes[gsc[k].i].c_str(), gsc[k].f);
				kk++;
			}
			*/
			/*
			fprintf(stderr, "Query %d\n", i);
			fprintf(stderr, "Weights:");
			for(j=0; j<(int)(0.05* numExperiments); j++){
				if(j%10==0){
					fprintf(stderr, "\n");
				}
				fprintf(stderr, "%.2f ", b[j].f);
			}
			fprintf(stderr, "\n");
			vector<AResultFloat> gsc;
			gsc.resize(numGenes);
			for(j=0; j<numGenes; j++){
				gsc[j].i = j;
				gsc[j].f = gs[j];
			}
			sort(gsc.begin(), gsc.end());
			fprintf(stderr, "Results:\n");
			int kk=0;
			for(k=0; k<numGenes && kk<500; k++){
				if(is_query[gsc[k].i]==1) continue;
				if(gsc[k].f==-320) continue;
				fprintf(stderr, "%d %.5f\n", gsc[k].i, gsc[k].f);
				kk++;
			}
			*/
			gs.clear();
			b.clear();		

        }
				
		
	}
	
	
}
