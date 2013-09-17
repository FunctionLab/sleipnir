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

enum QUERY_MODE{
	SINGLE_QUERY, MULTI_QUERY, MULTI_DWEIGHT
};

enum METRIC{
	RBP, AVGP, PR, PR_ALL, AUC
};

bool GetRandom(gsl_rng *r, const vector<AResultFloat> &geneScore, 
	vector<AResultFloat> &random, vector<char> &excludeGene,
	vector<char> &includeGene, vector<utype> &queryGeneID, 
	const float nan){

	int i, j;
	int ss = 0;

	vector<char> q;
	q.clear();
	CSeekTools::InitVector(q, geneScore.size(), (char) 0);
	for(j=0; j<queryGeneID.size(); j++)
		q[queryGeneID[j]] = 1;

	for(i=0; i<geneScore.size(); i++){
		if(includeGene[geneScore[i].i]==0){
		}else if(q[geneScore[i].i]==1){
		}else if(excludeGene[geneScore[i].i]==1){
		}else{
			ss++;
		}
	}

	float *gs = (float*)malloc(ss*sizeof(float));
	random.clear();
	random.resize(geneScore.size());
	
	int ii = 0;
	for(i=0; i<geneScore.size(); i++){
		if(includeGene[geneScore[i].i]==0){
		}else if(q[geneScore[i].i]==1){
		}else if(excludeGene[geneScore[i].i]==1){
		}else{
			//gs[ii] = geneScore[i].f;
			gs[ii] = (float) ii;
			ii++;
		}
	}

	gsl_ran_shuffle(r, gs, ss, sizeof(float));

	ii = 0;
	for(i=0; i<geneScore.size(); i++){
		random[i].i = geneScore[i].i;
		random[i].f = nan;
		if(includeGene[geneScore[i].i]==0){
		}else if(q[geneScore[i].i]==1){
		}else if(excludeGene[geneScore[i].i]==1){
		}else{
			random[i].f = gs[ii];
			ii++;
		}
	}
	free(gs);
	return true;
}

float RankBiasedPrecision(const float &rate,
	const vector<AResultFloat> &sortedScore, const vector<char> &gold_std,
	const float &nanVal){

	vector<AResultFloat>::const_iterator iterScore = sortedScore.begin();
	float x = 0;
	utype i;
	for(i=0; iterScore!=sortedScore.end(); i++, iterScore++){
		//if(iterScore->f == nanVal) continue;
		if(gold_std[iterScore->i]==1) x += pow(rate, i);
	}
	x*=(1.0-rate);
	return x;
}

vector<float>* Precision(const vector<AResultFloat> &sortedScore,
	const vector<char> &gold_std, const float &nanVal){

	utype i, numPos = 1;
	vector<float> *r = new vector<float>();

	vector<AResultFloat>::const_iterator iterScore = sortedScore.begin();
	for(i=0; iterScore!=sortedScore.end(); i++, iterScore++){
		//if(iterScore->f == nanVal) continue;
		if(gold_std[iterScore->i]==1){
			//fprintf(stderr, "%d %d\n", numPos, i);
			r->push_back((float) numPos++ / (float) (i + 1));
		}
	}
	r->resize(r->size());
	return r;
}

bool MeanStandardDeviation(const vector<float> &f, float &avg, float &stdev){
	avg = 0;
	stdev = 0;
	avg = std::accumulate(f.begin(), f.end(), (float) 0);
	avg /= (float) f.size();
	vector<float>::const_iterator iterF = f.begin();
	for(; iterF!=f.end(); iterF++) stdev += (*iterF - avg) * (*iterF - avg);
	stdev /= (float) f.size();
	stdev = sqrt(stdev);
	return true;
}

bool MinMaxQuartile(vector<float> &f, float &min, float &max, float &Q1,
		float &Q2, float &Q3){
	min = *min_element(f.begin(), f.end());
	max = *max_element(f.begin(), f.end());
	Q1 = 0;
	Q2 = 0;
	Q3 = 0;

	if(f.size()==1){
		Q1 = f[0];
		Q2 = Q1;
		Q3 = Q2;
		return true;
	}
	if(f.size()==2){
		Q1 = f[0];
		Q2 = (f[0] + f[1] ) / 2.0;
		Q3 = f[1];
		return true;
	}
	if(f.size()==3){
		Q1 = (f[0] + f[1] ) / 2.0;
		Q2 = f[1];
		Q3 = (f[1] + f[2]) / 2.0;
		return true;
	}

	if(f.size()%2==0){
		int m1 = f.size()/2 - 1;
		int m2 = f.size()/2;
		std::nth_element(f.begin(), f.begin()+m1, f.end());
		float f1 = f[m1];
		std::nth_element(f.begin(), f.begin()+m2, f.end());
		float f2 = f[m2];
		Q2 = (f1 + f2) / 2.0;

		int s1 = m1 + 1;
		if(s1%2==0){
			int m3 = s1/2 - 1;
			int m4 = s1/2;
			std::nth_element(f.begin(), f.begin()+m3, f.end());
			float f3 = f[m3];
			std::nth_element(f.begin(), f.begin()+m4, f.end());
			float f4 = f[m4];
			Q1 = (f3 + f4) / 2.0;
		}else{
			int m3 = s1/2;
			std::nth_element(f.begin(), f.begin()+m3, f.end());
			float f3 = f[m3];
			Q1 = f3;
		}

		int s2 = (f.size()-1) - m2 + 1;
		if(s2%2==0){
			int m3 = m2 + s2/2 - 1;
			int m4 = m2 + s2/2;
			std::nth_element(f.begin(), f.begin()+m3, f.end());
			float f3 = f[m3];
			std::nth_element(f.begin(), f.begin()+m4, f.end());
			float f4 = f[m4];
			Q3 = (f3 + f4) / 2.0;
		}else{
			int m3 = m2 + s2/2;
			std::nth_element(f.begin(), f.begin()+m3, f.end());
			float f3 = f[m3];
			Q3 = f3;
		}

	}else{
		int m1 = f.size()/2;
		std::nth_element(f.begin(), f.begin()+m1, f.end());
		float f1 = f[m1];
		Q2 = f1;

		int s1 = m1 - 1 + 1;
		if(s1%2==0){
			int m3 = s1/2 - 1;
			int m4 = s1/2;
			std::nth_element(f.begin(), f.begin()+m3, f.end());
			float f3 = f[m3];
			std::nth_element(f.begin(), f.begin()+m4, f.end());
			float f4 = f[m4];
			Q1 = (f3 + f4) / 2.0;
		}else{
			int m3 = s1/2;
			std::nth_element(f.begin(), f.begin()+m3, f.end());
			float f3 = f[m3];
			Q1 = f3;
		}

		int s2 = (f.size()-1) - (m1 + 1) + 1;
		if(s2%2==0){
			int m3 = m1 + 1 + s2/2 - 1;
			int m4 = m1 + 1 + s2/2;
			std::nth_element(f.begin(), f.begin()+m3, f.end());
			float f3 = f[m3];
			std::nth_element(f.begin(), f.begin()+m4, f.end());
			float f4 = f[m4];
			Q3 = (f3 + f4) / 2.0;
		}else{
			int m3 = m1 + 1 + s2/2;
			std::nth_element(f.begin(), f.begin()+m3, f.end());
			float f3 = f[m3];
			Q3 = f3;
		}
	}
	return true;
}

bool EvaluateOneQuery(const gengetopt_args_info &sArgs, const enum METRIC &met,
	const vector<AResultFloat> &sortedGenes,
	const vector<char> &goldstdGenePresence, const float &nan,
	vector<float> &eval){
	size_t i;
	eval.clear();
	//metric
	if(met==PR_ALL){
		vector<float> *vf = Precision(sortedGenes, goldstdGenePresence, nan);
		for(i=0; i<vf->size(); i++) eval.push_back(vf->at(i));
		delete vf;
		return true;
	}
	fprintf(stderr, "Invalid option!\n");
	return false;
}


bool EvaluateOneQuery(const gengetopt_args_info &sArgs, const enum METRIC &met,
	const vector<AResultFloat> &sortedGenes,
	const vector<char> &goldstdGenePresence, const float &nan, float &eval){
	size_t i;
	//metric
	if(met==RBP){
		float rate = sArgs.rbp_p_arg;
		float rbp = RankBiasedPrecision(rate, sortedGenes,
			goldstdGenePresence, nan);
		eval = rbp;
		//fprintf(stdout, "%.5f\n", rbp);
		return true;
	}
	else if(met==AVGP || met==PR){
		vector<float> *vf = Precision(sortedGenes, goldstdGenePresence, nan);
		/*if(met==PR_ALL){
			for(i=0; i<vf->size()-1; i++){
				fprintf(stdout, "%.5f ", vf->at(i));
				}
				fprintf(stdout, "%.5f\n", vf->at(i));
				delete vf;
				return true;
			}*/

		if(sArgs.x_int_arg==-1 && sArgs.x_per_arg==0){
			fprintf(stderr, "Must specify --x_int or --x_per\n");
			delete vf;
			return false;
		}

		int X = -1;
		if(sArgs.x_int_arg!=-1 && sArgs.x_per_arg==0){
			X = sArgs.x_int_arg;
			if(X>vf->size()){
				fprintf(stderr, "Error: X is too large (>%d)\n", vf->size());
				delete vf;
				return false;
			}
		}
		else if(sArgs.x_int_arg==-1 && sArgs.x_per_arg>0){
			float per = sArgs.x_per_arg;
			if(per>1.0){
				fprintf(stderr, "Error: percentage is above 1.0\n");
				delete vf;
				return false;
			}
			X = (int) ( (float) per * (float) vf->size() );
			//fprintf(stderr, "%.5f %d %d", per, vf->size(), X);
			X = std::max(1, X);
		}

		if(met==AVGP){
			float avg = std::accumulate(vf->begin(), vf->begin()+X, (float) 0);
			avg /= (float) X;
			//fprintf(stdout, "%.5f\n", avg);
			eval = avg;
			delete vf;
			return true;
		}
		if(met==PR){
			//fprintf(stdout, "%.5f\n", vf->at(X-1));
			eval = vf->at(X-1);
			delete vf;
			return true;
		}
	}

	fprintf(stderr, "Invalid option!\n");
	return false;
}

void PrintVector(const vector<float> &f){
	vector<float>::const_iterator iterF = f.begin();
	vector<float>::const_iterator end = f.begin() + f.size() - 1;
	for(; iterF!=end; iterF++) fprintf(stderr, "%.5f ", *iterF);
	fprintf(stderr, "%.5f\n", *iterF);
}

void PrintResult(vector< vector<float> > f){
	int i;
	for(i=0; i<f.size(); i++){
		PrintVector(f[i]);
	}
}

bool DoAggregate(const gengetopt_args_info &sArgs, const enum METRIC &met, 
	vector<AResultFloat> *sortedGenes, vector<utype> *queryGeneID, 
	int listSize, vector<char> *goldstdGenePresence, vector<char> *excludeGene,
	vector<char> *includeGene, 
	vector< vector<float> > &result

	){

	//fprintf(stderr, "Finished reading gene score list\n");
	vector<AResultFloat> *master = NULL;
	vector<AResultFloat> master_score;
	vector<AResultFloat> master_rank;
	utype i, j;
	float nan = sArgs.nan_arg;
	result.clear();

	vector<char> *q = new vector<char>[listSize];
	for(i=0; i<listSize; i++){
		q[i] = vector<char>();
		CSeekTools::InitVector(q[i], sortedGenes[i].size(), (char) 0);
		for(j=0; j<queryGeneID[i].size(); j++)
			q[i][queryGeneID[i][j]] = 1;
	}

	if(sArgs.agg_ranksum_flag==1 || sArgs.agg_scoresum_flag==1){
		//Gold standard must be the same across all queries for this mode!!
		//ASSUME THIS IS TRUE
		for(i=0; i<listSize; i++){
			for(j=0; j<sortedGenes[0].size(); j++){
				if(q[i][sortedGenes[i][j].i]==1){
					sortedGenes[i][j].f = nan;
				}
				if(excludeGene[i][sortedGenes[i][j].i]==1){
					sortedGenes[i][j].f = nan;
				}
			}
		}

		if(sArgs.agg_scoresum_flag==1){
			master_score.resize(sortedGenes[0].size());

			for(j=0; j<sortedGenes[0].size(); j++){
				master_score[j].i = j;
				master_score[j].f = 0.0;
			}

			for(j=0; j<sortedGenes[0].size(); j++)
				for(i=0; i<listSize; i++)
					master_score[sortedGenes[i][j].i].f += sortedGenes[i][j].f;

			for(j=0; j<sortedGenes[0].size(); j++)
				master_score[j].f /= (float)listSize;

			sort(master_score.begin(), master_score.end());

			master = &master_score;
		}

		else if(sArgs.agg_ranksum_flag==1){
			//fprintf(stderr, "Got here 2\n"); 
			master_rank.resize(sortedGenes[0].size());

			//fprintf(stderr, "Got here 2a\n"); 

			for(i=0; i<listSize; i++){
				sort(sortedGenes[i].begin(), sortedGenes[i].end());
				//fprintf(stderr, "Sort %d fine\n", i); 
			}	

			//fprintf(stderr, "Got here 2b\n"); 

			for(j=0; j<sortedGenes[0].size(); j++){
				master_rank[j].i = j;
				master_rank[j].f = 0;
			}

			for(i=0; i<listSize; i++)
				for(j=0; j<sortedGenes[i].size(); j++)
					master_rank[sortedGenes[i][j].i].f +=
						(float) sortedGenes[i].size() - (float) j;

			for(j=0; j<sortedGenes[0].size(); j++)
				master_rank[j].f /= (float) listSize;

			sort(master_rank.begin(), master_rank.end());

			master = &master_rank;

			//fprintf(stderr, "Got here 3\n"); 
		}

		if(met!=PR_ALL){
			float eval;
			bool ret = EvaluateOneQuery(sArgs, met, *master,
				goldstdGenePresence[0], CMeta::GetNaN(), eval);
			if(!ret) return 1;
			result.resize(1);
			result[0] = vector<float>();
			result[0].push_back(eval);
			return 0;
		}else{
			vector<float> evalAll;
			bool ret = EvaluateOneQuery(sArgs, met, *master,
				goldstdGenePresence[0], CMeta::GetNaN(), evalAll);
			if(!ret) return 1;
			result.resize(1);
			result[0] = vector<float>();
			for(i=0; i<evalAll.size(); i++) 
				result[0].push_back(evalAll[i]);
			return 0;
		}
	}

	//Not Rank Sum
	vector<float> eval;
	vector< vector<float> > vecevalAll;

	eval.resize(listSize);
	vecevalAll.resize(listSize);

	for(i=0; i<listSize; i++){
		for(j=0; j<sortedGenes[i].size(); j++){
			//NEW
			if(includeGene[i][sortedGenes[i][j].i]==0){
				sortedGenes[i][j].f = nan;
			}
			if(q[i][sortedGenes[i][j].i]==1){
				sortedGenes[i][j].f = nan;
			}
			if(excludeGene[i][sortedGenes[i][j].i]==1){
				sortedGenes[i][j].f = nan;
			}
		}

		sort(sortedGenes[i].begin(), sortedGenes[i].end());

		if(met!=PR_ALL){
			float fEval;
			bool ret = EvaluateOneQuery(sArgs, met, sortedGenes[i],
				goldstdGenePresence[i], nan, fEval);
			if(!ret) return 1;
			eval[i] = fEval;

		}else{
			vector<float> evalAll;
			bool ret = EvaluateOneQuery(sArgs, met, sortedGenes[i],
				goldstdGenePresence[i], nan, evalAll);
			if(!ret) return 1;
			vecevalAll[i] = evalAll;
		}
	}

	if(sArgs.fold_over_random_flag==1){
		const gsl_rng_type *T;
		gsl_rng *rnd;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rnd = gsl_rng_alloc(T);

		vector<AResultFloat> *randomScores = 
			new vector<AResultFloat>[listSize];

		vector<float> random_eval;
		vector< vector<float> > random_vecevalAll;
		random_eval.resize(listSize);
		random_vecevalAll.resize(listSize);

		//should shuffle only within annotated genes
		for(i=0; i<listSize; i++){
			GetRandom(rnd, sortedGenes[i], randomScores[i], excludeGene[i], 
				includeGene[i], queryGeneID[i], nan);
			sort(randomScores[i].begin(), randomScores[i].end());
			if(met!=PR_ALL){
				float fEval;
				bool ret = EvaluateOneQuery(sArgs, met, randomScores[i],
					goldstdGenePresence[i], nan, fEval);
				if(!ret) return 1;
				random_eval[i] = fEval;
				//calculate fold
				if(random_eval[i]==eval[i]) eval[i] = 1.0;
				else if(random_eval[i]==0) eval[i] = 1.0;
				else eval[i] = eval[i] / random_eval[i];
			}else{
				vector<float> evalAll;
				bool ret = EvaluateOneQuery(sArgs, met, randomScores[i],
					goldstdGenePresence[i], nan, evalAll);
				if(!ret) return 1;
				random_vecevalAll[i] = evalAll;
				int ii=0; 
				//calculate fold
				for(ii=0; ii<random_vecevalAll[i].size(); ii++){
					if(random_vecevalAll[i][ii]==vecevalAll[i][ii])
						vecevalAll[i][ii] = 1.0;
					else if(random_vecevalAll[i][ii]==0)
						vecevalAll[i][ii] = 1.0;
					else
						vecevalAll[i][ii] /= random_vecevalAll[i][ii];
				}
			}
			//fprintf(stderr, "Got here!\n");	
		}
		//fprintf(stderr, "Got here\n");
		gsl_rng_free(rnd);
		delete[] randomScores;		
		//fprintf(stderr, "Got here 2\n");
	}

	if(met!=PR_ALL){
		if(sArgs.agg_avg_flag==1){
			float avg, stdev;
			MeanStandardDeviation(eval, avg, stdev);
			result.resize(1);
			result[0] = vector<float>();
			result[0].push_back(avg);
			result[0].push_back(stdev);
			return 0;
		}
		if(sArgs.agg_quartile_flag==1){
			float min, max, Q1, Q2, Q3;
			MinMaxQuartile(eval, min, max, Q1, Q2, Q3);
			result.resize(2);
			result[0] = vector<float>();
			result[0].push_back(min);
			result[0].push_back(max);
			result[1] = vector<float>();
			result[1].push_back(Q1);
			result[1].push_back(Q2);
			result[1].push_back(Q3);
			return 0;
		}
		if(sArgs.display_all_flag==1){
			result.resize(eval.size());
			for(i=0; i<eval.size(); i++){
				result[i] = vector<float>();
				result[i].push_back(eval[i]);
			}
			return 0;
		}

	}else{
		vector< vector<float> > veceval;
		veceval.resize(vecevalAll[0].size());
		for(i=0; i<vecevalAll.size(); i++)
			for(j=0; j<vecevalAll[i].size(); j++)
				veceval[j].push_back(vecevalAll[i][j]);

		if(sArgs.agg_avg_flag==1){
			vector<float> avgAll, stdevAll;
			for(i=0; i<veceval.size(); i++){
				float avg, stdev;
				MeanStandardDeviation(veceval[i], avg, stdev);
				avgAll.push_back(avg);
				stdevAll.push_back(stdev);
			}
			result.resize(2);
			result[0] = vector<float>();
			result[1] = vector<float>();
			for(i=0; i<avgAll.size(); i++){
				result[0].push_back(avgAll[i]);
			}
			for(i=0; i<stdevAll.size(); i++){
				result[1].push_back(stdevAll[i]);
			}
			return 0;
		}
		else if(sArgs.agg_quartile_flag==1){
			vector<float> minAll, maxAll, Q1All, Q2All, Q3All;
			for(i=0; i<veceval.size(); i++){
				float min, max, Q1, Q2, Q3;
				MinMaxQuartile(veceval[i], min, max, Q1, Q2, Q3);
				minAll.push_back(min);
				maxAll.push_back(max);
				Q1All.push_back(Q1);
				Q2All.push_back(Q2);
				Q3All.push_back(Q3);
			}
			result.resize(5);
			result[0] = vector<float>();
			result[1] = vector<float>();
			result[2] = vector<float>();
			result[3] = vector<float>();
			result[4] = vector<float>();
			for(i=0; i<minAll.size(); i++){
				result[0].push_back(minAll[i]);
			}
			for(i=0; i<maxAll.size(); i++){
				result[1].push_back(maxAll[i]);
			}
			for(i=0; i<Q1All.size(); i++){
				result[2].push_back(Q1All[i]);
			}
			for(i=0; i<Q2All.size(); i++){
				result[3].push_back(Q2All[i]);
			}
			for(i=0; i<Q3All.size(); i++){
				result[4].push_back(Q3All[i]);
			}
			return 0;
		}
	}

	return true;
}

int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuffer	= 1024;
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	istream*			pistm;
	vector<string>		vecstrLine, vecstrGenes, vecstrDBs, vecstrQuery;
	char				acBuffer[ c_iBuffer ];
	size_t				i, j, k;


	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	/* reading gene-mapping file */
	ifsm.open(sArgs.input_arg);
	pistm = &ifsm;

	map<string, size_t> mapstriGenes;
	while( !pistm->eof( ) ) {
		pistm->getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) < 2 ) {
			//cerr << "Ignoring line: " << acBuffer << endl;
			continue;
		}
		if( !( i = atoi( vecstrLine[ 0 ].c_str( ) ) ) ) {
			cerr << "Illegal gene ID: " << vecstrLine[ 0 ] << " for "
				<< vecstrLine[ 1 ] << endl;
			return 1;
		}
		i--;
		if( vecstrGenes.size( ) <= i ) vecstrGenes.resize( i + 1 );
		vecstrGenes[ i ] = vecstrLine[ 1 ];
		mapstriGenes[vecstrGenes[i]] = i;
	}
	ifsm.close( );

	enum QUERY_MODE qmode;
	if(sArgs.single_flag==1) qmode = SINGLE_QUERY;
	else if(sArgs.aggregate_flag==1) qmode = MULTI_QUERY;
	else if(sArgs.multi_weight_flag==1) qmode = MULTI_DWEIGHT;

	enum METRIC met;
	if(sArgs.rbp_flag==1) met = RBP;
	else if(sArgs.avgp_flag==1) met = AVGP;
	else if(sArgs.pr_flag==1) met = PR;
	else if(sArgs.pr_all_flag==1) met = PR_ALL;
	else if(sArgs.auc_flag==1) met = AUC;

	//fprintf(stderr, "Query mode: %d\n", qmode);

	if(qmode==SINGLE_QUERY){
		if(sArgs.display_weight_flag==1){
			vector<float> ww;
			CSeekTools::ReadArray(sArgs.weight_arg, ww);
			vector<string> vecstrDatasets, vecstrP;
			CSeekTools::ReadListTwoColumns(sArgs.dataset_map_arg, 
				vecstrDatasets, vecstrP);
			vector<AResultFloat> sortedDatasets;
			sortedDatasets.resize(ww.size());
			for(i=0; i<sortedDatasets.size(); i++){
				sortedDatasets[i].i = i;
				sortedDatasets[i].f = ww[i];
			}
			sort(sortedDatasets.begin(), sortedDatasets.end());
			for(i=0; i<sortedDatasets.size(); i++){
				fprintf(stderr, "%.2e\t%s\n", sortedDatasets[i].f, 
					vecstrDatasets[sortedDatasets[i].i].c_str());
			}
			return 0;
		}

		string queryFile = sArgs.query_arg;
		vector<string> queryGenes;
		CSeekTools::ReadMultiGeneOneLine(queryFile, queryGenes);
		vector<utype> queryGeneID;
		for(i=0; i<queryGenes.size(); i++)
			queryGeneID.push_back(mapstriGenes[queryGenes[i]]);

		string genescoreFile = sArgs.gscore_arg;
		vector<float> geneScores;
		CSeekTools::ReadArray(genescoreFile.c_str(), geneScores);
		//float maxScore = *std::max_element(geneScores.begin(),
		//	geneScores.end());

		float nan = sArgs.nan_arg;
		vector<AResultFloat> sortedGenes;
		sortedGenes.resize(geneScores.size());
		for(i=0; i<sortedGenes.size(); i++){
			sortedGenes[i].i = i;
			sortedGenes[i].f = geneScores[i];
		}

		//Query genes themselves have lowest score, to prevent
		//them from being counted in PR
		for(i=0; i<queryGeneID.size(); i++)
			sortedGenes[queryGeneID[i]].f = nan;

		sort(sortedGenes.begin(), sortedGenes.end());

		if(sArgs.p_value_flag==1){
			string random_directory = sArgs.random_dir_arg;
			int num_random = sArgs.random_num_arg;
			int ii, jj;
			char ac[256];
			vector<vector<int> > randomRank;
			vector<vector<float> > randomSc;
			vector<int> geneRank;

			randomRank.resize(sortedGenes.size());
			randomSc.resize(sortedGenes.size());
			geneRank.resize(sortedGenes.size());
			for(ii=0; ii<sortedGenes.size(); ii++){
				randomRank[ii].resize(num_random);
				randomSc[ii].resize(num_random);
			}

			for(ii=0; ii<num_random; ii++){
				vector<float> randomScores;
				sprintf(ac, "%s/%d.gscore", random_directory.c_str(), ii);
				CSeekTools::ReadArray(ac, randomScores);
				vector<AResultFloat> sortedRandom;
				sortedRandom.resize(randomScores.size());
				for(jj=0; jj<randomScores.size(); jj++){
					sortedRandom[jj].i = jj;
					sortedRandom[jj].f = randomScores[jj];
				}
				sort(sortedRandom.begin(), sortedRandom.end());
				for(jj=0; jj<randomScores.size(); jj++){
					randomRank[sortedRandom[jj].i][ii] = jj;
					randomSc[sortedRandom[jj].i][ii] = sortedRandom[jj].f;
				}
			}

			for(jj=0; jj<geneScores.size(); jj++){
				sort(randomRank[jj].begin(), randomRank[jj].end());
				sort(randomSc[jj].begin(), randomSc[jj].end(), std::greater<float>());
				geneRank[sortedGenes[jj].i] = jj;
			}

			for(jj=0; jj<geneScores.size(); jj++){
				int gene = sortedGenes[jj].i;
				int gene_rank = jj;
				float gene_score = sortedGenes[jj].f;
				if(gene_score<0) 
					continue;
				vector<int> &rR = randomRank[gene];
				vector<float> &rF = randomSc[gene];
				int kk = 0;
				for(kk=0; kk<rR.size(); kk++){
					//if(gene_rank<=rR[kk]){
					if(gene_score>=rF[kk]){
						//float f1 = (float) kk / (float) rR.size();
						fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
							gene_rank, kk, gene_score, randomSc[gene][kk]);
						break;
					}else if(kk==rR.size()-1){
						//float f1 = (float) kk / (float) rR.size();
						fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
							gene_rank, kk, gene_score, randomSc[gene][kk]);
						//fprintf(stderr, "%d %.6f\n", gene_rank, f1);
						//fprintf(stderr, "%d %d / 100\n", gene_rank, kk);
					}
				}
			}

			return 0;
		}

		if(sArgs.dislay_only_flag==1){
			for(i=0; i<15000; i++)
				fprintf(stderr, "%s\t%.5f\n", 
					vecstrGenes[sortedGenes[i].i].c_str(), sortedGenes[i].f);
			return 0;
		}

		string goldstdFile = sArgs.goldstd_arg;
		vector<string> goldstdGenes;
		CSeekTools::ReadMultiGeneOneLine(goldstdFile, goldstdGenes);
		vector<char> goldstdGenePresence;
		CSeekTools::InitVector(goldstdGenePresence,
			vecstrGenes.size(), (char) 0);

		for(i=0; i<goldstdGenes.size(); i++)
			goldstdGenePresence[mapstriGenes[goldstdGenes[i]]] = 1;

		if(met!=PR_ALL){
			float eval;
			bool ret = EvaluateOneQuery(sArgs, met, sortedGenes,
				goldstdGenePresence, nan, eval);
			if(!ret) return 1;
			fprintf(stderr, "%.5f\n", eval);
			return 0;
		}else{
			vector<float> evalAll;
			//fprintf(stderr, "Got here");
			bool ret = EvaluateOneQuery(sArgs, met, sortedGenes,
				goldstdGenePresence, nan, evalAll);
			if(!ret) return 1;
			PrintVector(evalAll);
			return 0;
		}
	}

	if(qmode == MULTI_QUERY){
		string goldstdList = sArgs.goldstd_list_arg;
		vector<string> vecstrList;
		CSeekTools::ReadListOneColumn(goldstdList, vecstrList);
		vector<char> *goldstdGenePresence =
			new vector<char>[vecstrList.size()];
		for(i=0; i<vecstrList.size(); i++){
			vector<string> goldstdGenes;
			CSeekTools::ReadMultiGeneOneLine(vecstrList[i], goldstdGenes);
			CSeekTools::InitVector(goldstdGenePresence[i],
				vecstrGenes.size(), (char) 0);
			for(j=0; j<goldstdGenes.size(); j++)
				goldstdGenePresence[i][mapstriGenes[goldstdGenes[j]]] = 1;
		}

		//fprintf(stderr, "Finished reading gold standard list\n");

		string queryList = sArgs.query_list_arg;
		vecstrList.clear();
		CSeekTools::ReadListOneColumn(queryList, vecstrList);
		vector<utype> *queryGeneID = new vector<utype>[vecstrList.size()];
		for(i=0; i<vecstrList.size(); i++){
			vector<string> queryGenes;
			CSeekTools::ReadMultiGeneOneLine(vecstrList[i], queryGenes);
			for(j=0; j<queryGenes.size(); j++)
				queryGeneID[i].push_back(mapstriGenes[queryGenes[j]]);
		}

		//fprintf(stderr, "Finished reading query list\n");

		vector<string> exclude_list;
		string excl = sArgs.exclude_list_arg;
		exclude_list.clear();
		CSeekTools::ReadListOneColumn(excl, exclude_list);
		vector<char> *excludeGene = new vector<char>[vecstrList.size()];
		for(i=0; i<exclude_list.size(); i++){
			vector<string> ex;
			CSeekTools::ReadMultiGeneOneLine(exclude_list[i], ex);
			CSeekTools::InitVector(excludeGene[i], vecstrGenes.size(), (char) 0);
			for(j=0; j<ex.size(); j++)
				excludeGene[i][mapstriGenes[ex[j]]] = 1;
		}

		vector<string> include_list;
		string incl = sArgs.include_list_arg;
		include_list.clear();
		CSeekTools::ReadListOneColumn(incl, include_list);
		vector<char> *includeGene = new vector<char>[vecstrList.size()];
		for(i=0; i<include_list.size(); i++){
			vector<string> in;
			CSeekTools::ReadMultiGeneOneLine(include_list[i], in, 40000);
			CSeekTools::InitVector(includeGene[i], vecstrGenes.size(), (char) 0);
			for(j=0; j<in.size(); j++)
				includeGene[i][mapstriGenes[in[j]]] = 1;
		}


		string genescoreList = sArgs.gscore_list_arg;
		vecstrList.clear();
		CSeekTools::ReadListOneColumn(genescoreList, vecstrList);
		vector<float> *geneScores = new vector<float>[vecstrList.size()];
		vector<AResultFloat> *sortedGenes =
			new vector<AResultFloat>[vecstrList.size()];
		float nan = sArgs.nan_arg;
		for(i=0; i<vecstrList.size(); i++){
			CSeekTools::ReadArray(vecstrList[i].c_str(), geneScores[i]);
			//maxScore[i] = *std::max_element(geneScores[i].begin(),
			//	geneScores[i].end());
			sortedGenes[i].resize(geneScores[i].size());
			for(j=0; j<sortedGenes[i].size(); j++){
				sortedGenes[i][j].i = j;
				sortedGenes[i][j].f = geneScores[i][j];
				if(isnan(geneScores[i][j]) || isinf(geneScores[i][j])){
					sortedGenes[i][j].f = nan;
				}
				//fprintf(stderr, "%.5f\n", sortedGenes[i][j].f);
			}
		}

		//fprintf(stderr, "Got here"); 

		vector<vector<float> > result;
		DoAggregate(sArgs, met, sortedGenes, queryGeneID, vecstrList.size(),
			goldstdGenePresence, excludeGene, includeGene, result);


		/*
		if(sArgs.fold_over_random_flag==1){
			vector<AResultFloat> *randomScores = new vector<AResultFloat>[vecstrList.size()];
			//should shuffle only within annotated genes
			for(i=0; i<vecstrList.size(); i++){
				GetRandom(rnd, sortedGenes[i], randomScores[i], excludeGene[i], 
					includeGene[i], queryGeneID[i], nan);
			}
			vector<vector<float> > random_result;
			DoAggregate(sArgs, met, randomScores, queryGeneID, vecstrList.size(),
				goldstdGenePresence, excludeGene, includeGene, random_result);

			vector<vector<float> > fold;
			fold.resize(result.size());
			for(i=0; i<result.size(); i++){
				fold[i] = vector<float>();
				fold[i].resize(result[i].size());
				for(j=0; j<result[i].size(); j++){
					if(random_result[i][j]==result[i][j]){
						fold[i][j] = 1.0;
					}else if(random_result[i][j]==0){
						fold[i][j] = 1.0;
					}else{
						fold[i][j] = result[i][j]/random_result[i][j];
					}
					fprintf(stderr, "%.2f %d %d %.2f %.2f\n", fold[i][j], i, j, result[i][j], random_result[i][j]);
				}
			}
			PrintResult(fold);
		*/
		//}else{
		PrintResult(result);	
		//}

	}

	if(qmode == MULTI_DWEIGHT){
		string dweightList = sArgs.dweight_list_arg;
		vector<string> vecstrList;
		CSeekTools::ReadListOneColumn(dweightList, vecstrList);
		int numDataset = 0;

		for(j=0; j<vecstrList.size(); j++){
			vector<float> ww;
			CSeekTools::ReadArray(vecstrList[j].c_str(), ww);
			vector<AResultFloat> sortedDatasets;
			sortedDatasets.resize(ww.size());
			numDataset = ww.size();
			for(i=0; i<sortedDatasets.size(); i++){
				sortedDatasets[i].i = i;
				sortedDatasets[i].f = ww[i];
			}
			sort(sortedDatasets.begin(), sortedDatasets.end());
			/*for(i=0; i<sortedDatasets.size(); i++){
				fprintf(stderr, "%.2e\t%s\n", sortedDatasets[i].f, 
					vecstrDatasets[sortedDatasets[i].i].c_str());
			}*/
			vector<int> depth;
			depth.push_back((int) ((float)numDataset * 0.005));	
			depth.push_back((int) ((float)numDataset * 0.01));	
			depth.push_back((int) ((float)numDataset * 0.05));	
			depth.push_back((int) ((float)numDataset * 0.10));	
			depth.push_back((int) ((float)numDataset * 0.20));	
			depth.push_back((int) ((float)numDataset * 0.50));	
			vector<float> avg;
			for(i=0; i<depth.size(); i++){
				float a = 0;
				for(k=0; k<depth[i]; k++)
					a+=sortedDatasets[k].f;
				a /= (float) depth[i];
				avg.push_back(a);
			}
			for(i=0; i<depth.size(); i++){
				fprintf(stderr, "%.3e\t", avg[i]);
			}
			fprintf(stderr, "\n");
		}

	}


#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }
