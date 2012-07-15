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
	SINGLE_QUERY, MULTI_QUERY
};

enum METRIC{
	RBP, AVGP, PR, PR_ALL, AUC
};

float RankBiasedPrecision(const float &rate,
	const vector<AResultFloat> &sortedScore, const vector<char> &gold_std,
	const float &nanVal){

	vector<AResultFloat>::const_iterator iterScore = sortedScore.begin();
	float x = 0;
	ushort i;
	for(i=0; iterScore!=sortedScore.end(); i++, iterScore++){
		if(iterScore->f == nanVal) continue;
		if(gold_std[iterScore->i]==1) x += pow(rate, i);
	}
	x*=(1.0-rate);
	return x;
}

vector<float>* Precision(const vector<AResultFloat> &sortedScore,
	const vector<char> &gold_std, const float &nanVal){

	ushort i, numPos = 1;
	vector<float> *r = new vector<float>();

	vector<AResultFloat>::const_iterator iterScore = sortedScore.begin();
	for(i=0; iterScore!=sortedScore.end(); i++, iterScore++){
		if(iterScore->f == nanVal) continue;
		if(gold_std[iterScore->i]==1) r->push_back((float)
			numPos++ / (float) (i + 1));
	}
	r->resize(r->size());
	return r;
}

bool MeanStandardDeviation(const vector<float> &f, float &avg, float &stdev){
	avg = 0;
	stdev = 0;
	avg = std::accumulate(f.begin(), f.end(), 0);
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
		fprintf(stdout, "%.5f\n", rbp);
	}
	else if(met==AVGP || met==PR){
		vector<float> *vf = Precision(sortedGenes, goldstdGenePresence,
			nan);
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
			if(X>=vf->size()){
				fprintf(stderr, "Error: X is too large (>=%d)\n", vf->size());
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
			X = (int) ( per * (float) vf->size() );
			X = std::min(0, X - 1);
		}

		if(met==AVGP){
			float avg = std::accumulate(vf->begin(), vf->begin()+X, 0);
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
			return 0;
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
	size_t				i, j;

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
			cerr << "Ignoring line: " << acBuffer << endl;
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

	enum METRIC met;
	if(sArgs.rbp_flag==1) met = RBP;
	else if(sArgs.avgp_flag==1) met = AVGP;
	else if(sArgs.pr_flag==1) met = PR;
	else if(sArgs.pr_all_flag==1) met = PR_ALL;
	else if(sArgs.auc_flag==1) met = AUC;

	if(qmode==SINGLE_QUERY){
		string goldstdFile = sArgs.goldstd_arg;
		vector<string> goldstdGenes;
		CSeekTools::ReadMultiGeneOneLine(goldstdFile, goldstdGenes);
		vector<char> goldstdGenePresence;
		CSeekTools::InitVector(goldstdGenePresence,
			vecstrGenes.size(), (char) 0);

		for(i=0; i<goldstdGenes.size(); i++)
			goldstdGenePresence[mapstriGenes[goldstdGenes[i]]] = 1;

		string queryFile = sArgs.query_arg;
		vector<string> queryGenes;
		CSeekTools::ReadMultiGeneOneLine(queryFile, queryGenes);
		vector<ushort> queryGeneID;
		for(i=0; i<queryGenes.size(); i++)
			queryGeneID.push_back(mapstriGenes[queryGenes[i]]);

		string genescoreFile = sArgs.gscore_arg;
		vector<float> geneScores;
		CSeekTools::ReadArray(genescoreFile.c_str(), geneScores);
		float maxScore = *std::max_element(geneScores.begin(),
			geneScores.end());

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

		if(met!=PR_ALL){
			float eval;
			bool ret = EvaluateOneQuery(sArgs, met, sortedGenes,
				goldstdGenePresence, nan, eval);
			if(!ret) return 1;
			fprintf(stderr, "%.5f\n", eval);
			return 0;
		}else{
			vector<float> evalAll;
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

		string queryList = sArgs.query_list_arg;
		vecstrList.clear();
		CSeekTools::ReadListOneColumn(queryList, vecstrList);
		vector<ushort> *queryGeneID = new vector<ushort>[vecstrList.size()];
		for(i=0; i<vecstrList.size(); i++){
			vector<string> queryGenes;
			CSeekTools::ReadMultiGeneOneLine(vecstrList[i], queryGenes);
			for(j=0; j<queryGenes.size(); j++)
				queryGeneID[i].push_back(mapstriGenes[queryGenes[j]]);
		}

		string genescoreList = sArgs.gscore_list_arg;
		vecstrList.clear();
		CSeekTools::ReadListOneColumn(genescoreList, vecstrList);
		vector<float> *geneScores = new vector<float>[vecstrList.size()];
		float *maxScore = new float[vecstrList.size()];
		vector<AResultFloat> *sortedGenes =
			new vector<AResultFloat>[vecstrList.size()];
		float nan = sArgs.nan_arg;
		for(i=0; i<vecstrList.size(); i++){
			CSeekTools::ReadArray(vecstrList[i].c_str(), geneScores[i]);
			maxScore[i] = *std::max_element(geneScores[i].begin(),
				geneScores[i].end());
			sortedGenes[i].resize(geneScores[i].size());
			for(j=0; j<sortedGenes[i].size(); j++){
				sortedGenes[i][j].i = j;
				sortedGenes[i][j].f = geneScores[i][j];
			}
		}

		if(sArgs.agg_ranksum_flag==1 || sArgs.agg_scoresum_flag==1){
			//Gold standard must be the same across all queries for this mode!!
			//ASSUME THIS IS TRUE

			for(i=0; i<vecstrList.size(); i++)
				for(j=0; j<queryGeneID[i].size(); j++)
					sortedGenes[i][queryGeneID[i][j]].f = maxScore[i];

			vector<AResultFloat> *master = NULL;
			vector<AResultFloat> master_score;
			vector<AResultFloat> master_rank;

			if(sArgs.agg_scoresum_flag==1){
				master_score.resize(sortedGenes[0].size());
				for(j=0; j<sortedGenes[0].size(); j++){
					int count_nan = 0;
					master_score[j].i = j;
					master_score[j].f = 0;
					for(i=0; i<vecstrList.size(); i++){
						if(sortedGenes[i][j].f==nan){
							count_nan++;
							continue;
						}
						master_score[j].f += sortedGenes[i][j].f;
					}
					if(count_nan>0) master_score[j].f = -320;
					else master_score[j].f /= (float)vecstrList.size();
				}

				sort(master_score.begin(), master_score.end());

				master = &master_score;
			}

			else if(sArgs.agg_ranksum_flag==1){
				vector<AResultFloat> master_rank;
				master_rank.resize(sortedGenes[0].size());

				for(i=0; i<vecstrList.size(); i++)
					sort(sortedGenes[i].begin(), sortedGenes[i].end());

				for(j=0; j<sortedGenes[0].size(); j++){
					master_rank[j].i = j;
					master_rank[j].f = 0;
				}

				for(i=0; i<vecstrList.size(); i++)
					for(j=0; j<sortedGenes[i].size(); j++)
						master_rank[sortedGenes[i][j].i].f +=
							sortedGenes[i].size() - j;

				for(j=0; j<sortedGenes[0].size(); j++)
					master_rank[j].f /= (float) vecstrList.size();

				sort(master_rank.begin(), master_rank.end());
				master = &master_rank;
			}

			if(met!=PR_ALL){
				float eval;
				bool ret = EvaluateOneQuery(sArgs, met, *master,
					goldstdGenePresence[0], CMeta::GetNaN(), eval);
				if(!ret) return 1;
				fprintf(stderr, "%.5f\n", eval);
				return 0;
			}else{
				vector<float> evalAll;
				bool ret = EvaluateOneQuery(sArgs, met, *master,
					goldstdGenePresence[0], CMeta::GetNaN(), evalAll);
				if(!ret) return 1;
				PrintVector(evalAll);
				return 0;
			}
		}

		//Not Rank Sum
		vector<float> eval;
		vector< vector<float> > vecevalAll;

		eval.resize(vecstrList.size());
		vecevalAll.resize(vecstrList.size());

		for(i=0; i<vecstrList.size(); i++){
			for(j=0; j<queryGeneID[i].size(); j++)
				sortedGenes[i][queryGeneID[i][j]].f = nan;
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

		if(met!=PR_ALL){
			if(sArgs.agg_avg_flag==1){
				float avg, stdev;
				MeanStandardDeviation(eval, avg, stdev);
				fprintf(stderr, "%.5f %.5f\n", avg, stdev);
				return 0;
			}
			if(sArgs.agg_quartile_flag==1){
				float min, max, Q1, Q2, Q3;
				MinMaxQuartile(eval, min, max, Q1, Q2, Q3);
				fprintf(stderr, "%.5f %.5f\n", min, max);
				fprintf(stderr, "%.5f %.5f %.5f\n", Q1, Q2, Q3);
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
				PrintVector(avgAll);
				PrintVector(stdevAll);
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
				PrintVector(minAll);
				PrintVector(maxAll);
				PrintVector(Q1All);
				PrintVector(Q2All);
				PrintVector(Q3All);
				return 0;
			}
		}
	}


#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }
