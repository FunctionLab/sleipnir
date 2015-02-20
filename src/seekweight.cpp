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
#include "seekweight.h"


namespace Sleipnir {

//MIN_REQUIRED is for a given gene A, the minimum query genes required that
//correlate with A in order to count A's query score
bool CSeekWeighter::LinearCombine(vector<utype> &rank,
	const vector<utype> &cv_query, CSeekDataset &sDataset,
	const bool bNegativeCor,
	const utype &MIN_REQUIRED, const bool &bSquareZ){

	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	if(mapQ==NULL) return true;

	if(cv_query.size()==0){
		cerr << "cv_query empty" << endl;
		return true;
	}
	utype i, j, g, q;
	vector<utype>::const_iterator iter;

	utype iNumGenes = sDataset.GetNumGenes();
	utype q_size = cv_query.size();
	utype **f = sDataset.GetDataMatrix();

	utype DEFAULT_NA = 0;
	if(bNegativeCor){
		DEFAULT_NA = 640;
	}

	CSeekTools::InitVector(rank, iNumGenes, (utype) DEFAULT_NA);

	vector<utype> queryPos;
	queryPos.resize(q_size);
	for(i=0; i<q_size; i++) queryPos[i] = mapQ->GetForward(cv_query[i]);

	sort(queryPos.begin(), queryPos.end());

	vector<utype> offset;
	offset.push_back(0);
	for(i=1; i<q_size; i++) offset.push_back(queryPos[i] - queryPos[i-1]);
	offset.resize(offset.size());

	vector<utype>::iterator iter_g;
	vector<utype>::const_iterator iterOffset;
	utype **pf;
	utype *pp;
	utype totNonZero, tmpScore;

	//special case, no MIN_REQUIRED
	if(MIN_REQUIRED==1){
		if(bSquareZ){ //must be used with cutoff (or else we will mix
					//positive and negative correlations together)
			for(iter_g=rank.begin(), pf = &f[0]; iter_g!=rank.end(); 
			iter_g++, pf++){
				for(totNonZero=0, tmpScore = 0, pp = &(*pf)[queryPos[0]],
				iterOffset = offset.begin(); iterOffset!=offset.end();
				iterOffset++){
					pp+=(*iterOffset);
					if((*pp)==0) continue;
					float sc = (float) ((*pp) - 320) / 100.0; 
					sc = fabs(sc) + 1.0; //add one adjustment, suitable if cutoff=0
					tmpScore += (utype) (sc * sc * 100.0 + 320.0);
					++totNonZero;
				}
				if(totNonZero==0) continue;
				//(*iter_g) = tmpScore / totNonZero;
				(*iter_g) = tmpScore / q_size;
			}		
		}else{
			for(iter_g=rank.begin(), pf = &f[0]; iter_g!=rank.end(); 
			iter_g++, pf++){
				for(totNonZero=0, tmpScore = 0, pp = &(*pf)[queryPos[0]],
				iterOffset = offset.begin(); iterOffset!=offset.end();
				iterOffset++){
					pp+=(*iterOffset);
					if((*pp)==0) continue;
					tmpScore += *pp;
					++totNonZero;
				}
				if(totNonZero==0) continue;
				//(*iter_g) = tmpScore / totNonZero;
				(*iter_g) = tmpScore / q_size;
			}
		}
	}
	else{ //MIN_REQUIRED ENABLED

		if(bSquareZ){
			//if the score of a gene to a query (*pp) is 0, it could 
			//either mean the gene is absent, or the gene-to-query correlation 
			//is below cutoff
			for(iter_g=rank.begin(), pf = &f[0]; iter_g!=rank.end(); 
			iter_g++, pf++){
				for(totNonZero=0, tmpScore = 0, pp = &(*pf)[queryPos[0]],
				iterOffset = offset.begin(); iterOffset!=offset.end();
				iterOffset++){
					pp+=(*iterOffset);
					if((*pp)==0) continue;
					float sc = (float) ((*pp) - 320) / 100.0;
					sc = fabs(sc) + 1.0; //add one adjustment, suitable for cutoff=0
					tmpScore += (utype) (sc * sc * 100.0 + 320.0);
					++totNonZero;
				}
				//if enough query edges passed the cut off 
				if(totNonZero >= MIN_REQUIRED)
					(*iter_g) = tmpScore / totNonZero;
				else
					(*iter_g) = DEFAULT_NA;
			}
		}
		else{
			for(iter_g=rank.begin(), pf = &f[0]; iter_g!=rank.end(); 
			iter_g++, pf++){
				for(totNonZero=0, tmpScore = 0, pp = &(*pf)[queryPos[0]],
				iterOffset = offset.begin(); iterOffset!=offset.end();
				iterOffset++){
					pp+=(*iterOffset);
					if((*pp)==0) continue;
					tmpScore += *pp;
					++totNonZero;
				}
				if(totNonZero >= MIN_REQUIRED)
					(*iter_g) = tmpScore / totNonZero;
				else
					(*iter_g) = DEFAULT_NA;
			}
		}
	}

	return true;
}

bool CSeekWeighter::OrderStatisticsPreCompute(){
	utype cc, dd, kk;
	vector<float> all;
	all.resize(1100*600*600);

	utype i = 0;
	for(kk=0; kk<22000; kk+=20){
		float p = (float) (kk + 1) / 22000;
		for(dd=0; dd<3000; dd+=5){
			unsigned int rrk = dd + 1;
			for(cc=0; cc<3000; cc+=5){
				double gg = gsl_cdf_binomial_Q(dd+1, p, cc+1);
				all[i] = -1.0*log(gg);
				i++;
			}
		}
	}
	//fprintf(stderr, "Done calculating\n"); //getchar();

	CSeekTools::WriteArray("/tmp/order_stats.binomial.bin", all);

	//fprintf(stderr, "Done saving\n"); //getchar();
}

bool CSeekWeighter::OrderStatisticsRankAggregation(const utype &iDatasets,
	const utype &iGenes, utype **rank_d, const vector<utype> &counts,
	vector<float> &master_rank, const utype &numThreads, const bool bNegativeCor){

	//bNegativeCor: 
	//do integration normally (ie based on positive correlations)
	//then reverse the final ranking to get negative correlated gene ranking
	float DEFAULT_NA = -320; //CAUTION!!

	//vector<float> precompute;
	//CSeekTools::ReadArray("/tmp/order_stats.binomial.bin", precompute);
	//fprintf(stderr, "Before\n"); getchar();
	//OrderStatisticsTest();

	if(rank_d==NULL){
		fprintf(stderr, "rank_d is null");
		return false;
	}

	master_rank.clear();
	master_rank.resize(iGenes);
	utype i, j, k, dd, d;

	//Zero out genes that are present in few datasets (<50%)
	for(k=0; k<iGenes; k++)
		if(counts[k]<(int)(0.5*iDatasets))
			for(j=0; j<iDatasets; j++) rank_d[j][k] = 0;

	//Hold the normalized rank
	float **rank_f =
		CSeekTools::Init2DArray(iDatasets, iGenes, (float) 1.1);

	for(j=0; j<iDatasets; j++){
		vector<AResult> this_d;
		utype numNonZero = 0;
		for(k=0; k<iGenes; k++){
			if(rank_d[j][k]>0) numNonZero++;
		}
		if(numNonZero==0){
			continue;
		}
		this_d.resize(numNonZero);
		int kk=0;
		for(k=0; k<iGenes; k++){
			if(rank_d[j][k]<=0) continue;
			this_d[kk].i = k;
			this_d[kk].f = rank_d[j][k];
			kk++;
		}
		sort(this_d.begin(), this_d.end());
		for(k=0; k<numNonZero; k++){
			rank_f[j][this_d[k].i] =
				(float) (k+1) / (float) numNonZero;
		}
		this_d.clear();
	}

	//fprintf(stderr, "Done1\n");
	vector<gsl_vector_float *> gss;
	vector<gsl_permutation *> perms;
	vector<gsl_permutation *> rks;

	gss.resize(numThreads);
	perms.resize(numThreads);
	rks.resize(numThreads);
	for(i=0; i<numThreads; i++){
		gss[i] = gsl_vector_float_calloc(iDatasets);
		perms[i] = gsl_permutation_alloc(iDatasets);
		rks[i] = gsl_permutation_alloc(iDatasets);
	}

	//gsl_vector_float *gs = gsl_vector_float_calloc(iDatasets);
	//gsl_permutation *perm = gsl_permutation_alloc(iDatasets);
	//gsl_permutation *rk = gsl_permutation_alloc(iDatasets);
	//fprintf(stderr, "Finished allocating\n");

	#pragma omp parallel for \
	shared(rank_f, gss, perms, rks) private(k, dd) \
	schedule(dynamic)
	for(k=0; k<iGenes; k++){
		utype tid = omp_get_thread_num();
		gsl_vector_float *gs = gss[tid];
		gsl_permutation *perm = perms[tid];
		gsl_permutation *rk = rks[tid];

		master_rank[k] = DEFAULT_NA;
		if(counts[k]<(int)(0.5*iDatasets)) continue;

		for(dd=0; dd<iDatasets; dd++)
			gsl_vector_float_set(gs, dd, rank_f[dd][k]);

		gsl_sort_vector_float_index(perm, gs);
		gsl_permutation_inverse(rk, perm);

		float max = DEFAULT_NA;
		int max_rank = -1;
		float max_p = -1;
		for(dd=0; dd<iDatasets; dd++){
			//get the prob of the gene in dset dd
			float p = gsl_vector_float_get(gs, dd);
			if(p<0 || p>1.0){
				continue;
			//	fprintf(stderr, "D%d %.5f\n", dd, p);
			}
			//get the rank of dset dd across all dsets
			unsigned int rrk = rk->data[dd];
			double gg = gsl_cdf_binomial_Q(rrk+1, p, counts[k]);
			float tmp = -1.0*log(gg);

			//load precomputed value==============
			/*int ind_p = (int) (p / (20.0 / 22000.0));
			int ind_r = (int) (rrk / 5.0);
			int ind_c = (int) (counts[k] / 5.0);
			float tmp = precompute[ind_p*600*600 + ind_r*600 + ind_c];
			*/
			//end=================================

			if(isinf(tmp)) tmp = DEFAULT_NA;
			if(tmp>max){
				max = tmp;
				max_rank = rrk;
				max_p = p;
			}
		}
		if(max!=DEFAULT_NA){
			master_rank[k] = max;
			//fprintf(stderr, "rank %.5f %.5f\n", max_p, max);
		}
	}

	CSeekTools::Free2DArray(rank_f, iDatasets);
	for(i=0; i<numThreads; i++){
		gsl_permutation_free(perms[i]);
		gsl_permutation_free(rks[i]);
		gsl_vector_float_free(gss[i]);
	}
	perms.clear();
	rks.clear();
	gss.clear();


	//REVERSE FINAL RANKING IF SORTING BY NEGATIVE CORRELATIONS
	if(bNegativeCor){
		DEFAULT_NA = 320;
		float max = -320;
		for(k=0; k<iGenes; k++){
			if(master_rank[k]==max){
				master_rank[k] = DEFAULT_NA;
			}
		}
	}

	//gsl_permutation_free(perm);
	//gsl_permutation_free(rk);
	//gsl_vector_float_free(gs);

	return true;
}

//for equal weighting, in case user wants to still see
//ordering of datasets, based on distance from the average ranking
bool CSeekWeighter::OneGeneWeighting(CSeekQuery &sQuery, 
	CSeekDataset &sDataset, const float &rate, 
	const float &percent_required, const bool &bSquareZ,
	vector<utype> *rrank, const CSeekQuery *goldStd, const bool bNegativeCor){

	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	if(mapQ==NULL) return true;

	sDataset.InitializeCVWeight(1);

	utype i, j, qi, qj;

	vector<char> is_query, is_gold;
	CSeekTools::InitVector(is_query, sDataset.GetNumGenes(), (char) 0);
	CSeekTools::InitVector(is_gold, sDataset.GetNumGenes(), (char) 0);

	vector<utype> &rank = *rrank;

	utype TOP = 1000;
	//utype TOP = 0;
	vector<AResult> ar;
	ar.resize(rank.size());

	const vector<utype> &allQ = sQuery.GetQuery();
	vector<utype> query;
	utype num_q = 0;
	utype num_v = 0;

	const vector<utype> &vi = sQuery.GetCVQuery(qi);

	/* Set up query and gold standard */
	for(i=0; i<allQ.size(); i++){
		if(CSeekTools::IsNaN(mapQ->GetForward(allQ[i]))) continue;
		is_query[allQ[i]] = 1;
		query.push_back(allQ[i]);
		num_q++;
	}

	//assume goldStd must not be null
	assert(goldStd!=NULL);
	const vector<utype> &allGoldStd = goldStd->GetQuery();
	for(i=0; i<allGoldStd.size(); i++){
		if(CSeekTools::IsNaN(mapG->GetForward(allGoldStd[i])))
			continue;
		if(is_query[allGoldStd[i]]==1) continue;
		is_gold[allGoldStd[i]] = 1;
		num_v++;
	}

	if(num_q==0 || num_v==0){
		sDataset.SetCVWeight(0, -1);
	}else{
		/* actual weighting */
		float w = 0;
		const utype MIN_QUERY_REQUIRED =
			max((utype) 1, (utype) (percent_required * query.size()));
		bool ret = LinearCombine(rank, query, sDataset,
			bNegativeCor, MIN_QUERY_REQUIRED, bSquareZ);
		ret = CSeekPerformanceMeasure::RankBiasedPrecision(rate,
			rank, w, is_query, is_gold, *mapG, &ar, bNegativeCor, TOP);
		if(!ret) sDataset.SetCVWeight(0, -1);
		else sDataset.SetCVWeight(0, w);
	}

	ar.clear();
	return true;
}

bool CSeekWeighter::AverageWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset,
	const float &percent_required, const bool &bSquareZ, float &w, 
	const bool bNegativeCor){

	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	if(mapQ==NULL) return true;

	utype i, j, qi, qj;

	const vector<utype> &allQ = sQuery.GetQuery();
	vector<utype> presentQ;
	utype num_q = 0;
	for(qi=0; qi<allQ.size(); qi++){
		if(CSeekTools::IsNaN(mapQ->GetForward(allQ[qi]))) continue;
		num_q++;
		presentQ.push_back(allQ[qi]);
	}

	w = 0;
	const utype MIN_QUERY_REQUIRED =
		max((utype) 2, (utype) (percent_required * allQ.size()));

	if(num_q<MIN_QUERY_REQUIRED){
		w = -1;
		return true;
	}

	utype **f = sDataset.GetDataMatrix();
	/* as long as rank[g] does not overflow, due to too many queries, we are fine
	 * should control query size to be <100. */
	vector<utype> queryPos;
	queryPos.resize(num_q);
	for(i=0; i<num_q; i++)
		queryPos[i] = mapQ->GetForward(presentQ[i]);

	int pairs = 0;
	if(bSquareZ){
		for(qi=0; qi<presentQ.size(); qi++){
			for(qj=0; qj<presentQ.size(); qj++){
				if(qi==qj) continue;
				float t = (float) f[presentQ[qj]][queryPos[qi]];
				t = (t-320)/100.0;
				t = t * t;
				w += t;
				pairs++;
			}
		}
		w /= (float) pairs;

	}else{
		for(qi=0; qi<presentQ.size(); qi++){
			for(qj=0; qj<presentQ.size(); qj++){
				if(qi==qj) continue;
				w += (float) f[presentQ[qj]][queryPos[qi]];
				pairs++;
			}
		}
		w /= (float) pairs;
		w /= (float) 640;
	}

	if(bNegativeCor){
		w = w * -1.0;
	}

	return true;
}

bool CSeekWeighter::CVWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset,
	const float &rate, const float &percent_required, const bool &bSquareZ,
	const bool bNegativeCor,
	vector<utype> *rrank, const CSeekQuery *goldStd){

	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	if(mapQ==NULL) return true;

	utype iFold = sQuery.GetNumFold();
	sDataset.InitializeCVWeight(iFold);

	utype i, j, qi, qj;

	vector<char> is_query_cross, is_gold;
	CSeekTools::InitVector(is_query_cross, sDataset.GetNumGenes(), (char) 0);
	CSeekTools::InitVector(is_gold, sDataset.GetNumGenes(), (char) 0);

	vector<utype> &rank = *rrank;

	utype TOP = 1000;
	//utype TOP = 0; //disable TOP approximation
	vector<AResult> ar;
	ar.resize(rank.size());

	const vector<utype> &allQ = sQuery.GetQuery();

	for(qi=0; qi<iFold; qi++){
		vector<utype> cv_query;
		utype num_q = 0;
		utype num_v = 0;

		const vector<utype> &vi = sQuery.GetCVQuery(qi);

		/* Set up query and gold standard */
		for(i=0; i<vi.size(); i++){
			if(CSeekTools::IsNaN(mapQ->GetForward(vi[i]))) continue;
			is_query_cross[vi[i]] = 1;
			cv_query.push_back(vi[i]);
			num_q++;
		}

		if(goldStd==NULL){ //if no custom gold standard
			//Use query itself as gold standard
			for(i=0; i<allQ.size(); i++){
				if(CSeekTools::IsNaN(mapQ->GetForward(allQ[i]))) continue;
				if(is_query_cross[allQ[i]]==1) continue;
				is_gold[allQ[i]] = 1;
				num_v++;
			}
		}else{
			const vector<utype> &allGoldStd = goldStd->GetQuery();
			//Use custom gene-set as gold standard
			for(i=0; i<allGoldStd.size(); i++){
				if(CSeekTools::IsNaN(mapG->GetForward(allGoldStd[i])))
					continue;
				if(is_query_cross[allGoldStd[i]]==1) continue;
				is_gold[allGoldStd[i]] = 1;
				num_v++;
			}
		}

		if(num_q==0 || num_v==0){
			sDataset.SetCVWeight(qi, -1);
			//printf("num_q %d or num_v %d\n", num_q, num_v);
		}else{
			/* actual weighting */
			float w = 0;
			const utype MIN_QUERY_REQUIRED =
				max((utype) 1, (utype) (percent_required * cv_query.size()));
			bool ret = LinearCombine(rank, cv_query, sDataset, bNegativeCor,
				MIN_QUERY_REQUIRED, bSquareZ);
			ret = CSeekPerformanceMeasure::RankBiasedPrecision(rate,
				rank, w, is_query_cross, is_gold, *mapG, &ar, bNegativeCor, TOP);
			//fprintf(stderr, "Weight %.5f\n", w);
			//ret = CSeekPerformanceMeasure::AveragePrecision(
			//	rank, w, is_query_cross, is_gold, *mapG, &ar);
			if(!ret) sDataset.SetCVWeight(qi, -1);
			else sDataset.SetCVWeight(qi, w);
			//printf("Weight: %.5f\n", w);
		}

		/* Reset query and gold standard */
		for(i=0; i<vi.size(); i++){
			if(CSeekTools::IsNaN(mapQ->GetForward(vi[i]))) continue;
			is_query_cross[vi[i]] = 0;
		}

		if(goldStd==NULL){
			for(i=0; i<allQ.size(); i++){
				if(CSeekTools::IsNaN(mapQ->GetForward(allQ[i]))) continue;
				is_gold[allQ[i]]=0;
			}
		}else{
			const vector<utype> &allGoldStd = goldStd->GetQuery();
			for(i=0; i<allGoldStd.size(); i++){
				if(CSeekTools::IsNaN(mapG->GetForward(allGoldStd[i])))
					continue;
				is_gold[allGoldStd[i]] = 0;
			}
		}

	}

	ar.clear();

	return true;
}

}
