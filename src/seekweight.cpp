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

#include "seekweight.h"


namespace Sleipnir {

//MIN_REQUIRED is for a given gene A, the minimum query genes required that
//correlate with A in order to count A's query score
bool CSeekWeighter::LinearCombine(vector<ushort> &rank,
	const vector<ushort> &cv_query, CSeekDataset &sDataset,
	const ushort &MIN_REQUIRED, const bool &bSquareZ){

	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	if(mapQ==NULL) return true;

	if(cv_query.size()==0){
		cerr << "cv_query empty" << endl;
		return true;
	}
	ushort i, j, g, q;
	vector<ushort>::const_iterator iter;

	ushort iNumGenes = sDataset.GetNumGenes();
	ushort q_size = cv_query.size();
	ushort **f = sDataset.GetDataMatrix();

	rank.resize(iNumGenes);

	/* as long as rank[g] does not overflow, due to too many queries, we are fine
	 * should control query size to be <100. */
	vector<ushort> queryPos;
	queryPos.resize(q_size);
	for(i=0; i<q_size; i++) queryPos[i] = mapQ->GetForward(cv_query[i]);

	sort(queryPos.begin(), queryPos.end());

	vector<ushort> offset;
	offset.push_back(0);
	for(i=1; i<q_size; i++) offset.push_back(queryPos[i] - queryPos[i-1]);
	offset.resize(offset.size());

	vector<ushort>::iterator iter_g;
	vector<ushort>::const_iterator iterOffset;
	ushort **pf;
	ushort *pp;
	ushort totNonZero, tmpScore;

	//special case, no MIN_REQUIRED
	if(MIN_REQUIRED==1){
		if(bSquareZ){ //must be used with cutoff (or else we will mix
					//positive and negative correlations together)
			for(iter_g=rank.begin(), pf = &f[0]; iter_g!=rank.end(); 
			iter_g++, pf++){
				for(totNonZero=0, tmpScore = 0, pp = &(*pf)[queryPos[0]],
				iterOffset = offset.begin(); iterOffset!=offset.end();
				iterOffset++, pp+=(*iterOffset)){
					if((*pp)==0) continue;
					float sc = (float) ((*pp) - 320) / 100.0; 
					sc = fabs(sc) + 1.0; //add one adjustment, suitable if cutoff=0
					tmpScore += (ushort) (sc * sc * 100.0 + 320.0);
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
				iterOffset++, pp+=(*iterOffset)){
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
				iterOffset++, pp+=(*iterOffset)){
					if((*pp)==0) continue;
					float sc = (float) ((*pp) - 320) / 100.0;
					sc = fabs(sc) + 1.0; //add one adjustment, suitable for cutoff=0
					tmpScore += (ushort) (sc * sc * 100.0 + 320.0);
					++totNonZero;
				}
				//if enough query edges passed the cut off 
				if(totNonZero >= MIN_REQUIRED)
					(*iter_g) = tmpScore / totNonZero;
				else
					(*iter_g) = 0;
			}
		}
		else{
			for(iter_g=rank.begin(), pf = &f[0]; iter_g!=rank.end(); 
			iter_g++, pf++){
				for(totNonZero=0, tmpScore = 0, pp = &(*pf)[queryPos[0]],
				iterOffset = offset.begin(); iterOffset!=offset.end();
				iterOffset++, pp+=(*iterOffset)){
					if((*pp)==0) continue;
					tmpScore += *pp;
					++totNonZero;
				}
				if(totNonZero >= MIN_REQUIRED)
					(*iter_g) = tmpScore / totNonZero;
				else
					(*iter_g) = 0;
			}
		}
	}

	return true;
}

bool CSeekWeighter::OrderStatisticsPreCompute(){
	ushort cc, dd, kk;
	vector<float> all;
	all.resize(1100*600*600);

	ushort i = 0;
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

bool CSeekWeighter::OrderStatisticsRankAggregation(const ushort &iDatasets,
	const ushort &iGenes, ushort **rank_d, const vector<ushort> &counts,
	vector<float> &master_rank, const ushort &numThreads){

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
	ushort i, j, k, dd, d;

	//Zero out genes that are present in few datasets (<50%)
	for(k=0; k<iGenes; k++)
		if(counts[k]<(int)(0.5*iDatasets))
			for(j=0; j<iDatasets; j++) rank_d[j][k] = 0;

	//Hold the normalized rank
	float **rank_f =
		CSeekTools::Init2DArray(iDatasets, iGenes, (float) 1.1);

	const float DEFAULT_NULL = -320;

	for(j=0; j<iDatasets; j++){
		vector<AResult> this_d;
		ushort numNonZero = 0;
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
	shared(master_rank, counts, rank_f, gss, perms, rks, iGenes, \
	iDatasets) private(k, dd) \
	schedule(dynamic)
	for(k=0; k<iGenes; k++){
		ushort tid = omp_get_thread_num();
		gsl_vector_float *gs = gss[tid];
		gsl_permutation *perm = perms[tid];
		gsl_permutation *rk = rks[tid];

		master_rank[k] = DEFAULT_NULL;
		if(counts[k]<(int)(0.5*iDatasets)) continue;

		for(dd=0; dd<iDatasets; dd++)
			gsl_vector_float_set(gs, dd, rank_f[dd][k]);

		gsl_sort_vector_float_index(perm, gs);
		gsl_permutation_inverse(rk, perm);

		float max = DEFAULT_NULL;
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

			if(isinf(tmp)) tmp = DEFAULT_NULL;
			if(tmp>max){
				max = tmp;
				max_rank = rrk;
				max_p = p;
			}
		}
		if(max!=DEFAULT_NULL){
			master_rank[k] = max;
			//fprintf(stderr, "rank %.5f %.5f\n", max_p, max);
		}
	}

	CSeekTools::Free2DArray(rank_f);
	for(i=0; i<numThreads; i++){
		gsl_permutation_free(perms[i]);
		gsl_permutation_free(rks[i]);
		gsl_vector_float_free(gss[i]);
	}
	perms.clear();
	rks.clear();
	gss.clear();

	//gsl_permutation_free(perm);
	//gsl_permutation_free(rk);
	//gsl_vector_float_free(gs);

	return true;
}


bool CSeekWeighter::CVWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset,
	const float &rate, const float &percent_required, const bool &bSquareZ,
	vector<ushort> *rrank, const CSeekQuery *goldStd){

	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	if(mapQ==NULL) return true;

	ushort iFold = sQuery.GetNumFold();
	sDataset.InitializeCVWeight(iFold);

	ushort i, j, qi, qj;

	vector<char> is_query_cross, is_gold;
	CSeekTools::InitVector(is_query_cross, sDataset.GetNumGenes(), (char) 0);
	CSeekTools::InitVector(is_gold, sDataset.GetNumGenes(), (char) 0);

	vector<ushort> &rank = *rrank;

	ushort TOP = 1000;
	vector<AResult> ar;
	ar.resize(rank.size());

	const vector<ushort> &allQ = sQuery.GetQuery();

	for(qi=0; qi<iFold; qi++){
		vector<ushort> cv_query;
		ushort num_q = 0;
		ushort num_v = 0;

		const vector<ushort> &vi = sQuery.GetCVQuery(qi);

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
			const vector<ushort> &allGoldStd = goldStd->GetQuery();
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
			const ushort MIN_QUERY_REQUIRED =
				max((ushort) 1, (ushort) (percent_required * cv_query.size()));
			bool ret = LinearCombine(rank, cv_query, sDataset,
				MIN_QUERY_REQUIRED, bSquareZ);
			ret = CSeekPerformanceMeasure::RankBiasedPrecision(rate,
				rank, w, is_query_cross, is_gold, *mapG, &ar, TOP);
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
			const vector<ushort> &allGoldStd = goldStd->GetQuery();
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
