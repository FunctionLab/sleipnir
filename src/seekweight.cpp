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
#include "seekweight.h"
#include "seekreader.h"
#include "seekquery.h"
#include "seekevaluate.h"

namespace Sleipnir {

//MIN_REQUIRED is for a given gene A, the minimum query genes required that
//correlate with A in order to count A's query score
bool CSeekWeighter::LinearCombine(vector<ushort> &rank,
	const vector<ushort> &cv_query, CSeekDataset &sDataset,
	const ushort &MIN_REQUIRED){

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

	//if the score of a gene to a query (*pp) is 0, it could either mean
	//the gene is absent, or the gene-to-query correlation is below cutoff
	for(iter_g=rank.begin(), pf = &f[0]; iter_g!=rank.end(); iter_g++, pf++){
		for(totNonZero=0, tmpScore = 0, pp = &(*pf)[queryPos[0]],
			iterOffset = offset.begin(); iterOffset!=offset.end();
			iterOffset++, pp+=(*iterOffset)){
			tmpScore += *pp;
			if((*pp)>0) ++totNonZero;
		}
		//if enough query edges passed the cut off
		if(totNonZero >= MIN_REQUIRED)
			(*iter_g) = tmpScore / totNonZero;
		else
			(*iter_g) = 0;
		//(*iter_g) = tmpScore / q_size;
	}

	return true;
}

bool CSeekWeighter::OrderStatisticsRankAggregation(const ushort &iDatasets,
	const ushort &iGenes, ushort **rank_d, const vector<ushort> &counts,
	vector<float> &master_rank){
	if(rank_d==NULL){
		fprintf(stderr, "rank_d is null");
		return false;
	}

	gsl_vector_float *gs = gsl_vector_float_calloc(iDatasets);
	gsl_permutation *perm = gsl_permutation_alloc(iDatasets);
	gsl_permutation *rk = gsl_permutation_alloc(iDatasets);

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
		this_d.resize(iGenes);
		for(k=0; k<iGenes; k++){
			this_d[k].i = k;
			this_d[k].f = rank_d[j][k];
			if(rank_d[j][k]>0) numNonZero++;
		}
		if(numNonZero==0){
			this_d.clear();
			continue;
		}
		sort(this_d.begin(), this_d.end());
		for(k=0; k<iGenes; k++){
			if(this_d[k].f==0) break;
			rank_f[j][this_d[k].i] =
				(float) (k+1) / (float) numNonZero;
		}
		this_d.clear();
	}

	for(k=0; k<iGenes; k++){
		master_rank[k] = DEFAULT_NULL;
		if(counts[k]<(int)(0.5*iDatasets)) continue;

		for(dd=0; dd<iDatasets; dd++)
			gsl_vector_float_set(gs, dd, rank_f[dd][k]);

		gsl_sort_vector_float_index(perm, gs);
		gsl_permutation_inverse(rk, perm);

		float max = DEFAULT_NULL;
		for(dd=0; dd<iDatasets; dd++){
			if(rank_f[dd][k]==1.1) continue;
			//get the prob of the gene in dset dd
			float p = gsl_vector_float_get(gs, dd);
			//get the rank of dset dd across all dsets
			unsigned int rrk = rk->data[dd];
			double gg = gsl_cdf_binomial_Q(rrk+1, p, counts[k]);
			float tmp = -1.0*log(gg);
			if(isinf(tmp)) tmp = DEFAULT_NULL;
			if(tmp>max) max = tmp;
		}
		if(max!=DEFAULT_NULL)
			master_rank[k] = max;
	}

	CSeekTools::Free2DArray(rank_f);
	gsl_permutation_free(perm);
	gsl_permutation_free(rk);
	gsl_vector_float_free(gs);

	return true;
}


bool CSeekWeighter::CVWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset,
	const float &rate, const float &percent_required,
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
				MIN_QUERY_REQUIRED);
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
