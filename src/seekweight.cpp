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

bool CSeekWeighter::LinearCombine(vector<short> &rank, vector<int> &cv_query,
	CSeekDataset &sDataset){
	if(cv_query.size()==0){
		cerr << "cv_query empty" << endl;
		return true;
	}
	size_t iNumGenes = sDataset.GetNumGenes();
	CSeekTools::InitVector(rank, iNumGenes, (short)-32768);
	size_t i, j, k;

	int q_size = cv_query.size();
	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	CFullMatrix<short> *f = sDataset.GetDataMatrix();

	size_t iGenesPresent = mapG->GetNumSet();
	/* as long as rank[g] does not overflow, due to too many queries, we are fine
	 * should control query size to be <100. */
	for(i=0; i<iGenesPresent; i++){
		size_t g = mapG->GetReverse(i);
		rank[g] = 0;
		for(j=0; j<q_size; j++){
			int qq = cv_query[j];
			if(g==qq) continue;
			size_t q = mapQ->GetForward(qq);
			rank[g] += f->Get(g, q);
		}
		rank[g] /= (short) q_size;
	}
	return true;
}


bool CSeekWeighter::CVWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset){
	size_t iFold = sQuery.GetNumFold();
	sDataset.InitializeCVWeight(iFold);

	int i, j, qi, qj;

	vector<char> is_query_cross, is_gold;
	CSeekTools::InitVector(is_query_cross, sDataset.GetNumGenes(), (char) 0);
	CSeekTools::InitVector(is_gold, sDataset.GetNumGenes(), (char) 0);

	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	for(qi=0; qi<iFold; qi++){
		vector<int> vi = sQuery.GetCVQuery(qi);
		vector<int> cv_query;
		CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
		int num_q = 0;
		int num_v = 0;

		/* Set query and gold standard */
		for(i=0; i<vi.size(); i++){
			if(mapQ->GetForward(vi[i])==-1) continue;
			is_query_cross[vi[i]] = 1;
			cv_query.push_back(vi[i]);
			num_q++;
		}

		vector<int> allQ = sQuery.GetQuery();
		for(i=0; i<allQ.size(); i++){
			if(mapQ->GetForward(allQ[i])==-1) continue;
			if(is_query_cross[allQ[i]]==1) continue;
			is_gold[allQ[i]] = 1;
			num_v++;
		}

		if(num_q==0 || num_v==0){
			sDataset.SetCVWeight(qi, -1);
		}else{
			/* actual weighting */
			vector<short> rank;
			float w = 0;
			bool ret = LinearCombine(rank, cv_query, sDataset);
			ret = CSeekPerformanceMeasure::RankBiasedPrecision(0.95,
				rank, w, is_query_cross, is_gold, *mapG);
			if(!ret){
				sDataset.SetCVWeight(qi, -1);
			}else{
				sDataset.SetCVWeight(qi, w);
			}
			//printf("Weight: %.5f\n", w);
		}
		/* Reset query and gold standard */
		for(i=0; i<vi.size(); i++){
			if(mapQ->GetForward(vi[i])==-1) continue;
			is_query_cross[vi[i]] = 0;
		}
		for(i=0; i<allQ.size(); i++){
			if(mapQ->GetForward(allQ[i])==-1) continue;
			is_gold[allQ[i]] = 0;
		}

	}
	return true;
}

}
