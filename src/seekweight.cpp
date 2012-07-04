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

bool CSeekWeighter::LinearCombine(vector<ushort> &rank, const vector<ushort> &cv_query,
	CSeekDataset &sDataset, const bool bAllocate){
	if(cv_query.size()==0){
		cerr << "cv_query empty" << endl;
		return true;
	}
	ushort i, j;
	ushort g, q;
	vector<ushort>::const_iterator iter;

	ushort iNumGenes = sDataset.GetNumGenes();
	ushort q_size = cv_query.size();
	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	ushort **f = sDataset.GetDataMatrix();

	if(bAllocate){
		CSeekTools::InitVector(rank, iNumGenes, (ushort) 0);
	}else{
		fill(rank.begin(), rank.end(), 0);
	}

	/* as long as rank[g] does not overflow, due to too many queries, we are fine
	 * should control query size to be <100. */
	vector<ushort> queryPos;
	queryPos.resize(q_size);
	for(i=0; i<q_size; i++){
		queryPos[i] = mapQ->GetForward(cv_query[i]);
	}

	for(g=0; g<iNumGenes; g++){
		for(iter=queryPos.begin(); iter!=queryPos.end(); iter++){
			rank[g] += f[g][*iter];
		}
		rank[g] /= q_size;
	}
	return true;
}


bool CSeekWeighter::CVWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset){
	ushort iFold = sQuery.GetNumFold();
	sDataset.InitializeCVWeight(iFold);

	ushort i, j, qi, qj;

	vector<char> is_query_cross, is_gold;
	CSeekTools::InitVector(is_query_cross, sDataset.GetNumGenes(), (char) 0);
	CSeekTools::InitVector(is_gold, sDataset.GetNumGenes(), (char) 0);

	vector<ushort> rank;
	CSeekTools::InitVector(rank, sDataset.GetNumGenes(), (ushort) 0);

	ushort TOP = 1000;
	vector<AResult> ar;
	ar.resize(rank.size());

	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
	vector<ushort> &allQ = sQuery.GetQuery();

	for(qi=0; qi<iFold; qi++){
		vector<ushort> cv_query;
		ushort num_q = 0;
		ushort num_v = 0;

		vector<ushort> &vi = sQuery.GetCVQuery(qi);

		/* Set query and gold standard */
		for(i=0; i<vi.size(); i++){
			if(CSeekTools::IsNaN(mapQ->GetForward(vi[i]))) continue;
			is_query_cross[vi[i]] = 1;
			cv_query.push_back(vi[i]);
			num_q++;
		}

		for(i=0; i<allQ.size(); i++){
			if(CSeekTools::IsNaN(mapQ->GetForward(allQ[i]))) continue;
			if(is_query_cross[allQ[i]]==1) continue;
			is_gold[allQ[i]] = 1;
			num_v++;
		}

		if(num_q==0 || num_v==0){
			sDataset.SetCVWeight(qi, -1);
			//printf("num_q %d or num_v %d\n", num_q, num_v);
		}else{
			/* actual weighting */
			//vector<unsigned short> rank;
			float w = 0;
			bool ret = LinearCombine(rank, cv_query, sDataset, false);
			ret = CSeekPerformanceMeasure::RankBiasedPrecision(0.95,
				rank, w, is_query_cross, is_gold, *mapG, false, &ar, TOP);
			if(!ret){
				sDataset.SetCVWeight(qi, -1);
			}else{
				sDataset.SetCVWeight(qi, w);
			}
			//printf("Weight: %.5f\n", w);
		}
		/* Reset query and gold standard */
		for(i=0; i<vi.size(); i++){
			if(CSeekTools::IsNaN(mapQ->GetForward(vi[i]))) continue;
			is_query_cross[vi[i]] = 0;
		}
		for(i=0; i<allQ.size(); i++){
			if(CSeekTools::IsNaN(mapQ->GetForward(allQ[i]))) continue;
			is_gold[allQ[i]] = 0;
		}

	}
	return true;
}

}
