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

bool CSeekWeighter::LinearCombine(vector<float> &rank, vector<int> &cv_query,
	CSeekDataset &sDataset){
	if(cv_query.size()==0){
		cerr << "cv_query empty" << endl;
		return true;
	}
	size_t iNumGenes = sDataset.GetNumGenes();

	vector<float> new_rank;
	CSeekTools::InitVector(rank, iNumGenes, (float)0);
	CSeekTools::InitVector(new_rank, iNumGenes, (float)0);
	size_t i, j, k;

	int q_size = cv_query.size();
	for(i=0; i<q_size; i++){
		rank[cv_query[i]] = 1.0 / q_size;
	}

	/*if(q_size==0){
		printf("Bad!\n");
		getchar();
	}*/

	CSeekIntIntMap *mapG = sDataset.GetGeneMap();
	CSeekIntIntMap *mapQ = sDataset.GetQueryMap();

	CFullMatrix<float> *f = sDataset.GetFloatMatrix();

	size_t iGenesPresent = mapG->GetNumSet();
	for(i=0; i<iGenesPresent; i++){
		size_t g = mapG->GetReverse(i);
		for(j=0; j<q_size; j++){
			int qq = cv_query[j];
			if(g==qq) continue;
			size_t q = mapQ->GetForward(qq);
			/*if(f->Get(g,q)<-50.0 || f->Get(g,q)>50.0){
				printf("Bad %.5f\n", f->Get(g,q));
				getchar();
			}*/
			new_rank[g] += rank[qq] * f->Get(g, q);
		}
	}

	for(i=0; i<iGenesPresent; i++){
		size_t g = mapG->GetReverse(i);
		rank[g] = new_rank[g];
		//printf("Gene %d %.5f\n", g, rank[g]);
	}

	//getchar();

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

		/*printf("Cross Val %d %d %d\n", qi, num_q, num_v);
		printf("Cross Val %d\n", qi);
		for(i=0; i<vi.size(); i++){
			printf("%d ", vi[i]);
		}
		printf("\n");
		*/
		if(num_q==0 || num_v==0){
			sDataset.SetCVWeight(qi, -1);
		}else{
			/* actual weighting */
			vector<float> rank;
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
