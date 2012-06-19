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
#ifndef SEEKWEIGHT_H
#define SEEKWEIGHT_H

#include "stdafx.h"
#include "seekreader.h"
#include "seekquery.h"

namespace Sleipnir {


class CSeekWeighter{
public:
	CSeekWeighter(){

	}
	~CSeekWeighter(){

	}
	static bool CVWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset){
		sDataset.InitializeFloatMatrix();
		size_t iFold = sQuery.GetNumFold();
		sDataset.InitializeCVWeight(iFold);

		int i, j, qi, qj;

		vector<char> is_query_cross, is_gold;
		CSeekTools::InitVector(is_query_cross, sDataset.GetNumGenes(), (char) 0);
		CSeekTools::InitVector(is_gold, sDataset.GetNumGenes(), (char) 0);

		for(qi=0; qi<iFold; qi++){
			vector<int> vi = sQuery.GetCVQuery(qi);
			CSeekIntIntMap *mapQ = sDataset.GetQueryMap();
			int num_q = 0;
			int num_v = 0;
			for(i=0; i<vi.size(); i++){
				if(mapQ->GetForward(vi[i])==-1) continue;
				is_query_cross[vi[i]] = 1;
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
				sDataset.SetCVWeight(qi, 0);
				continue;
			}

		}
		sDataset.FreeFloatMatrix();
		return true;
	}



};


}
#endif
