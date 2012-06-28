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
#include "seekevaluate.h"
#include "stdafx.h"
#include "seekmap.h"

namespace Sleipnir {

bool CSeekPerformanceMeasure::SortRankVector(vector<short> &rank,
	CSeekIntIntMap &mapG, vector<AResult> &a){
	a.clear();
	int numGenesD = mapG.GetNumSet();
	float old_target = -32769;
	float new_target = -32769;
	float prev_target = -32769;
	int prev_numNonZero = 0;
	int numNonZero = 0;
	int ii, i, jj;

	while(1){
		numNonZero = 0;
		for(ii=0; ii<numGenesD; ii++){
			i = mapG.GetReverse(ii);
			if(rank[i]<=old_target) continue;
			new_target += rank[i];
			numNonZero++;
		}

		//printf("Non Zero %d %d\n", numNonZero, new_target);

		/* 1000 is adjustable, this is the top number of items to sort */
		if(numNonZero==0 || numNonZero<1000){
			old_target = prev_target;
			numNonZero = prev_numNonZero;
			break;
		}

		new_target /= (float) numNonZero;

		if(new_target <= old_target){
			numNonZero = prev_numNonZero;
			break;
		}
		prev_target = old_target;
		old_target = new_target;
		prev_numNonZero = numNonZero;
	}

	if(numNonZero==0){
		cerr << "This dataset is all zero!" << endl;
		return false;
	}

	a.resize(numNonZero);
	jj = 0;
	for(ii=0; ii<numGenesD; ii++){
		i = mapG.GetReverse(ii);
		if(rank[i]<=old_target) continue;
		a[jj].i = i;
		a[jj].f = rank[i];
		jj++;
	}
	sort(a.begin(), a.end());
	return true;
}

/* designed specifically for a CSeekDataset */
/* mask: the query genes which are not included in RBP calcualtion */
bool CSeekPerformanceMeasure::RankBiasedPrecision(float rate, vector<short> &rank, float &rbp,
	vector<char> &mask, vector<char> &gold, CSeekIntIntMap &mapG){

	int i, ii, j, jj;
	vector<AResult> sing;
	bool ret = CSeekPerformanceMeasure::SortRankVector(rank, mapG, sing);
	if(!ret){
		rbp = -1;
		return false;
	}

	float x = 0;
	jj = 0;
	int numNonZero = sing.size();
	for(i=0; i<numNonZero; i++){
		if(sing[i].f<=-32768) break;
		if(mask[sing[i].i]==1) continue;
		if(gold[sing[i].i]==1){
			x+=pow(rate, jj);
		}
		jj++;
	}
	x*=(1.0-rate);
	rbp = x;
	return true;
}


}
