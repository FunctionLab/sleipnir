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
#ifndef SEEKEVALUATE_H
#define SEEKEVALUATE_H

#include "stdafx.h"
#include "seekmap.h"

namespace Sleipnir {

struct AResult{
	int i;
	//float f;
	short f;
	bool operator<(const AResult& val) const{
		if(f <= val.f){
			return false;
		}else{
			return true;
		}
	}
};

struct AResultFloat{
	int i;
	float f;
	bool operator<(const AResultFloat& val) const{
		if(f <= val.f){
			return false;
		}else{
			return true;
		}
	}
};

class CSeekPerformanceMeasure{
public:
	static bool SortRankVector(vector<short> &rank,
		CSeekIntIntMap &mapG, vector<AResult> &a);
	/* designed specifically for a CSeekDataset */
	/* mask: the query genes which are not included in RBP calcualtion */
	static bool RankBiasedPrecision(float rate, vector<short> &rank, float &rbp,
		vector<char> &mask, vector<char> &gold, CSeekIntIntMap &mapG);
};


}
#endif
