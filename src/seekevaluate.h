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
	ushort i;
	unsigned short f;
	bool operator<(const AResult& val) const{
		if(f <= val.f){
			return false;
		}else{
			return true;
		}
	}
};

struct Ascending
{
    bool operator()( const AResult& lx, const AResult& rx ) const {
        return lx.f <= rx.f;
    }
};



struct AResultFloat{
	ushort i;
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
	static bool SortRankVector(const vector<unsigned short> &rank,
		const CSeekIntIntMap &mapG, vector<AResult> &a, const bool bAllocate = true,
		const ushort top = 0);
	/* designed specifically for a CSeekDataset */
	/* mask: the query genes which are not included in RBP calcualtion */
	static bool RankBiasedPrecision(const float &rate,
		const vector<unsigned short> &rank, float &rbp,
		const vector<char> &mask, const vector<char> &gold,
		const CSeekIntIntMap &mapG,
		/* optional */
		const bool bAllocate = true,
		vector<AResult> *sing= NULL,
		const ushort top = 0);
};


}
#endif
