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

#include "seekbasic.h"
#include "seekreader.h"
#include "seekquery.h"
#include "seekevaluate.h"

namespace Sleipnir {

class CSeekWeighter{
public:
	/*cv_query must be present in sDataset */
	static bool LinearCombine(vector<ushort> &rank,
		const vector<ushort> &cv_query, CSeekDataset &sDataset,
		const ushort &);
	static bool CVWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset,
		const float &rate, const float &percent_required,
		vector<ushort> *rrank, const CSeekQuery *goldStd = NULL);
	static bool OrderStatisticsRankAggregation(const ushort&, const ushort&,
		ushort**, const vector<ushort> &, vector<float>&, const ushort&);
	static bool OrderStatisticsPreCompute();

};


}
#endif
