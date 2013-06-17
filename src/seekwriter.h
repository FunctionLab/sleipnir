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
#ifndef SEEKWRITER_H
#define SEEKWRITER_H

#include "seekbasic.h"
#include "seekmap.h"
#include "seekevaluate.h"
#include "datapair.h"

namespace Sleipnir {

class CSeekWriter{
public:
	static bool ReadSparseMatrix(const char *fileName, 
		vector<vector<float> > &mat, 
		CSeekIntIntMap &m, const int maxRank, const float rbp_p,
		const vector<string> &vecstrGenes);

	static bool ProductNorm(const vector<vector<float> > &mat1,
		const vector<vector<float> > &mat2, const CSeekIntIntMap &m1, 
		const CSeekIntIntMap &m2, vector<vector<float> > &re);

	static bool WriteSparseMatrix(vector<vector<unsigned short> > &umat,
		int maxRank, const vector<string> &vecstrGenes, const char *fileName);

	static bool GetSparseRankMatrix(CDat &Dat,
		vector<vector<unsigned short> > &umat, const unsigned short nullValue,
		int maxRank, const vector<string> &vecstrGenes);

	static bool RankNormalizeDAB(CDat &Dat,
		const vector<string> &vecstrGenes, int max_rank, float rbp_p);

	static bool NormalizeDAB(CDat &Dat,
		const vector<string> &vecstrGenes,
		bool cutoff, bool expTransform, bool divideNorm, bool subtractNorm);

	static bool GetGeneAverage(CDataPair &Dat,
		const vector<string> &vecstrGenes,
		vector<float> &vecResult, bool logit=false, float top_percent=1.0);
	static bool GetGenePresence(CDataPair &Dat,
		const vector<string> &vecstrGenes,
		vector<char> &vecResult);
	static bool GetDatasetSinfo(CDataPair &Dat, float &mean,
		float &stdev);
};

}
#endif
