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
namespace Sleipnir {

bool CSeekPerformanceMeasure::SortRankVector(
	const vector<unsigned short> &rank,
	const CSeekIntIntMap &mapG, vector<AResult> &a,
	//optional argument
	const ushort top){

	ushort numGenesD = mapG.GetNumSet();
	ushort TOP = 0;
	ushort numNonZero = 0;
	ushort i;
	bool DEBUG = false;

	//a should be the same size as rank
	if(top==0) TOP = rank.size();
	else TOP = top;

	vector<ushort>::const_iterator itRank = rank.begin();
	vector<AResult>::iterator itA = a.begin();
	for(i = 0; itRank!=rank.end(); itRank++, itA++, i++){
		itA->i = i;
		itA->f = *itRank;
		if(*itRank>0) numNonZero++;
	}

	if(numNonZero==0){
		if(DEBUG) cerr << "This dataset is all zero!" << endl;
		return false;
	}

	//printf("Top is %d", TOP); getchar();
	if(TOP==rank.size()){
		sort(a.begin(), a.end());
	}else{
		nth_element(a.begin(), a.begin()+TOP, a.end());
		sort(a.begin(), a.begin()+TOP);
	}

	return true;
}

/* designed specifically for a CSeekDataset */
/* mask: the query genes which are not included in RBP calcualtion */
bool CSeekPerformanceMeasure::RankBiasedPrecision(const float &rate,
	const vector<unsigned short> &rank, float &rbp,
	const vector<char> &mask, const vector<char> &gold,
	const CSeekIntIntMap &mapG, vector<AResult> *ar,
	/* optional arguments */
	const ushort top){

	ushort i, ii, j, jj;
	float x;

	ushort TOP = top;
	if(top==0) TOP = rank.size();

	//ar should be same size as rank
	vector<AResult> *sing = ar;
	bool ret =
		CSeekPerformanceMeasure::SortRankVector(rank, mapG, *sing, top);

	if(!ret){
		rbp = -1;
		return false;
	}

	x = 0;
	jj = 0;

	AResult *aa;
	for(i=0; i<TOP; i++){
		aa = &(*sing)[i];
		if(aa->f==0) break;
		if(mask[aa->i]==1) continue;
		if(gold[aa->i]==1) x+=pow(rate, jj);
		jj++;
	}
	x *= (1.0-rate);
	rbp = x;
	return true;
}

bool CSeekPerformanceMeasure::AveragePrecision(
	const vector<unsigned short> &rank, float &ap,
	const vector<char> &mask, const vector<char> &gold,
	const CSeekIntIntMap &mapG, vector<AResult> *ar){

	ushort i, ii, j, jj;
	float x;

	ushort TOP = rank.size();

	//ar should be same size as rank
	vector<AResult> *sing = ar;
	bool ret =
		CSeekPerformanceMeasure::SortRankVector(rank, mapG, *sing, TOP);

	if(!ret){
		ap = -1;
		return false;
	}

	jj = 0;

	AResult *aa;
	int num = 0;
	float sum = 0;
	for(i=0; i<TOP; i++){
		aa = &(*sing)[i];
		if(aa->f==0) break;
		if(mask[aa->i]==1) continue;
		if(gold[aa->i]==1){
			sum+=(float) (num+1) / (float) (jj+1);
			num++;
		}
		jj++;
	}
	if(num==0){
		ap = 0;
	}else{
		ap = sum / num;
	}
	return true;
}

}
