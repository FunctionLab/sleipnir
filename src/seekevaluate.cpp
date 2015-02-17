/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "seekevaluate.h"
namespace Sleipnir {

bool CSeekPerformanceMeasure::SortRankVector(
	const vector<utype> &rank,
	const CSeekIntIntMap &mapG, vector<AResult> &a,
	const bool bNegativeCor,
	//optional argument
	const utype top){

	utype numGenesD = mapG.GetNumSet();
	utype TOP = 0;
	utype numNonZero = 0;
	utype i;
	bool DEBUG = false;

	//a should be the same size as rank
	if(top==0) TOP = rank.size();
	else TOP = top;

	vector<utype>::const_iterator itRank = rank.begin();
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
	if(bNegativeCor){
		if(TOP==rank.size()){
			sort(a.begin(), a.end(), Ascending());
		}else{
			nth_element(a.begin(), a.begin()+TOP, a.end(), Ascending());
			sort(a.begin(), a.begin()+TOP, Ascending());
		}

	}else{
		if(TOP==rank.size()){
			sort(a.begin(), a.end());
		}else{
			nth_element(a.begin(), a.begin()+TOP, a.end());
			sort(a.begin(), a.begin()+TOP);
		}
	}

	return true;
}

/* designed specifically for a CSeekDataset */
/* mask: the query genes which are not included in RBP calcualtion */
bool CSeekPerformanceMeasure::RankBiasedPrecision(const float &rate,
	const vector<utype> &rank, float &rbp,
	const vector<char> &mask, const vector<char> &gold,
	const CSeekIntIntMap &mapG, vector<AResult> *ar,
	const bool bNegativeCor,
	/* optional arguments */
	const utype top){

	utype i, ii, j, jj;
	float x;

	utype TOP = top;
	if(top==0) TOP = rank.size();

	//ar should be same size as rank
	vector<AResult> *sing = ar;
	bool ret =
		CSeekPerformanceMeasure::SortRankVector(rank, mapG, *sing, bNegativeCor, top);

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
		if(gold[aa->i]==1){
			x+=pow(rate, (float)jj);
			//fprintf(stderr, "Sorted %d %d %.5f\n", jj, aa->i, (aa->f-320)/100.0);
		}
		jj++;
	}
	x *= (1.0-rate);
	rbp = x;
	//fprintf(stderr, "%.3e\n", 0, rbp);

	return true;
}

bool CSeekPerformanceMeasure::AveragePrecision(
	const vector<utype> &rank, float &ap,
	const vector<char> &mask, const vector<char> &gold,
	const CSeekIntIntMap &mapG, vector<AResult> *ar, const bool bNegativeCor){

	utype i, ii, j, jj;
	float x;

	utype TOP = rank.size();

	//ar should be same size as rank
	vector<AResult> *sing = ar;
	bool ret =
		CSeekPerformanceMeasure::SortRankVector(rank, mapG, *sing, bNegativeCor, TOP);

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
