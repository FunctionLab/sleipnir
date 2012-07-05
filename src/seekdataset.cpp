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
#include "seekmap.h"
#include "seekreader.h"
#include "seekdataset.h"
#include "stdafx.h"
#include "datapair.h"

namespace Sleipnir {

CSeekDataset::CSeekDataset(){
	r = NULL;
	rData = NULL;
	geneAverage.clear();
	geneVariance.clear();
	genePresence.clear();
	m_fDsetAverage = CMeta::GetNaN();
	m_fDsetStdev = CMeta::GetNaN();
	weight.clear();
	sum_weight = -1;
}

CSeekDataset::~CSeekDataset(){
	if(r!=NULL){
		CSeekTools::Free2DArray(r);
	}
	if(rData!=NULL){
		CSeekTools::Free2DArray(rData);
	}
	geneAverage.clear();
	geneVariance.clear();
	genePresence.clear();
}

bool CSeekDataset::ReadGeneAverage(const string &strFileName){
	return CSeekTools::ReadArray(strFileName.c_str(), geneAverage);
}

bool CSeekDataset::ReadGeneVariance(const string &strFileName){
	return CSeekTools::ReadArray(strFileName.c_str(), geneVariance);
}

bool CSeekDataset::ReadGenePresence(const string &strFileName){
	return CSeekTools::ReadArray(strFileName.c_str(), genePresence);
}

bool CSeekDataset::InitializeGeneMap(){
	if(geneAverage.empty() || genePresence.empty()) {
		cerr << "Gene average or gene presence unread" << endl;
		return false;
	}
	ushort i;
	ushort iSize = genePresence.size();
	geneMap = new CSeekIntIntMap(iSize);
	vector<char>::const_iterator iterGenePresence = genePresence.begin();
	vector<float>::const_iterator iterGeneAverage = geneAverage.begin();

	for(i=0; iterGenePresence!=genePresence.end(); i++, iterGenePresence++, iterGeneAverage++){
		if(*iterGenePresence==1 && !CMeta::IsNaN(*iterGeneAverage)){
			geneMap->Add(i);
		}
	}
	return true;
}

/* requires presence vector */
bool CSeekDataset::InitializeQuery(const vector<char> &query){
	ushort iSize = query.size();
	queryMap = new CSeekIntIntMap(iSize);
	ushort iGenes = geneMap->GetNumSet();
	const vector<ushort> &allRGenes = geneMap->GetAllReverse();
	vector<ushort>::const_iterator iterR = allRGenes.begin();
	vector<ushort>::const_iterator endR = allRGenes.begin() + iGenes;
	for(; iterR!=endR; iterR++){
		if(query[*iterR]==0) continue;
		queryMap->Add(*iterR);
	}
	iQuerySize = queryMap->GetNumSet();
	iNumGenes = iSize;
	bool DEBUG = false;

	if(iQuerySize==0){
		if(DEBUG) cerr << "Dataset will be skipped" << endl;
		r = NULL;
		return true;
	}
	r = CSeekTools::Init2DArray(iQuerySize, iNumGenes, (unsigned char) 255);
	return true;
}

bool CSeekDataset::DeleteQuery(){
	if(queryMap!=NULL){
		delete queryMap;
		queryMap = NULL;
	}
	iQuerySize = 0;
	iNumGenes = 0;
	if(r!=NULL){
		CSeekTools::Free2DArray(r);
		r = NULL;
	}
	return true;
}


bool CSeekDataset::SetQuery(const ushort &i, const ushort &j, const unsigned char &c){
	ushort query = queryMap->GetForward(i);
	if(CSeekTools::IsNaN(query)){
		return false;
	}
	r[query][j] = c;
	return true;
}

bool CSeekDataset::SetQueryNoMapping(const ushort &i, const ushort &j, const unsigned char &c){
	r[i][j] = c;
	return true;
}

bool CSeekDataset::SetQuery(const ushort &i, const vector<unsigned char> &c){
	ushort query = queryMap->GetForward(i);
	if(CSeekTools::IsNaN(query)){
		return false;
	}
	ushort j = 0;
	for(j=0; j<c.size(); j++){
		r[query][j] = c[j];
	}
	return true;
}

ushort** CSeekDataset::GetDataMatrix(){
	return rData;
}

bool CSeekDataset::InitializeDataMatrix(ushort **rD, const ushort &iRows, const ushort &iColumns,
	const bool bSubtractAvg, const bool bSubtractPlatformAvg){
	/* assume platform is already set */

	//hard coded quant file
	vector<float> quant;
	float w = -5.0;
	while(w<5.01){
		quant.push_back(w);
		w+=0.1;
	}
	quant.resize(quant.size());

	//Assuming rData is not NULL
	/* transpose */
	/* numGenes * numQueries */
	ushort i, j;
	rData = rD;
	ushort ii;
	ushort iNumGenes = geneMap->GetNumSet();
	ushort iNumQueries = queryMap->GetNumSet();

	//iRows is the gene id, iColumns is the query id
	memset(&rData[0][0], 0, sizeof(ushort)*iRows*iColumns);

	if(bSubtractAvg){

		if(bSubtractPlatformAvg){
			float *platform_avg = new float[iColumns];
			float *platform_stdev = new float[iColumns];

			//GetColumns() is numQuery
			for(j=0; j<queryMap->GetNumSet(); j++){
				ushort jj = queryMap->GetReverse(j);
				platform_avg[j] = platform->GetPlatformAvg(jj);
				platform_stdev[j] = platform->GetPlatformStdev(jj);
			}

			const vector<ushort> &allRGenes = geneMap->GetAllReverse();
			vector<ushort>::const_iterator iterRGenes = allRGenes.begin();
			vector<ushort>::const_iterator lastRGenes = allRGenes.begin() + geneMap->GetNumSet();
			for(; iterRGenes != lastRGenes; iterRGenes++){
				i = *iterRGenes;
				/* numGenes */
				float a = GetGeneAverage(i);
				/* numQueries */
				ushort *rDataP = &rData[i][0];
				unsigned char *rP = &r[0][i];
				float *plAvgP = &platform_avg[0];
				float *plStdevP = &platform_stdev[0];
				ushort *end = &rData[i][0] + iNumQueries;

				for(; rDataP!=end; plAvgP++, plStdevP++, rDataP++, rP+=iRows){
					if(*rP==255){
						continue;
					}
					float vv = (quant[*rP] - a - *plAvgP) / *plStdevP;
					//Do not remove the (float) cast in front of min
					vv = max((float) min(vv, (float)3.2), (float)-3.2);
					*rDataP = (ushort) (vv*100.0) + 320;
				}
			}
			delete[] platform_avg;
			delete[] platform_stdev;

		}else{
			for(ii=0; ii<iNumGenes; ii++){
				i = geneMap->GetReverse(ii);
				float a = GetGeneAverage(i);
				/* numQueries */
				for(j=0; j<iNumQueries; j++){
					unsigned char x = r[j][i];
					if(x==255){
						continue;
					}
					float vv = quant[x] - a;
					vv = max((float) min(vv, (float)3.2), (float)-3.2);
					rData[i][j]= (ushort) (vv*100.0) + 320;
				}
			}
		}

		return true;
	}

	/* numGenes */
	for(ii=0; ii<iNumGenes; ii++){
		i = geneMap->GetReverse(ii);
		/* numQueries */
		for(j=0; j<iNumQueries; j++){
			float vv = quant[r[j][i]];
			vv = max((float) min(vv, (float)3.2), (float)-3.2);
			rData[i][j] = (ushort) (vv*100.0) + 320;
		}
	}

	return true;
}


unsigned char** CSeekDataset::GetMatrix(){
	return r;
}

CSeekIntIntMap* CSeekDataset::GetGeneMap(){
	return geneMap;
}

CSeekIntIntMap* CSeekDataset::GetQueryMap(){
	return queryMap;
}

float CSeekDataset::GetGeneVariance(const ushort &i) const{
	return geneVariance[i];
}

float CSeekDataset::GetGeneAverage(const ushort &i) const{
	return geneAverage[i];
}

ushort CSeekDataset::GetNumGenes() const{
	return iNumGenes;
}

bool CSeekDataset::InitializeCVWeight(const ushort &i){
	weight.clear();
	weight.resize(i);
	return true;
}

bool CSeekDataset::SetCVWeight(const ushort &i, const float &f){
	weight[i] = f;
	return true;
}

float CSeekDataset::GetDatasetSumWeight(){
	ushort i;
	ushort num = 0;
	if(sum_weight==-1){
		sum_weight = 0;
		for(i=0; i<weight.size(); i++){
			if(weight[i]==-1) continue;
			sum_weight+=weight[i];
			num++;
		}
		if(num>0){
			sum_weight/=(float)num;
		}else{
			sum_weight = -1;
		}
	}
	return sum_weight;
}

void CSeekDataset::SetPlatform(CSeekPlatform &cp){
	platform = &cp;
}

CSeekPlatform& CSeekDataset::GetPlatform() const{
	return *platform;
}




}
