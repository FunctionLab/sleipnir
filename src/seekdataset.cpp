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
		delete r;
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
	bool ret = CSeekTools::ReadArray(strFileName.c_str(), genePresence);
	if(!ret) return ret;
	geneMap = new CSeekIntIntMap(genePresence);
	return true;
}

/* requires presence vector */
bool CSeekDataset::InitializeQuery(vector<char> &query){
	size_t iSize = query.size();
	size_t i, j;
	queryMap = new CSeekIntIntMap(iSize);
	for(i=0; i<geneMap->GetNumSet(); i++){
		j = geneMap->GetReverse(i);
		if(query[j]==0) continue;
		queryMap->Add(j);
	}
	iQuerySize = queryMap->GetNumSet();
	iNumGenes = iSize;

	if(iQuerySize==0){
		//cerr << "Dataset will be skipped" << endl;
		r = NULL;
		return true;
	}
	r = new CFullMatrix<unsigned char>();
	r->Initialize(iQuerySize, iNumGenes);
	for(i=0; i<iQuerySize; i++){
		for(j=0; j<iNumGenes; j++){
			r->Set(i, j, 255);
		}
	}
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
		delete r;
		r = NULL;
	}
	return true;
}

bool CSeekDataset::SetQuery(size_t &i, size_t &j, unsigned char &c){
	size_t query = queryMap->GetForward(i);
	if(query==-1){
		return false;
	}
	r->Set(query, j, c);
	return true;
}

bool CSeekDataset::SetQueryNoMapping(size_t &i, size_t &j, unsigned char &c){
	r->Set(i, j, c);
	return true;
}

bool CSeekDataset::SetQuery(size_t &i, vector<unsigned char> &c){		
	size_t query = queryMap->GetForward(i);
	if(query==-1){
		return false;
	}
	size_t j = 0;
	for(j=0; j<c.size(); j++){
		r->Set(query, j, c[j]);
	}
	return true;
}

CFullMatrix<short>* CSeekDataset::GetDataMatrix(){
	return rData;
}

bool CSeekDataset::InitializeDataMatrix(bool bSubtractAvg,
	bool bSubtractPlatformAvg){
	/* assume platform is already set */

	//hard coded quant file
	vector<float> quant;
	float w = -5.0;
	while(w<5.01){
		quant.push_back(w);
		w+=0.1;
	}
	quant.resize(quant.size());

	//rData = new CFullMatrix<float>();
	rData = new CFullMatrix<short>();

	/* transpose */
	/* numGenes * numQueries */
	rData->Initialize(r->GetColumns(), r->GetRows());

	size_t i,j;
	if(bSubtractAvg){

		if(bSubtractPlatformAvg){
			float *platform_avg = new float[rData->GetColumns()];
			float *platform_stdev = new float[rData->GetColumns()];

			for(j=0; j<rData->GetColumns(); j++){
				size_t jj = queryMap->GetReverse(j);
				platform_avg[j] = platform->GetPlatformAvg(jj);
				platform_stdev[j] = platform->GetPlatformStdev(jj);
			}
			for(i=0; i<rData->GetRows(); i++){
				/* numGenes */
				float a = GetGeneAverage(i);
				if(CMeta::IsNaN(a)){
					for(j=0; j<rData->GetColumns(); j++){
						rData->Set(i, j, -32768);
					}
					continue;
				}
				/* numQueries */
				for(j=0; j<rData->GetColumns(); j++){
					unsigned char x = r->Get(j, i);
					if(x==255){
						rData->Set(i, j, -32768);
						continue;
					}
					/*if(CMeta::IsNaN(platform_avg[j]) ||
						CMeta::IsNaN(platform_stdev[j])){
						printf("platform average or stdev is NaN\n");
						getchar();
						continue;
					}*/
					float vv = (quant[x] - a - platform_avg[j]) / platform_stdev[j];
					rData->Set(i, j, (short)(vv*100.0));
				}
			}
			delete[] platform_avg;
			delete[] platform_stdev;

		}else{
			for(i=0; i<rData->GetRows(); i++){
				float a = GetGeneAverage(i);
				if(CMeta::IsNaN(a)){
					for(j=0; j<rData->GetColumns(); j++){
						rData->Set(i, j, -32768);
					}
					continue;
				}
				/* numQueries */
				for(j=0; j<rData->GetColumns(); j++){
					unsigned char x = r->Get(j, i);
					if(x==255){
						rData->Set(i, j, -32768);
						continue;
					}
					float v = quant[x] - a;
					rData->Set(i, j, (short)(v*100.0));
				}
			}

		}


		return true;
	}

	/* numGenes */
	for(i=0; i<rData->GetRows(); i++){
		/* numQueries */
		for(j=0; j<rData->GetColumns(); j++){
			rData->Set(i, j, (short) (quant[r->Get(j, i)] * 100.0));
		}
	}

	return true;
}

bool CSeekDataset::FreeDataMatrix(){
	delete rData;
	return true;
}

CFullMatrix<unsigned char>* CSeekDataset::GetMatrix(){
	return r;
}

CSeekIntIntMap* CSeekDataset::GetGeneMap(){
	return geneMap;
}

CSeekIntIntMap* CSeekDataset::GetQueryMap(){
	return queryMap;
}

float CSeekDataset::GetGeneVariance(size_t i){
	return geneVariance[i];
}

float CSeekDataset::GetGeneAverage(size_t i){
	return geneAverage[i];
}

size_t CSeekDataset::GetNumGenes(){
	return iNumGenes;
}

bool CSeekDataset::InitializeCVWeight(size_t i){
	weight.clear();
	weight.resize(i);
	return true;
}

bool CSeekDataset::SetCVWeight(size_t i, float f){
	weight[i] = f;
	return true;
}

float CSeekDataset::GetDatasetSumWeight(){	
	size_t i;
	size_t num = 0;
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

CSeekPlatform& CSeekDataset::GetPlatform(){
	return *platform;
}




}
