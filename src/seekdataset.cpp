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
	geneMap = new CSeekIntIntMap(genePresence.size());
	ushort i;
	ushort iSize = genePresence.size();
	for(i=0; i<iSize; i++){
		if(genePresence[i]==1 && !isnan(geneAverage[i]) && !isinf(geneAverage[i])){
			geneMap->Add(i);
		}
	}
	return true;
}

/* requires presence vector */
bool CSeekDataset::InitializeQuery(const vector<char> &query){
	ushort iSize = query.size();
	ushort i, j;
	queryMap = new CSeekIntIntMap(iSize);
	ushort iGenes = geneMap->GetNumSet();
	for(i=0; i<iGenes; i++){
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
	r = CSeekTools::Init2DArray(iQuerySize, iNumGenes, (unsigned char) 255);
	/*r = new CFullMatrix<unsigned char>();
	r->Initialize(iQuerySize, iNumGenes);
	for(i=0; i<iQuerySize; i++){
		for(j=0; j<iNumGenes; j++){
			r->Set(i, j, 255);
		}
	}*/
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
	//r->Set(query, j, c);
	return true;
}

bool CSeekDataset::SetQueryNoMapping(const ushort &i, const ushort &j, const unsigned char &c){
	r[i][j] = c;
	//r->Set(i, j, c);
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
		//r->Set(query, j, c[j]);
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
	//rData->Initialize(r->GetColumns(), r->GetRows());
	ushort i, j;
	rData = rD;
	ushort ii;
	ushort iNumGenes = geneMap->GetNumSet();
	ushort iNumQueries = queryMap->GetNumSet();

	//iRows is the gene id, iColumns is the query id
	memset(&rData[0][0], 0, sizeof(ushort)*iRows*iColumns);
	/*for(i=0; i<iRows; i++){
		for(j=0; j<iColumns; j++){
			rData[i][j] = 0;
		}
	}
	printf("iRows and iColumns %d %d\n", iRows, iColumns);
	 */
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

			for(ii=0; ii<iNumGenes; ii++){
				i = geneMap->GetReverse(ii);
				/* numGenes */
				float a = GetGeneAverage(i);
				//if(isnan(a) || isinf(a)){
					/*for(j=0; j<rData->GetColumns(); j++){
						rData->Set(i, j, 0);
					}*/
				//	continue;
				//}
				/* numQueries */
				for(j=0; j<iNumQueries; j++){
					unsigned char x = r[j][i];
					if(x==255){
						//rData->Set(i, j, 0);
						continue;
					}
					/*if(CMeta::IsNaN(platform_avg[j]) ||
						CMeta::IsNaN(platform_stdev[j])){
						printf("platform average or stdev is NaN\n");
						getchar();
						continue;
					}*/
					float vv = (quant[x] - a - platform_avg[j]) / platform_stdev[j];
					vv = max(min(vv, (float)3.2), (float)-3.2);
					/*if(vv>=3.2){
						vv = 3.2;
					}
					else if(vv<=-3.2){
						vv = -3.2;
					}*/
					rData[i][j] = (ushort) (vv*100.0) + 320;
				}
			}
			delete[] platform_avg;
			delete[] platform_stdev;

		}else{
			for(ii=0; ii<iNumGenes; ii++){
				i = geneMap->GetReverse(ii);
				float a = GetGeneAverage(i);
				//if(isnan(a) || isinf(a)){
					/*for(j=0; j<rData->GetColumns(); j++){
						rData->Set(i, j, 0);
					}*/
				//	continue;
				//}
				/* numQueries */
				for(j=0; j<iNumQueries; j++){
					unsigned char x = r[j][i];
					if(x==255){
						//rData->Set(i, j, 0);
						continue;
					}
					float vv = quant[x] - a;
					if(vv>=3.2){
						vv = 3.2;
					}
					else if(vv<=-3.2){
						vv = -3.2;
					}
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
			if(vv>=3.2){
				vv = 3.2;
			} else if(vv<=-3.2){
				vv = -3.2;
			}
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
