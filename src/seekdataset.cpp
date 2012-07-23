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
#include "seekevaluate.h"
#include "stdafx.h"
#include "datapair.h"

namespace Sleipnir {

CSeekDataset::CSeekDataset(){
	r = NULL;
	rData = NULL;
	dbMap = NULL;
	geneMap = NULL;
	queryMap = NULL;
	platform = NULL;
	geneAverage.clear();
	geneVariance.clear();
	genePresence.clear();
	query.clear();
	queryIndex.clear();
	iQuerySize = 0;
	iNumGenes = 0;
	iDBSize = 0;
	m_bIsNibble = false;
	m_fDsetAverage = CMeta::GetNaN();
	m_fDsetStdev = CMeta::GetNaN();
	weight.clear();
	sum_weight = -1;
}

CSeekDataset::~CSeekDataset(){
	DeleteQuery();
	DeleteQueryBlock();
	if(geneMap!=NULL){
		delete geneMap;
		geneMap = NULL;
		iNumGenes = 0;
	}
	geneAverage.clear();
	geneVariance.clear();
	genePresence.clear();
	m_fDsetAverage = CMeta::GetNaN();
	m_fDsetStdev = CMeta::GetNaN();
	platform = NULL;
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
	iNumGenes = iSize;
	geneMap = new CSeekIntIntMap(iSize);
	vector<char>::const_iterator iterGenePresence = genePresence.begin();
	vector<float>::const_iterator iterGeneAverage = geneAverage.begin();

	for(i=0; iterGenePresence!=genePresence.end(); i++,
		iterGenePresence++, iterGeneAverage++){
		if(*iterGenePresence==1 && !CMeta::IsNaN(*iterGeneAverage)){
			geneMap->Add(i);
		}
	}
	return true;
}

/* requires presence vector */
bool CSeekDataset::InitializeQueryBlock(const vector<ushort> &queryBlock){
	DeleteQueryBlock();

	dbMap = new CSeekIntIntMap(iNumGenes);

	vector<ushort>::const_iterator iterQ = queryBlock.begin();
	for(; iterQ!=queryBlock.end(); iterQ++){
		if(CSeekTools::IsNaN(geneMap->GetForward(*iterQ))){
			//this query gene is not present in the dataset
			continue;
		}
		dbMap->Add(*iterQ);
	}
	iDBSize = dbMap->GetNumSet();

	bool DEBUG = false;

	if(iDBSize==0){
		if(DEBUG) cerr << "Dataset will be skipped" << endl;
		DeleteQueryBlock();
		return true;
	}

	r = CSeekTools::Init2DArray(iDBSize, iNumGenes, (unsigned char) 255);

	return true;
}

bool CSeekDataset::InitializeQuery(const vector<ushort> &query){
	DeleteQuery();

	if(iDBSize==0 || dbMap==NULL) return true;

	//must require initializequeryBlock be executed first

	queryMap = new CSeekIntIntMap(iNumGenes);

	ushort i;
	vector<ushort>::const_iterator iterQ = query.begin();

	vector<AResult> a;
	a.resize(query.size());
	vector<AResult>::iterator iterA = a.begin();
	for(iQuerySize = 0; iterQ!=query.end(); iterQ++){
		if(CSeekTools::IsNaN(i = dbMap->GetForward(*iterQ))) continue;
		(*iterA).i = *iterQ;
		(*iterA).f = i;
		iterA++;
		iQuerySize++;
	}

	bool DEBUG = false;
	if(iQuerySize==0){
		if(DEBUG) cerr << "Dataset will be skipped" << endl;
		DeleteQuery();
		return true;
	}

	a.resize(iQuerySize);
	sort(a.begin(), a.end(), Ascending());

	for(iterA = a.begin(); iterA!=a.end(); iterA++){
		queryMap->Add((*iterA).i);
		this->query.push_back((*iterA).i);
		this->queryIndex.push_back((*iterA).f);
	}

	this->queryIndex.resize(this->queryIndex.size());
	this->query.resize(this->query.size());

	return true;
}

bool CSeekDataset::DeleteQuery(){
	if(queryMap!=NULL){
		delete queryMap;
		queryMap = NULL;
	}
	iQuerySize = 0;
	rData = NULL;
	weight.clear();
	query.clear();
	queryIndex.clear();
	sum_weight = -1;
	return true;
}

bool CSeekDataset::DeleteQueryBlock(){
	DeleteQuery();
	if(dbMap!=NULL){
		delete dbMap;
		dbMap = NULL;
	}
	if(r!=NULL){
		CSeekTools::Free2DArray(r);
		r = NULL;
	}
	iDBSize = 0;
	return true;
}

const vector<ushort>& CSeekDataset::GetQuery() const{
	return this->query;
}

const vector<ushort>& CSeekDataset::GetQueryIndex() const{
	return this->queryIndex;
}

ushort** CSeekDataset::GetDataMatrix(){
	return rData;
}

bool CSeekDataset::InitializeDataMatrix(ushort **rD,
	const vector<float> &quant, const ushort &iRows,
	const ushort &iColumns, const bool bSubtractAvg,
	const bool bSubtractPlatformAvg, const bool logit){
	/* assume platform is already set */

	//hard coded quant file
	/*vector<float> quant;
	float w = -5.0;
	while(w<5.01){
		quant.push_back(w);
		w+=0.04;
		//w+=0.1;
	}
	quant.resize(quant.size());
	 */
	//Assuming rData is not NULL
	/* transpose */
	/* numGenes * numQueries */
	ushort i, j;
	rData = rD;
	ushort ii;
	ushort iNumGenes = geneMap->GetNumSet();
	ushort iNumQueries = iQuerySize;

	//fprintf(stderr, "iNumQueries is %d\n", iNumQueries);
	//iRows is the gene id, iColumns is the query id
	memset(&rData[0][0], 0, sizeof(ushort)*iRows*iColumns);

	//assume queryIndex is already sorted
	vector<ushort> offset;
	offset.push_back(0);
	for(i=1; i<queryIndex.size(); i++)
		offset.push_back(queryIndex[i] - queryIndex[i-1]);

	if(bSubtractAvg){

		if(bSubtractPlatformAvg){
			float *platform_avg = new float[iColumns];
			float *platform_stdev = new float[iColumns];

			//GetColumns() is numQuery
			for(j=0; j<iNumQueries; j++){
				ushort jj = this->query[j];
				platform_avg[j] = platform->GetPlatformAvg(jj);
				platform_stdev[j] = platform->GetPlatformStdev(jj);
				//fprintf(stderr, "Q%d ", jj);
			}
			//fprintf(stderr, "\n");

			const vector<ushort> &allRGenes = geneMap->GetAllReverse();
			float a = 0;
			float vv = 0;
			unsigned char x = 0;
			ushort iGeneMapSize = geneMap->GetNumSet();
			if(logit){
				for(ii=0; ii<iGeneMapSize; ii++){
					for(j=0, i = allRGenes[ii], a=GetGeneAverage(i);
						j<iNumQueries; j++){
						if((x = r[queryIndex[j]][i])==255) continue;
						vv = ((log(quant[x]) - log((float) 1.0 - quant[x]))
							- a - platform_avg[j]) / platform_stdev[j];
						vv = max((float) min(vv, (float)3.2), (float)-3.2);
						//fprintf(stderr, "%.5f %.5f %.5f %.5f\n", quant[x], vv, a, platform_avg[j]);
						rData[i][j]= (ushort) (vv*100.0) + 320;
					}
				}
			}else{
				for(ii=0; ii<iGeneMapSize; ii++){
					for(j=0, i = allRGenes[ii], a=GetGeneAverage(i);
						j<iNumQueries; j++){
						if((x = r[queryIndex[j]][i])==255) continue;
						vv = (quant[x] - a - platform_avg[j])
							/ platform_stdev[j];
						vv = max((float) min(vv, (float)3.2), (float)-3.2);
						rData[i][j]= (ushort) (vv*100.0) + 320;
						//fprintf(stderr, "r %.2f\n", quant[x]);
					}
				}
			}

			/*const vector<ushort> &allRGenes = geneMap->GetAllReverse();
			vector<ushort>::const_iterator iterRGenes = allRGenes.begin();
			vector<ushort>::const_iterator lastRGenes = allRGenes.begin() + geneMap->GetNumSet();

			for(; iterRGenes != lastRGenes; iterRGenes++){
				i = *iterRGenes;
				float a = GetGeneAverage(i);

				// numQueries
				ushort *rDataP = &rData[i][0];
				unsigned char *rP = &r[queryIndex[0]][i];
				float *plAvgP = &platform_avg[0];
				float *plStdevP = &platform_stdev[0];
				ushort *end = &rData[i][0] + iNumQueries;
				vector<ushort>::const_iterator iterOffset = offset.begin();
				for(; rDataP!=end; plAvgP++, plStdevP++, rDataP++, iterOffset++, rP+=iRows*(*iterOffset)){
					if(*rP==255){
						continue;
					}
					float vv = (quant[*rP] - a - *plAvgP) / *plStdevP;
					//fprintf(stderr, "%.3f\n", vv);
					//Do not remove the (float) cast in front of min
					vv = max((float) min(vv, (float)3.2), (float)-3.2);
					*rDataP = (ushort) (vv*100.0) + 320;
				}
			}*/

			delete[] platform_avg;
			delete[] platform_stdev;

		}else{
			if(logit){
				for(ii=0; ii<iNumGenes; ii++){
					float a, vv; unsigned char x;
					for(i = geneMap->GetReverse(ii),
						a = GetGeneAverage(i), j=0; j<iNumQueries; j++){
						if((x = r[queryIndex[j]][i])==255) continue;
						vv = max((float)-3.2, (float)min((float) 3.2, (float)
							(log(quant[x]) - log((float)(1.0 - quant[x]))
							- a)));
						//fprintf(stderr, "%.5f %.5f %.5f\n", quant[x], vv, a);
						rData[i][j]= (ushort) (vv*100.0) + 320;
					}
				}
			}
			else{
				for(ii=0; ii<iNumGenes; ii++){
					float a, vv; unsigned char x;
					/* numQueries */
					for(i = geneMap->GetReverse(ii),
						a = GetGeneAverage(i), j=0; j<iNumQueries; j++){
						if((x = r[queryIndex[j]][i])==255) continue;
						vv = max((float) min((float)(quant[x] - a),
							(float)3.2), (float)-3.2);
						rData[i][j]= (ushort) (vv*100.0) + 320;
					}
				}
			}
		}

		return true;
	}

	/* numGenes */
	if(logit){
		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				float vv = log(quant[x]) - log((float) 1.0 - quant[x]);
				vv = max((float) min(vv, (float)3.2), (float)-3.2);
				rData[i][j] = (ushort) (vv*100.0) + 320;
			}
		}
	}else{
		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				float vv = quant[x];
				vv = max((float) min(vv, (float)3.2), (float)-3.2);
				rData[i][j] = (ushort) (vv*100.0) + 320;
			}
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

CSeekIntIntMap* CSeekDataset::GetDBMap(){
	return dbMap;
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
	sum_weight = -1;
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
		for(sum_weight=0, i=0; i<weight.size(); i++){
			if(weight[i]==-1) continue;
			sum_weight+=weight[i];
			num++;
		}
		if(num>0) sum_weight/=(float)num;
		else sum_weight = -1;
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
