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
#include "seekdataset.h"
#include "seekreader.h"
#include "seekevaluate.h"

namespace Sleipnir {

CSeekDataset::CSeekDataset(){
	r = NULL; //if declared, will be iDBSize * iNumGenes
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
	//fprintf(stderr, "Destructor called!\n");
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

//Copy CSeekDataset from a given object
bool CSeekDataset::Copy(CSeekDataset *src){
	//copies the following data if available: geneAverage, genePresence
	//geneVariance, dsetAverage, dsetStdev, geneMap, platform
	
	if(src->geneAverage.size()>0){
		//fprintf(stderr, "Great a!\n");
		geneAverage.resize(src->geneAverage.size());
		utype i;
		for(i=0; i<src->geneAverage.size(); i++){
			geneAverage[i] = src->geneAverage[i];
		}
	}
	if(src->genePresence.size()>0){
		//fprintf(stderr, "Great b!\n");
		genePresence.resize(src->genePresence.size());
		utype i;
		for(i=0; i<src->genePresence.size(); i++){
			genePresence[i] = src->genePresence[i];
		}
	}
	if(src->geneVariance.size()>0){
		geneVariance.resize(src->geneVariance.size());
		utype i;
		for(i=0; i<src->geneVariance.size(); i++){
			geneVariance[i] = src->geneVariance[i];
		}
		//copy(src->geneVariance.begin(), src->geneVariance.end(),
		//	geneVariance.begin());
	}
	m_fDsetAverage = src->m_fDsetAverage;
	m_fDsetStdev = src->m_fDsetStdev;

	if(geneMap!=NULL){
		delete geneMap;
		geneMap = NULL;
		iNumGenes = 0;
	}

	utype iSize = genePresence.size();
	iNumGenes = iSize;
	geneMap = new CSeekIntIntMap(src->geneMap);

	return true;	
}

bool CSeekDataset::ReadDatasetAverageStdev(const string &strFileName){
	vector<float> t;
	CSeekTools::ReadArray(strFileName.c_str(), t);
	m_fDsetAverage = t[0];
	m_fDsetStdev = t[1];
	return true;
}

bool CSeekDataset::ReadGeneAverage(const string &strFileName){
	return CSeekTools::ReadArray(strFileName.c_str(), geneAverage);
}

bool CSeekDataset::ReadGenePresence(const string &strFileName){
	return CSeekTools::ReadArray(strFileName.c_str(), genePresence);
}

bool CSeekDataset::ReadGeneVariance(const string &strFileName){
	return CSeekTools::ReadArray(strFileName.c_str(), geneVariance);
}


bool CSeekDataset::InitializeGeneMap(){
	if(geneAverage.empty() || genePresence.empty()) {
		cerr << "Gene average or gene presence unread" << endl;
		return false;
	}
	utype i;
	utype iSize = genePresence.size();
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


bool CSeekDataset::InitializeQueryBlock(const vector<utype> &queryBlock){

	DeleteQueryBlock();

	dbMap = new CSeekIntIntMap(iNumGenes);

	vector<utype>::const_iterator iterQ = queryBlock.begin();
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

	//fprintf(stderr, "0x1 %lu %d %d\n", CMeta::GetMemoryUsage(), iDBSize, iNumGenes);
	r = CSeekTools::Init2DArray(iDBSize, iNumGenes, (unsigned char) 255);
	//fprintf(stderr, "0x2 %lu\n", CMeta::GetMemoryUsage());
	//free(r[0]);
	//free(r);
	//fprintf(stderr, "0x3 %lu\n", CMeta::GetMemoryUsage());


	return true;
}


bool CSeekDataset::InitializeQuery(const vector<utype> &query){
	DeleteQuery();

	if(iDBSize==0 || dbMap==NULL) return true;
	//must require initializequeryBlock be executed first

	queryMap = new CSeekIntIntMap(iNumGenes);

	utype i;
	vector<utype>::const_iterator iterQ = query.begin();

	vector<AResult> a;
	a.resize(query.size());
	vector<AResult>::iterator iterA = a.begin();

	for(iQuerySize = 0; iterQ!=query.end(); iterQ++){
		//if the query does not exist in this data-set, continue
		if(CSeekTools::IsNaN(i = dbMap->GetForward(*iterQ))) continue;
		//.i: query genes
		(*iterA).i = *iterQ;
		//.f: query gene position in dbMap
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
		//add gene to queryMap
		queryMap->Add((*iterA).i);
		this->query.push_back((*iterA).i);
		//add (gene-position in dbMap) to queryIndex
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
		//fprintf(stderr, "iDBSize is %d\n", iDBSize);
		CSeekTools::Free2DArray(r);
		r = NULL;
	}
	iDBSize = 0;
	return true;
}

const vector<utype>& CSeekDataset::GetQuery() const{
	return this->query;
}

const vector<utype>& CSeekDataset::GetQueryIndex() const{
	return this->queryIndex;
}

utype** CSeekDataset::GetDataMatrix(){
	return rData;
}

bool CSeekDataset::InitializeDataMatrix(utype **rD,
	const vector<float> &quant, const utype &iRows,
	const utype &iColumns, const bool bSubtractAvg,
	const bool bNormPlatform, const bool logit,
	const enum CSeekDataset::DistanceMeasure dist_measure,
	const float cutoff,
	const bool bRandom, gsl_rng *rand){
	/* assume platform is already set */
	//Assuming rData is not NULL

	bool bCorrelation = false;
	if(dist_measure==CSeekDataset::CORRELATION)
		bCorrelation = true;

	if(bCorrelation && (bSubtractAvg || bNormPlatform)){
		fprintf(stderr, "%s, %s\n", "correlation mode enabled",
			"please set subtract_avg, subtract_platform to false");
		return false;
	}

	utype i, j;
	rData = rD;
	utype ii;
	utype iNumGenes = geneMap->GetNumSet();
	utype iNumQueries = iQuerySize;

	//fprintf(stderr, "iNumQueries is %d\n", iNumQueries);
	//iRows is the gene id, iColumns is the query id
	//default value for all rData entries is 0
	memset(&rData[0][0], 0, sizeof(utype)*iRows*iColumns);

	//assume queryIndex is already sorted
	vector<utype> offset;
	offset.push_back(0);
	for(i=1; i<queryIndex.size(); i++)
		offset.push_back(queryIndex[i] - queryIndex[i-1]);

	if(bSubtractAvg){ //z-score (Z) + average z score subtraction (A)

		if(bNormPlatform){ //subtract platform avg / platform stddev
			float *platform_avg = new float[iColumns];
			float *platform_stdev = new float[iColumns];

			//GetColumns() is numQuery
			for(j=0; j<iNumQueries; j++){
				utype jj = this->query[j];
				platform_avg[j] = platform->GetPlatformAvg(jj);
				platform_stdev[j] = platform->GetPlatformStdev(jj);
			}

			const vector<utype> &allRGenes = geneMap->GetAllReverse();
			float a = 0;
			float vv = 0;
			unsigned char x = 0;
			utype iGeneMapSize = geneMap->GetNumSet();
			if(logit){ //NOT checked
				for(ii=0; ii<iGeneMapSize; ii++){
					for(j=0, i = allRGenes[ii], a=GetGeneAverage(i);
						j<iNumQueries; j++){
						if((x = r[queryIndex[j]][i])==255) continue;
						vv = ((log(quant[x]) - log((float) 1.0 - quant[x]))
							- a - platform_avg[j]) / platform_stdev[j];
						vv = max((float) min(vv, (float)3.2), (float)-3.2);
						//By default, cutoff = -nan (i.e., always true)
						if(vv>cutoff){
							rData[i][j]= (utype) (vv*100.0) + 320;
							//fprintf(stderr, "%.5f %.5f %.5f %.5f\n", quant[x], vv, a, platform_avg[j]);
						}else{
							rData[i][j] = 0;
						}
					}
				}
			}else{ //if just normal score
				for(ii=0; ii<iGeneMapSize; ii++){
					for(j=0, i = allRGenes[ii], a=GetGeneAverage(i);
						j<iNumQueries; j++){
						if((x = r[queryIndex[j]][i])==255) continue;
						vv = (quant[x] - a - platform_avg[j])
							/ platform_stdev[j];
						vv = max((float) min(vv, (float)3.2), (float)-3.2);
						if(vv>cutoff){
							rData[i][j]= (utype) (vv*100.0) + 320;
							//fprintf(stderr, "r %.2f\n", quant[x]);
						}else{
							rData[i][j]= 0;
						}
					}
				}
			}

			delete[] platform_avg;
			delete[] platform_stdev;

		}else{ //do not normalize by platform, just subtract gene avg 
			if(logit){
				for(ii=0; ii<iNumGenes; ii++){
					float a, vv; unsigned char x;
					for(i = geneMap->GetReverse(ii),
						a = GetGeneAverage(i), j=0; j<iNumQueries; j++){
						if((x = r[queryIndex[j]][i])==255) continue;
						vv = log(quant[x]) - log((float)(1.0-quant[x])) - a;
						vv = max((float)-3.2, (float)min((float) 3.2, (float) vv));
						//fprintf(stderr, "%.5f %.5f %.5f\n", quant[x], vv, a);
						if(vv>cutoff){
							rData[i][j]= (utype) (vv*100.0) + 320;
						}else{
							rData[i][j] = 0;
						}
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
						vv = quant[x] - a;
						vv = max((float) min((float)vv, (float)3.2), (float)-3.2);
						if(vv>cutoff){
							rData[i][j]= (utype) (vv*100.0) + 320;
						}else{
							rData[i][j] = 0;
						}
					}
				}
			}
		}

		//return true;
	}
	/* numGenes */
	else if(logit){ //logit on z-scores
		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				float vv = log(quant[x]) - log((float) 1.0 - quant[x]);
				vv = max((float) min(vv, (float)3.2), (float)-3.2);
				if(vv>cutoff){
					rData[i][j] = (utype) (vv*100.0) + 320;
				}else{
					rData[i][j] = 0;
				}
			}
		}
	}else if(bCorrelation){ //correlation mode
		//assumed to be from -1 to +1

		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				float vv = quant[x];
				vv = vv * m_fDsetStdev + m_fDsetAverage;
				float e = exp(2.0*vv);
				vv = (e - 1.0) / (e + 1.0);
				if(vv>cutoff){
					vv = vv * 3.0; //scale up by factor of 3 for retaining 
								   //precision, should put value to range 
								   //(-3.0 to 3.0)
					vv = max((float) min(vv, (float)3.2), (float)-3.2);
					rData[i][j] = (utype) (vv*100.0) + 320;
				}else{
					rData[i][j] = (utype) 0; 	//default value for 
												//not meeting cutoff
				}
			}
		}

	}else{ //just simple z-scores, no subtraction of avg z scores
		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				float vv = quant[x];

				//for functional network=================================
				//vv = vv * 6.0 - 3.0; //transform values from 0-1 to values from -3 to +3

				vv = max((float) min(vv, (float)3.2), (float)-3.2);
				if(vv>cutoff){
					rData[i][j] = (utype) (vv*100.0) + 320;
				}else{
					rData[i][j] = 0;
				}
			}
		}
	}

	if(bRandom){
		vector<vector<utype> > allRandom;
		allRandom.resize(iNumQueries);
		for(i=0; i<iNumQueries; i++){
			allRandom[i] = vector<utype>();
		}

		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				allRandom[j].push_back(rData[i][j]);
			}
		}

		int max_size = 0;
		for(i=0; i<iNumQueries; i++){
			if(allRandom[i].size()>max_size){
				max_size = allRandom[i].size();
			}
		}

		utype **a = CSeekTools::Init2DArray(iNumQueries, max_size, (utype)0);

		for(i=0; i<iNumQueries; i++){
			for(j=0; j<allRandom[i].size(); j++){
				a[i][j] = allRandom[i][j];
			}
			gsl_ran_shuffle(rand, a[i], allRandom[i].size(), sizeof(utype));
		}

		vector<int> k;
		CSeekTools::InitVector(k, iNumQueries, (int) 0);

		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				//fprintf(stderr, "%d %d\n", rData[i][j], a[j][k[j]]);
				rData[i][j] = a[j][k[j]];
				k[j]++;
			}
		}
		
		CSeekTools::Free2DArray(a);
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

float CSeekDataset::GetDatasetAverage() const{
	return m_fDsetAverage;
}

float CSeekDataset::GetDatasetStdev() const{
	return m_fDsetStdev;
}

float CSeekDataset::GetGeneVariance(const utype &i) const{
	return geneVariance[i];
}

float CSeekDataset::GetGeneAverage(const utype &i) const{
	return geneAverage[i];
}

utype CSeekDataset::GetNumGenes() const{
	return iNumGenes;
}

bool CSeekDataset::InitializeCVWeight(const utype &i){
	weight.clear();
	weight.resize(i);
	sum_weight = -1;
	return true;
}

bool CSeekDataset::SetCVWeight(const utype &i, const float &f){
	weight[i] = f;
	return true;
}

float CSeekDataset::GetCVWeight(const utype &i){
	return weight[i];
}

const vector<float>& CSeekDataset::GetCVWeight() const{
	return weight;
}

float CSeekDataset::GetDatasetSumWeight(){
	utype i;
	utype num = 0;
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
