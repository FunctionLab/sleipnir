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
#include "seekdataset.h"
#include "seekreader.h"
#include "seekevaluate.h"

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

//Copy CSeekDataset from a given object
bool CSeekDataset::Copy(CSeekDataset *src){
	//copies the following data if available: geneAverage, genePresence
	//geneVariance, dsetAverage, dsetStdev, geneMap, platform
	
	if(src->geneAverage.size()>0){
		//fprintf(stderr, "Great a!\n");
		geneAverage.resize(src->geneAverage.size());
		ushort i;
		for(i=0; i<src->geneAverage.size(); i++){
			geneAverage[i] = src->geneAverage[i];
		}
		//copy(src->geneAverage.begin(), src->geneAverage.end(), 
		//	geneAverage.begin());
	}
	if(src->genePresence.size()>0){
		//fprintf(stderr, "Great b!\n");
		genePresence.resize(src->genePresence.size());
		ushort i;
		for(i=0; i<src->genePresence.size(); i++){
			genePresence[i] = src->genePresence[i];
		}
		//copy(src->genePresence.begin(), src->genePresence.end(),
		//	genePresence.begin());
	}
	if(src->geneVariance.size()>0){
		geneVariance.resize(src->geneVariance.size());
		ushort i;
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
		iNumGenes = 0;
	}

	ushort iSize = genePresence.size();
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

	//fprintf(stderr, "0x1 %lu %d %d\n", CMeta::GetMemoryUsage(), iDBSize, iNumGenes);
	r = CSeekTools::Init2DArray(iDBSize, iNumGenes, (unsigned char) 255);
	//fprintf(stderr, "0x2 %lu\n", CMeta::GetMemoryUsage());
	//CSeekTools::Free2DArray((unsigned short**)r);
	//free(r[0]);
	//free(r);
	//fprintf(stderr, "0x3 %lu\n", CMeta::GetMemoryUsage());


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
		CSeekTools::Free2DArray(r);
		r = NULL;
	}
	iDBSize = 0;
	return true;
}

//test function only
bool CSeekDataset::DeleteR(){
	delete dbMap;
	CSeekTools::Free2DArray(r);
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
	const bool bSubtractPlatformAvg, const bool logit,
	const bool bCorrelation,
	const float cutoff, 
	const bool bRandom, gsl_rng *rand){
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

	if(bCorrelation && (bSubtractAvg || bSubtractPlatformAvg)){
		fprintf(stderr, "%s, %s\n", "correlation mode enabled",
			"please set subtract_avg, subtract_platform to false");
		return false;
	}

	ushort i, j;
	rData = rD;
	ushort ii;
	ushort iNumGenes = geneMap->GetNumSet();
	ushort iNumQueries = iQuerySize;

	//fprintf(stderr, "iNumQueries is %d\n", iNumQueries);
	//iRows is the gene id, iColumns is the query id
	//default value for all rData entries is 0
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
			}

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
						//By default, cutoff = -nan (i.e., always true)
						if(vv>cutoff){
							rData[i][j]= (ushort) (vv*100.0) + 320;
							//fprintf(stderr, "%.5f %.5f %.5f %.5f\n", quant[x], vv, a, platform_avg[j]);
						}else{
							rData[i][j] = 0;
						}
					}
				}
			}else{
				for(ii=0; ii<iGeneMapSize; ii++){
					for(j=0, i = allRGenes[ii], a=GetGeneAverage(i);
						j<iNumQueries; j++){
						if((x = r[queryIndex[j]][i])==255) continue;
						//fprintf(stderr, "Correlation %d %d %.5f %.5f %.5f %.5f\n", queryIndex[j], i, quant[x], a, platform_avg[j], platform_stdev[j]);
						vv = (quant[x] - a - platform_avg[j])
							/ platform_stdev[j];
						vv = max((float) min(vv, (float)3.2), (float)-3.2);

						if(vv>cutoff){
							rData[i][j]= (ushort) (vv*100.0) + 320;
							//fprintf(stderr, "r %.2f\n", quant[x]);
						}else{
							rData[i][j]= 0;
						}
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
						vv = log(quant[x]) - log((float)(1.0-quant[x])) - a;
						vv = max((float)-3.2, (float)min((float) 3.2, (float) vv));
						//fprintf(stderr, "%.5f %.5f %.5f\n", quant[x], vv, a);
						if(vv>cutoff){
							rData[i][j]= (ushort) (vv*100.0) + 320;
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
							rData[i][j]= (ushort) (vv*100.0) + 320;
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
	else if(logit){
		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				float vv = log(quant[x]) - log((float) 1.0 - quant[x]);
				vv = max((float) min(vv, (float)3.2), (float)-3.2);
				if(vv>cutoff){
					rData[i][j] = (ushort) (vv*100.0) + 320;
				}else{
					rData[i][j] = 0;
				}
			}
		}
	}else if(bCorrelation){

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
					rData[i][j] = (ushort) (vv*100.0) + 320;
				}else{
					rData[i][j] = (ushort) 0; 	//default value for 
												//not meeting cutoff
				}
			}
		}

	}else{ //raw z-scores
		for(ii=0; ii<iNumGenes; ii++){
			unsigned char x;
			for(i = geneMap->GetReverse(ii), j=0; j<iNumQueries; j++){
				if((x = r[queryIndex[j]][i])==255) continue;
				float vv = quant[x];
				vv = max((float) min(vv, (float)3.2), (float)-3.2);
				if(vv>cutoff){
					rData[i][j] = (ushort) (vv*100.0) + 320;
				}else{
					rData[i][j] = 0;
				}
			}
		}
	}

	if(bRandom){
		vector<vector<ushort> > allRandom;
		allRandom.resize(iNumQueries);
		for(i=0; i<iNumQueries; i++){
			allRandom[i] = vector<ushort>();
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

		ushort **a = CSeekTools::Init2DArray(iNumQueries, max_size, (ushort)0);

		for(i=0; i<iNumQueries; i++){
			for(j=0; j<allRandom[i].size(); j++){
				a[i][j] = allRandom[i][j];
			}
			gsl_ran_shuffle(rand, a[i], allRandom[i].size(), sizeof(ushort));
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

float CSeekDataset::GetCVWeight(const ushort &i){
	return weight[i];
}

const vector<float>& CSeekDataset::GetCVWeight() const{
	return weight;
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
