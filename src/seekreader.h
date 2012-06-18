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
#ifndef SEEKREADER_H
#define SEEKREADER_H

#include "seekmap.h"
#include "stdafx.h"
#include "datapair.h"


namespace Sleipnir {

class CSeekTools{
public:
	template<class tType>
	static bool ReadArray(const char *fileName, vector<tType> &vData){
		FILE *f = fopen(fileName, "rb");
		if(f==NULL){
			cerr << "File not found" << endl;
			return false;
		}
		size_t iSize;
		int ret;
		ret = fread((char*) (&iSize), 1, sizeof(iSize), f);
		vData.clear();
		vData.resize(iSize);
		tType *m_Data = (tType*)malloc(iSize*sizeof(tType));
		ret = fread((char*)m_Data, 1, iSize*sizeof(tType), f);
		size_t i;
		for(i=0; i<iSize; i++){
			vData[i] = m_Data[i];
		}
		free(m_Data);
		fclose(f);
		return true;
	}

	template<class tType>
	static bool WriteArray(const char *fileName, vector<tType> &vData){
		FILE *f = fopen(fileName, "wb");
		if(f==NULL){
			cerr << "File not found" << endl;
			return false;
		}
		size_t i;
		tType *m_Data = (tType*)malloc(vData.size()*sizeof(tType));
		for(i=0; i<vData.size(); i++){
			m_Data[i] = vData[i];
		}
		size_t iSize = vData.size();
		fwrite((char*) (&iSize), 1, sizeof(iSize), f);
		fwrite((char*) (m_Data), 1, iSize*sizeof(tType), f);
		free(m_Data);
		fclose(f);
		return true;
	}

	template<class tType>
	static bool InitVector(vector<tType> &vData, size_t iSize, tType tValue){
		size_t i;
		vData.clear();
		vData.resize(iSize);
		for(i=0; i<iSize; i++){
			vData[i] = tValue;
		}
		return true;
	}

	static bool CreatePresenceVector(vector<int> &srcData, vector<char> &destData, size_t iSize){
		size_t i;
		destData.clear();
		destData.resize(iSize);
		for(i=0; i<iSize; i++){
			destData[i] = 0;
		}
		for(i=0; i<srcData.size(); i++){
			destData[srcData[i]] = 1;
		}
		return true;
	}

};

class CSeekDataset{
public:
	CSeekDataset(){
		r = NULL;
		geneAverage.clear();
		geneVariance.clear();
		genePresence.clear();
		m_fDsetAverage = CMeta::GetNaN();
		m_fDsetStdev = CMeta::GetNaN();
	}
	~CSeekDataset(){
		if(r!=NULL){
			delete r;
		}
		geneAverage.clear();
		geneVariance.clear();
		genePresence.clear();
	}
	bool ReadGeneAverage(const string &strFileName){
		return CSeekTools::ReadArray(strFileName.c_str(), geneAverage);
	}
	bool ReadGeneVariance(const string &strFileName){
		return CSeekTools::ReadArray(strFileName.c_str(), geneVariance);
	}
	bool ReadGenePresence(const string &strFileName){
		bool ret = CSeekTools::ReadArray(strFileName.c_str(), genePresence);
		if(!ret) return ret;
		geneMap = new CSeekIntIntMap(genePresence);
		return true;
	}

	/* requires presence vector */
	bool InitializeQuery(vector<char> &query){
		size_t iSize = query.size();
		size_t i, j;
		queryMap = new CSeekIntIntMap(iSize);
		for(i=0; i<geneMap->GetNumSet(); i++){
			size_t j = geneMap->GetReverse(i);
			if(query[j]==0) continue;
			queryMap->Add(j);
		}
		iQuerySize = queryMap->GetNumSet();
		iNumGenes = iSize;

		if(iQuerySize==0){
			cerr << "Dataset will be skipped" << endl;
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

	bool DeleteQuery(){
		if(queryMap!=NULL){
			delete queryMap;
		}
		iQuerySize = 0;
		iNumGenes = 0;
		if(r!=NULL){
			delete r;
		}
		return true;
	}

	bool SetQuery(size_t &i, size_t &j, unsigned char &c){
		size_t query = queryMap->GetForward(i);
		if(query==-1){
			return false;
		}
		r->Set(query, j, c);
		return true;
	}

	bool SetQueryNoMapping(size_t &i, size_t &j, unsigned char &c){
		r->Set(i, j, c);
		return true;
	}

	bool SetQuery(size_t &i, vector<unsigned char> &c){
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

	CFullMatrix<float> *GetFloatMatrix(){
		return rData;
	}

	bool InitializeFloatMatrix(bool bSubtractAvg = true){
		//hard coded quant file
		vector<float> quant;
		float w = -5.0;
		while(w<5.01){
			quant.push_back(w);
			w+=1.0;
		}
		quant.resize(quant.size());
		rData = new CFullMatrix<float>();
		rData->Initialize(r->GetRows(), r->GetColumns());
		size_t i,j;
		if(bSubtractAvg){
			for(i=0; i<rData->GetRows(); i++){
				for(j=0; j<rData->GetColumns(); j++){
					float a = GetGeneAverage(j);
					rData->Set(i, j, quant[r->Get(i, j)] - a);
				}
			}
		}else{
			for(i=0; i<rData->GetRows(); i++){
				for(j=0; j<rData->GetColumns(); j++){
					rData->Set(i, j, quant[r->Get(i, j)]);
				}
			}
		}
		return true;
	}

	bool FreeFloatMatrix(){
		delete rData;
		return true;
	}

	CFullMatrix<unsigned char> *GetMatrix(){
		return r;
	}

	CSeekIntIntMap* GetGeneMap(){
		return geneMap;
	}

	CSeekIntIntMap* GetQueryMap(){
		return queryMap;
	}

	float GetGeneVariance(size_t i){
		return geneVariance[i];
	}

	float GetGeneAverage(size_t i){
		return geneAverage[i];
	}


private:
	string strName;
	string strPlatform;
	CFullMatrix<unsigned char> *r;
	vector<float> geneAverage;
	vector<float> geneVariance;

	vector<char> genePresence;
	CSeekIntIntMap *geneMap;
	CSeekIntIntMap *queryMap;

	/* previously known as sinfo file */
	float m_fDsetAverage;
	float m_fDsetStdev;

	size_t iQuerySize;
	size_t iNumGenes;

	vector<float> weight;
	CFullMatrix<float> *rData;

	bool m_bIsNibble;
};



}
#endif
