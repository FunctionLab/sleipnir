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
#ifndef SEEKDATASET_H
#define SEEKDATASET_H

#include "seekmap.h"
#include "stdafx.h"
#include "datapair.h"
#include "seekplatform.h"


namespace Sleipnir {

class CSeekDataset{
public:
	CSeekDataset();
	~CSeekDataset();
	bool ReadGeneAverage(const string &);
	bool ReadGeneVariance(const string &);
	bool ReadGenePresence(const string &);
	bool InitializeQuery(vector<char> &);
	bool DeleteQuery();
	bool SetQuery(size_t &, size_t &, unsigned char &);
	bool SetQueryNoMapping(size_t &, size_t &, unsigned char &);
	bool SetQuery(size_t &, vector<unsigned char> &);

	CFullMatrix<short> *GetDataMatrix();
	bool InitializeDataMatrix(bool=true, bool=true);
	bool FreeDataMatrix();

	CFullMatrix<unsigned char> *GetMatrix();
	CSeekIntIntMap* GetGeneMap();
	CSeekIntIntMap* GetQueryMap();
	float GetGeneVariance(size_t);
	float GetGeneAverage(size_t);
	size_t GetNumGenes();
	bool InitializeCVWeight(size_t);
	bool SetCVWeight(size_t, float);
	float GetDatasetSumWeight();
	void SetPlatform(CSeekPlatform &);
	CSeekPlatform& GetPlatform();

private:
	string strName;
	CSeekPlatform *platform;
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
	float sum_weight;
	//CFullMatrix<float> *rData;
	CFullMatrix<short> *rData;

	bool m_bIsNibble;

};



}
#endif
