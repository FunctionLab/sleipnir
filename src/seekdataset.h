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
	bool InitializeQuery(const vector<ushort> &);
	bool InitializeQueryBlock(const vector<ushort> &);
	bool DeleteQuery();
	bool DeleteQueryBlock();

	/*bool SetQuery(const ushort &, const ushort &, const unsigned char &);
	bool SetQueryNoMapping(const ushort &, const ushort &, const unsigned char &);
	bool SetQuery(const ushort &, const vector<unsigned char> &);
	 */
	bool InitializeDataMatrix(ushort**, const ushort&, const ushort&,
		const bool=true, const bool=true);
	ushort** GetDataMatrix();

	unsigned char** GetMatrix();
	CSeekIntIntMap* GetGeneMap();
	CSeekIntIntMap* GetDBMap();
	CSeekIntIntMap* GetQueryMap();

	const vector<ushort>& GetQuery() const;
	const vector<ushort>& GetQueryIndex() const;

	float GetGeneVariance(const ushort&) const;
	float GetGeneAverage(const ushort&) const;
	ushort GetNumGenes() const;
	bool InitializeCVWeight(const ushort&);
	bool SetCVWeight(const ushort&, const float&);
	float GetDatasetSumWeight();
	void SetPlatform(CSeekPlatform &);
	CSeekPlatform& GetPlatform() const;
	bool InitializeGeneMap();


private:
	string strName;
	CSeekPlatform *platform;
	vector<float> geneAverage;
	vector<float> geneVariance;

	vector<char> genePresence;
	CSeekIntIntMap *geneMap;
	CSeekIntIntMap *dbMap;
	CSeekIntIntMap *queryMap;
	vector<ushort> query;
	vector<ushort> queryIndex;

	/* previously known as sinfo file */
	float m_fDsetAverage;
	float m_fDsetStdev;

	ushort iQuerySize;
	ushort iNumGenes;
	ushort iDBSize;

	vector<float> weight;

	ushort **rData;
	unsigned char **r;

	float sum_weight;
	bool m_bIsNibble;

};



}
#endif
