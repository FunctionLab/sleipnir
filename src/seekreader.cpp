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
#include "stdafx.h"
#include "datapair.h"
#include "seekdataset.h"
#include "seekreader.h"
#include "database.h"

namespace Sleipnir {

bool CSeekTools::CreatePresenceVector(vector<int> &srcData, vector<char> &destData, size_t iSize){
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

bool CSeekTools::LoadDatabase(CDatabase &DB, string &strInputDirectory, string &strPrepInputDirectory, 
	vector<char> &cQuery, vector<string> &vecstrQuery, vector<string> &vecstrDatasets, 
	vector<CSeekDataset*> &vc){
		
	DB.Open(strInputDirectory);
	size_t iDatasets = DB.GetDatasets();
	size_t iGenes = DB.GetGenes();
	size_t i, j,k;
	vc.clear();
	vc.resize(iDatasets);
	for(i=0; i<iDatasets; i++){
		vc[i] = new CSeekDataset();
		string strFileStem = vecstrDatasets[i];
		//string strFileStem = CMeta::Deextension(CMeta::Basename(vecstrDatasets[i].c_str()));
		string strAvgPath = strPrepInputDirectory + "/" + strFileStem + ".gavg";
		string strPresencePath = strPrepInputDirectory + "/" + strFileStem + ".gpres";
		vc[i]->ReadGeneAverage(strAvgPath);
		vc[i]->ReadGenePresence(strPresencePath);
	}

	CSeekTools::InitVector(cQuery, iGenes, (char) 0);

	for(i=0; i<vecstrQuery.size(); i++){
		k = DB.GetGene(vecstrQuery[i]);
		if(k==-1) continue;
		cQuery[k] = 1;
	}
	for(i=0; i<iDatasets; i++){
		vc[i]->InitializeQuery(cQuery);
	}

	vector<unsigned char> *Q =
		new vector<unsigned char>[vecstrQuery.size()];

	for(i=0; i<vecstrQuery.size(); i++){
		if(!DB.GetGene(vecstrQuery[i], Q[i])){
			cerr << "Gene does not exist" << endl;
		}
	}

	for(i=0; i<vecstrQuery.size(); i++){
		if(DB.GetGene(vecstrQuery[i])==-1){
			continue;
		}
		size_t m = DB.GetGene(vecstrQuery[i]);
		size_t l = 0;
		for(j=0; j<iDatasets; j++){
			CSeekIntIntMap *qu = vc[j]->GetQueryMap();
			size_t query = qu->GetForward(m);
			if(query==-1) continue;
			for(k=0; k<iGenes; k++){
				unsigned char c = Q[i][k*iDatasets + j];
				vc[j]->SetQueryNoMapping(query, k, c);
			}
		}
	}

	delete[] Q;

	return true;
}

}
