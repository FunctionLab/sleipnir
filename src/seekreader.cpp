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

bool CSeekTools::LoadDatabase(CDatabase &DB, string &strPrepInputDirectory,
	vector<char> &cQuery, vector<string> &vecstrQuery, vector<string> &vecstrDatasets, 
	map<string, string> &mapstrstrDatasetPlatform, map<string, size_t> &mapstriPlatform,
	vector<CSeekPlatform> &vp, vector<CSeekDataset*> &vc){
		
	size_t iDatasets = DB.GetDatasets();
	size_t iGenes = DB.GetGenes();
	size_t i, j, k;
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
		string strPlatform = mapstrstrDatasetPlatform[strFileStem];
		size_t platform_id = mapstriPlatform[strPlatform];
		//printf("Platform id %s %d\n", strPlatform.c_str(), platform_id);
		vc[i]->SetPlatform(vp[platform_id]);
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

bool CSeekTools::ReadPlatforms(string &strPlatformDirectory, vector<CSeekPlatform> &plat,
		vector<string> &vecstrPlatforms, map<string, size_t> &mapstriPlatforms){

	string strAvgFile = strPlatformDirectory + "/" + "all_platforms.gplatavg";
	string strStdevFile = strPlatformDirectory + "/" + "all_platforms.gplatstdev";
	string strPlatformOrderFile = strPlatformDirectory + "/" + "all_platforms.gplatorder";

	CFullMatrix<float> plat_avg;
	plat_avg.Open(strAvgFile.c_str());
	CFullMatrix<float> plat_stdev;
	plat_stdev.Open(strStdevFile.c_str());
	plat.clear();
	plat.resize(plat_avg.GetRows());
	size_t i, j;

	/*for(i=0; i<plat_avg.GetRows(); i++){
		int c = 0;
		for(j=0; j<plat_avg.GetColumns(); j++){
			if(CMeta::IsNaN(plat_avg.Get(i,j))){
				continue;
			}
			c++;
			//printf("Gene %d %d: %.5f %.5f\n", i, j, plat_avg.Get(i, j), plat_stdev.Get(i,j));
		}
		printf("Platform %d, %d\n", i, c);
	}
	printf("Done");
	 */
	vecstrPlatforms.clear();
	mapstriPlatforms.clear();
	ifstream ifsm;
	ifsm.open(strPlatformOrderFile.c_str());
	char acBuffer[1024];
	int c_iBuffer = 1024;
	i = 0;
	while(!ifsm.eof()){
		ifsm.getline(acBuffer, c_iBuffer -1);
		if(acBuffer[0]==0){
			break;
		}
		acBuffer[c_iBuffer-1] = 0;
		vecstrPlatforms.push_back(acBuffer);
		mapstriPlatforms[acBuffer] = i;
		i++;
	}
	vecstrPlatforms.resize(vecstrPlatforms.size());
	ifsm.close();

	for(i=0; i<plat_avg.GetRows(); i++){
		plat[i].InitializePlatform(plat_avg.GetColumns(), vecstrPlatforms[i]);
		for(j=0; j<plat_avg.GetColumns(); j++){
			plat[i].SetPlatformAvg(j, plat_avg.Get(i, j));
			plat[i].SetPlatformStdev(j, plat_stdev.Get(i, j));
		}
	}

	return true;
}

bool CSeekTools::ReadListTwoColumns(string &strFile, vector<string> &vecstrList1, vector<string> &vecstrList2){
	ifstream ifsm;
	ifsm.open(strFile.c_str());
	if(!ifsm.is_open()){
		cerr << "Error opening file " << strFile << endl;
		return false;
	}
	char acBuffer[1024];
	int c_iBuffer = 1024;
	vecstrList1.clear();
	vecstrList2.clear();

	while(!ifsm.eof()){
		ifsm.getline(acBuffer, c_iBuffer -1);
		if(acBuffer[0]==0){
			break;
		}
		acBuffer[c_iBuffer-1] = 0;
		vector<string> tok;
		CMeta::Tokenize(acBuffer, tok);
		vecstrList1.push_back(tok[0]);
		vecstrList2.push_back(tok[1]);
	}
	vecstrList1.resize(vecstrList1.size());
	vecstrList2.resize(vecstrList2.size());
	ifsm.close();
	return true;
}

bool CSeekTools::ReadListOneColumn(string &strFile, vector<string> &vecstrList, CSeekStrIntMap &mapstriList){
	ifstream ifsm;
	ifsm.open(strFile.c_str());
	if(!ifsm.is_open()){
		cerr << "Error opening file " << strFile << endl;
		return false;
	}

	char acBuffer[1024];
	int c_iBuffer = 1024;
	vecstrList.clear();

	size_t i = 0;
	while(!ifsm.eof()){
		ifsm.getline(acBuffer, c_iBuffer -1);
		if(acBuffer[0]==0){
			break;
		}
		acBuffer[c_iBuffer-1] = 0;
		string line = acBuffer;
		vecstrList.push_back(line);
		mapstriList.Set(line, i);
		i++;
	}
	vecstrList.resize(vecstrList.size());
	ifsm.close();
	return true;
}


}
