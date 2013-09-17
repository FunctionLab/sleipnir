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
#include "seekreader.h"

namespace Sleipnir {

string CSeekTools::ConvertInt(const int &number){
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

bool CSeekTools::IsNaN(const utype &v){
	if(v==MAX_UTYPE) return true;
	return false;
}

utype CSeekTools::GetNaN(){
	return MAX_UTYPE;
}

bool CSeekTools::ReadDatabaselets(const vector<CDatabase*> &DB,
	const size_t &iGenes, const size_t &iDatasets,
	const vector< vector<string> > &vecstrAllQuery,
	vector<CSeekDataset*> &vc, const map<string,utype> &mapstriGenes,
	const vector<vector<string> > &dbDatasets,
	const map<string,utype> &mapstriDatasets,
	//network mode (data sent to client)
	const int &iClient, const bool &bNetwork){

	//requires LoadDatabase to be called beforehand
	size_t i, j, k;
	vector<char> cAllQuery;

	CSeekTools::InitVector(cAllQuery, iGenes, (char) 0);

	for(i=0; i<vecstrAllQuery.size(); i++){
		for(j=0; j<vecstrAllQuery[i].size(); j++){
			if(mapstriGenes.find(vecstrAllQuery[i][j])==mapstriGenes.end()) continue;
			utype k = mapstriGenes.find(vecstrAllQuery[i][j])->second;
			cAllQuery[k] = 1;
		}
	}

	vector<utype> allQ;
	for(i=0; i<cAllQuery.size(); i++) if(cAllQuery[i]==1) allQ.push_back(i);
	allQ.resize(allQ.size());

	//for now
	for(i=0; i<iDatasets; i++){
		if(vc[i]->GetDBMap()!=NULL){
			vc[i]->DeleteQueryBlock();
		}
	}

	int ret; //system call return

	fprintf(stderr, "Initializing query map\n"); 
	ret = system("date +%s%N 1>&2");
	if(bNetwork && CSeekNetwork::Send(iClient, "Initializing query map")==-1){
		fprintf(stderr, "Error sending client message\n");
		return false;
	}
	
	//fprintf(stderr, "Here\n");	
	#pragma omp parallel for \
	shared(allQ) private(i) schedule(dynamic)
	for(i=0; i<iDatasets; i++){
		vc[i]->InitializeQueryBlock(allQ);
	}

	fprintf(stderr, "Done initializing query map\n");
	ret = system("date +%s%N 1>&2");
	if(bNetwork && CSeekNetwork::Send(iClient, 
		"Done initializing query map")==-1){
		fprintf(stderr, "Error sending client message\n");
		return false;
	}

	fprintf(stderr, "Reading %lu query genes' correlations\n",
		allQ.size());
	ret = system("date +%s%N 1>&2");
	if(bNetwork && CSeekNetwork::Send(iClient, "Reading " + 
		CSeekTools::ConvertInt(allQ.size()) + 
		" query genes' correlations")==-1){
		fprintf(stderr, "Error sending client message\n");
		return false;
	}

	size_t m, d;

	for(i=0; i<allQ.size(); i++){
		m = allQ[i];
		for(d=0; d<DB.size(); d++){ //number of CDatabase collections
			vector<unsigned char> Qi;
			if(!DB[d]->GetGene(m, Qi)){
				cerr << "Gene does not exist" << endl;
				continue;
			}
			utype db;
			CSeekIntIntMap *qu = NULL;
			unsigned char **r = NULL;
			vector<utype> vecDatasetID;
			for(j=0; j<dbDatasets[d].size(); j++){
				utype qq = mapstriDatasets.find(dbDatasets[d][j])->second;
				vecDatasetID.push_back(qq);
			}
			#pragma omp parallel for \
			shared(Qi) private(j, k) \
			firstprivate(m, qu, r, db) schedule(dynamic)
			for(j=0; j<vecDatasetID.size(); j++){
				if((qu=vc[vecDatasetID[j]]->GetDBMap())==NULL) continue;
				if(CSeekTools::IsNaN(db = (qu->GetForward(m)))) continue;
				for(r = vc[vecDatasetID[j]]->GetMatrix(), k=0; k<iGenes; k++)
					r[db][k] = Qi[k*vecDatasetID.size()+j];
			}
			Qi.clear();
		}
	}

	fprintf(stderr, "Finished reading query genes' correlations\n");
	ret = system("date +%s%N 1>&2");
	if(bNetwork && CSeekNetwork::Send(iClient, 
		"Finished reading query genes' correlations")==-1){
		fprintf(stderr, "Error sending client message\n");
		return false;
	}

	return true;
}
	
bool CSeekTools::ReadQuantFile(const string &strFile, vector<float> &quant,
	const int lineSize){
	return CSeekTools::ReadQuantFile(strFile.c_str(), quant, lineSize);
}

bool CSeekTools::ReadQuantFile(const char *file, vector<float> &quant, const int lineSize){
	ifstream ifsm;
	ifsm.open(file);
	char acBuffer[lineSize];
	utype c_iBuffer = lineSize;
	vector<string> vecstrLine;

	ifsm.getline(acBuffer, c_iBuffer -1);
	//fprintf(stderr, "%s\n", acBuffer);
	CMeta::Tokenize( acBuffer, vecstrLine, " ", false);
	quant.clear();
	utype i;
	for(i=0; i<vecstrLine.size(); i++){
		quant.push_back(atof(vecstrLine[i].c_str()));
		//fprintf(stderr, "%.5f\n", atof(vecstrLine[i].c_str()));
	}
	quant.resize(quant.size());
	ifsm.close();
	return true;
}

bool CSeekTools::LoadDatabase(const vector<CDatabase*> &DB,
	const size_t &iGenes, const size_t &iDatasets,
	vector<CSeekDataset*> &vc, const vector<CSeekDataset*> &vc_src, 
	vector<CSeekPlatform> &vp, const vector<CSeekPlatform> &vp_src, 
	const vector<string> &vecstrDatasets,
	const map<string, string> &mapstrstrDatasetPlatform, 
	const map<string, utype> &mapstriPlatform){

	size_t i, j, k;

	vc.clear();
	vc.resize(iDatasets);
	vp.clear();
	vp.resize(vp_src.size());

	for(i=0; i<vp.size(); i++)
		vp[i].Copy(vp_src[i]);

	int ret; //system call returns

	fprintf(stderr, "Start reading average and presence files\n");
	ret = system("date +%s%N 1>&2");
	fprintf(stderr, "Done reading average and presence files\n");
	ret = system("date +%s%N 1>&2");

	fprintf(stderr, "Initializing gene map\n"); ret = system("date +%s%N 1>&2");
	#pragma omp parallel for \
	private(i) schedule(dynamic)
	for(i=0; i<iDatasets; i++){
		vc[i] = new CSeekDataset();
		vc[i]->Copy(vc_src[i]);
		string strFileStem = vecstrDatasets[i];
		string strPlatform =
			mapstrstrDatasetPlatform.find(strFileStem)->second;
		utype platform_id = mapstriPlatform.find(strPlatform)->second;
		vc[i]->SetPlatform(vp[platform_id]);
	}

	fprintf(stderr, "Done initializing gene map\n"); ret = system("date +%s%N 1>&2");
	return true;
}

bool CSeekTools::LoadDatabase(const vector<CDatabase*> &DB,
	const size_t &iGenes, const size_t &iDatasets,
	const vector<CSeekDBSetting*> &DBSetting,
	const vector<string> &vecstrDatasets,
	const map<string, string> &mapstrstrDatasetPlatform,
	const map<string, utype> &mapstriPlatform, vector<CSeekPlatform> &vp,
	vector<CSeekDataset*> &vc, const vector<vector<string> > &dbDataset,
	const map<string,utype> &mapstriDataset,
	const bool bVariance, const bool bCorrelation){

	size_t i, j, k;
	vc.clear();
	vc.resize(iDatasets);

	if(bCorrelation){
		for(i=0; i<DB.size(); i++){
			if(DBSetting[i]->GetValue("sinfo")=="NA"){
				fprintf(stderr, "sinfo parameter must be given.\n");
				return false;
			}
		}
	}

	if(bVariance){
		for(i=0; i<DB.size(); i++){
			if(DBSetting[i]->GetValue("gvar")=="NA"){
				fprintf(stderr, "gene variance parameter must be given.\n");
				return false;
			}
		}
	}

	int ret; //system call return

	fprintf(stderr, "Start reading average and presence files\n");
	ret = system("date +%s%N 1>&2");
	for(i=0; i<DB.size(); i++){
		const vector<string> &dset = dbDataset[i];
		string strPrepInputDirectory = DBSetting[i]->GetValue("prep");
		string strGvarInputDirectory = DBSetting[i]->GetValue("gvar");
		string strSinfoInputDirectory = DBSetting[i]->GetValue("sinfo");

		for(j=0; j<dset.size(); j++){
			utype d = mapstriDataset.find(dset[j])->second;
			vc[d] = new CSeekDataset();
			string strFileStem = dset[j];
			string strAvgPath = strPrepInputDirectory + "/" +
				strFileStem + ".gavg";
			string strPresencePath = strPrepInputDirectory + "/" +
				strFileStem + ".gpres";
			vc[d]->ReadGeneAverage(strAvgPath);
			vc[d]->ReadGenePresence(strPresencePath);
			if(bVariance){
				string strVariancePath = strGvarInputDirectory + "/" +
					strFileStem + ".gexpvar";
				vc[d]->ReadGeneVariance(strVariancePath);
			}
			if(bCorrelation){
				string strSinfoPath = strSinfoInputDirectory + "/" +
					strFileStem + ".sinfo";
				vc[d]->ReadDatasetAverageStdev(strSinfoPath);
			}
			string strPlatform =
				mapstrstrDatasetPlatform.find(strFileStem)->second;
			utype platform_id = mapstriPlatform.find(strPlatform)->second;
			vc[d]->SetPlatform(vp[platform_id]);
		}
	}

	fprintf(stderr, "Done reading average and presence files\n");
	ret = system("date +%s%N 1>&2");

	fprintf(stderr, "Initializing gene map\n"); ret = system("date +%s%N 1>&2");
	#pragma omp parallel for \
	private(i) schedule(dynamic)
	for(i=0; i<iDatasets; i++) vc[i]->InitializeGeneMap();

	fprintf(stderr, "Done initializing gene map\n"); ret = system("date +%s%N 1>&2");
	return true;
}

bool CSeekTools::ReadPlatforms(const string &strPlatformDirectory,
		vector<CSeekPlatform> &plat, vector<string> &vecstrPlatforms,
		map<string, utype> &mapstriPlatforms, const int lineSize){
	return CSeekTools::ReadPlatforms(strPlatformDirectory.c_str(), plat,
		vecstrPlatforms, mapstriPlatforms, lineSize);
}

bool CSeekTools::ReadPlatforms(const char *plat_dir,
		vector<CSeekPlatform> &plat, vector<string> &vecstrPlatforms,
		map<string, utype> &mapstriPlatforms, const int lineSize){

	string strPlatformDirectory = plat_dir;
	string strAvgFile = strPlatformDirectory + "/" +
		"all_platforms.gplatavg";
	string strStdevFile = strPlatformDirectory + "/" +
		"all_platforms.gplatstdev";
	string strPlatformOrderFile = strPlatformDirectory + "/" +
		"all_platforms.gplatorder";

	CFullMatrix<float> plat_avg;
	plat_avg.Open(strAvgFile.c_str());
	CFullMatrix<float> plat_stdev;
	plat_stdev.Open(strStdevFile.c_str());
	plat.clear();
	plat.resize(plat_avg.GetRows());
	utype i, j;

	vecstrPlatforms.clear();
	mapstriPlatforms.clear();
	ifstream ifsm;
	ifsm.open(strPlatformOrderFile.c_str());
	char acBuffer[lineSize];
	utype c_iBuffer = lineSize;
	i = 0;
	while(!ifsm.eof()){
		ifsm.getline(acBuffer, c_iBuffer -1);
		if(acBuffer[0]==0) break;
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

bool CSeekTools::ReadListTwoColumns(const string &strFile,
		vector<string> &vecstrList1, vector<string> &vecstrList2, const int lineSize){
	return CSeekTools::ReadListTwoColumns(strFile.c_str(),
		vecstrList1, vecstrList2, lineSize);
}

bool CSeekTools::ReadListTwoColumns(const char *file,
		vector<string> &vecstrList1, vector<string> &vecstrList2,
		const int lineSize){
	ifstream ifsm;
	ifsm.open(file);
	if(!ifsm.is_open()){
		fprintf(stderr, "Error opening file %s\n", file);
		return false;
	}
	char acBuffer[lineSize];
	utype c_iBuffer = lineSize;
	vecstrList1.clear();
	vecstrList2.clear();

	while(!ifsm.eof()){
		ifsm.getline(acBuffer, c_iBuffer -1);
		if(acBuffer[0]==0) break;
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

bool CSeekTools::ReadListOneColumn(const string &strFile,
	vector<string> &vecstrList, CSeekStrIntMap &mapstriList, const int lineSize){
	return CSeekTools::ReadListOneColumn(strFile.c_str(),
		vecstrList, mapstriList, lineSize);
}


bool CSeekTools::ReadListOneColumn(const char *file,
	vector<string> &vecstrList, CSeekStrIntMap &mapstriList, const int lineSize){
	ifstream ifsm;
	ifsm.open(file);
	if(!ifsm.is_open()){
		fprintf(stderr, "Error opening file %s\n", file);
		return false;
	}

	char acBuffer[lineSize];
	utype c_iBuffer = lineSize;
	vecstrList.clear();

	int i = 0;
	while(!ifsm.eof()){
		ifsm.getline(acBuffer, c_iBuffer -1);
		if(acBuffer[0]==0) break;
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

bool CSeekTools::ReadMultipleQueries(const string &strFile,
	vector< vector<string> > &qList, const int lineSize){
	return CSeekTools::ReadMultipleQueries(strFile.c_str(), qList, lineSize);
}

bool CSeekTools::ReadMultipleQueries(const char *file,
	vector< vector<string> > &qList, const int lineSize){
	qList.clear();
	FILE *infile;
	if((infile=fopen(file, "r"))==NULL){
		fprintf(stderr, "Error opening file %s\n", file);
		return false;
	}

	char *acBuffer;
	int MAX_CHAR_PER_LINE = lineSize;
	int lineLen = MAX_CHAR_PER_LINE;
	acBuffer = (char*)malloc(lineLen);
	while(fgets(acBuffer, lineLen, infile)!=NULL){
		while(strlen(acBuffer)==lineLen-1){
			int len = strlen(acBuffer);
			fseek(infile, -len, SEEK_CUR);
			lineLen+=MAX_CHAR_PER_LINE;
			acBuffer = (char*)realloc(acBuffer, lineLen);
			char *ret = fgets(acBuffer, lineLen, infile);
		}
	}
	rewind(infile);

	while(fgets(acBuffer, lineLen, infile)!=NULL){
		char *p = strtok(acBuffer, "\n");
		vector<string> tok;
		CMeta::Tokenize(p, tok, " ");
		qList.push_back(tok);
	}
	qList.resize(qList.size());
	free(acBuffer);

	fclose(infile);
	return true;
}

bool CSeekTools::ReadMultipleNotQueries(const char *file,
	vector<vector<vector<string> > > &qList, const int lineSize){
	qList.clear();
	FILE *infile;
	if((infile=fopen(file, "r"))==NULL){
		fprintf(stderr, "Error opening file %s\n", file);
		return false;
	}

	char *acBuffer;
	int MAX_CHAR_PER_LINE = lineSize;
	int lineLen = MAX_CHAR_PER_LINE;
	acBuffer = (char*)malloc(lineLen);
	while(fgets(acBuffer, lineLen, infile)!=NULL){
		while(strlen(acBuffer)==lineLen-1){
			int len = strlen(acBuffer);
			fseek(infile, -len, SEEK_CUR);
			lineLen+=MAX_CHAR_PER_LINE;
			acBuffer = (char*)realloc(acBuffer, lineLen);
			char *ret = fgets(acBuffer, lineLen, infile);
		}
	}
	rewind(infile);

	while(fgets(acBuffer, lineLen, infile)!=NULL){
		char *p = strtok(acBuffer, "\n");
		vector<vector<string> > aQ;
		vector<string> tok;
		CMeta::Tokenize(p, tok, "|");
		aQ.resize(tok.size());
		int i;
		for(i=0; i<tok.size(); i++){
			vector<string> tmp;
			CMeta::Tokenize(tok[i].c_str(), tmp, " ");
			aQ[i] = tmp;
		}
		qList.push_back(aQ);
	}
	qList.resize(qList.size());
	free(acBuffer);
	fclose(infile);
	return true;
}

bool CSeekTools::ReadMultiGeneOneLine(const string &strFile,
	vector<string> &list, const int lineSize){
	return CSeekTools::ReadMultiGeneOneLine(strFile.c_str(), list, lineSize);
}

bool CSeekTools::ReadMultiGeneOneLine(const char *file,
	vector<string> &list, const int lineSize){
	list.clear();
	ifstream ifsm;
	ifsm.open(file);
	if(!ifsm.is_open()){
		fprintf(stderr, "Error opening file %s\n", file);
		return false;
	}

	char *acBuffer= (char*)malloc(lineSize);
	int c_iBuffer = lineSize;
	int i = 0;

	ifsm.getline(acBuffer, c_iBuffer -1);
	acBuffer[c_iBuffer-1] = 0;
	vector<string> tok;
	CMeta::Tokenize(acBuffer, tok, " ");
	for(i = 0; i<tok.size(); i++){
		list.push_back(tok[i]);
	}

	list.resize(list.size());
	ifsm.close();
	free(acBuffer);
	return true;
}

bool CSeekTools::ReadListOneColumn(const string &strFile,
	vector<string> &vecstrList, const int lineSize){
	return CSeekTools::ReadListOneColumn(strFile.c_str(), vecstrList, lineSize);
}

bool CSeekTools::ReadListOneColumn(const char *file,
	vector<string> &vecstrList, const int lineSize){
	ifstream ifsm;
	ifsm.open(file);

	if(!ifsm.is_open()){
		fprintf(stderr, "Error opening file %s\n", file);
		return false;
	}
	char acBuffer[lineSize];
	utype c_iBuffer = lineSize;
	vecstrList.clear();
	int i = 0;
	while(!ifsm.eof()){
		ifsm.getline(acBuffer, c_iBuffer -1);
		if(acBuffer[0]==0) break;
		acBuffer[c_iBuffer-1] = 0;
		string line = acBuffer;
		vecstrList.push_back(line);
	}
	vecstrList.resize(vecstrList.size());
	ifsm.close();
	return true;
}

}
