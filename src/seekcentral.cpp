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
#include "seekcentral.h"

namespace Sleipnir {

CSeekCentral::CSeekCentral(){
	m_vecstrGenes.clear();
	m_vecstrDatasets.clear();
	m_vecstrSearchDatasets.clear();
	m_mapstrstrDatasetPlatform.clear();
	m_vc.clear();
	m_quant.clear();
	m_vecstrAllQuery.clear();
	m_vp.clear();
	m_mapstriPlatform.clear();
	m_vecstrPlatform.clear();
	m_vecstrDP.clear();
	m_mapstrintDataset.clear();
	m_mapstrintGene.clear();
	m_searchdsetMap.clear();
	m_DB = NULL;
	m_rData = NULL;
	m_maxNumDB = 50;

	m_master_rank_threads = NULL;
	m_sum_weight_threads = NULL;
	m_counts_threads = NULL;
	m_rank_normal_threads = NULL;
	m_rank_threads = NULL;
	m_rank_d = NULL;

	m_master_rank.clear();
	m_sum_weight.clear();
	m_weight.clear();
	m_counts.clear();
	m_mapLoadTime.clear();

	m_Query.clear();
	m_final.clear();

	m_iDatasets = 0;
	m_iGenes = 0;
	m_numThreads = 0;
	m_bSubtractGeneAvg = false;
	m_bSubtractPlatformAvg = false;
	m_bDividePlatformStdev = false;
	m_bLogit = false;
	m_bCorrelation = false;
	m_bOutputText = false;
	m_bSquareZ = false;
	m_bSharedDB = false;

	DEBUG = false;
	m_output_dir = "";

	m_iClient = -1;
	m_bEnableNetwork = false;
	m_bNetworkSendData = false;
	m_bNetworkSendStatus = false;
}

CSeekCentral::~CSeekCentral(){
	m_vecstrGenes.clear();
	m_vecstrDatasets.clear();
	m_vecstrSearchDatasets.clear();
	m_mapstrstrDatasetPlatform.clear();
	m_vecstrDP.clear();
	m_mapstrintDataset.clear();
	m_mapstrintGene.clear();

	ushort i, j;
	for(i=0; i<m_vc.size(); i++){
		if(m_vc[i]==NULL) continue;
		delete m_vc[i];
	}
	m_vc.clear();
	
	m_quant.clear();

	for(i=0; i<m_searchdsetMap.size(); i++){
		if(m_searchdsetMap[i]==NULL) continue;
		delete m_searchdsetMap[i];
		m_searchdsetMap[i] = NULL;
	}
	m_searchdsetMap.clear();

	//m_rData, m_master_rank_threads,
	//m_sum_weight_threads, m_counts_threads
	//m_rank_normal_threads, and m_rank_threads
	//should be already freed
	if(m_master_rank_threads!=NULL){
		CSeekTools::Free2DArray(m_master_rank_threads);
		m_master_rank_threads = NULL;		
	}
	if(m_sum_weight_threads != NULL){
		CSeekTools::Free2DArray(m_sum_weight_threads);
		m_sum_weight_threads = NULL;
	}
	if(m_counts_threads !=NULL){
		CSeekTools::Free2DArray(m_counts_threads);
		m_counts_threads = NULL;
	}
	if(m_rank_normal_threads!=NULL){
		for(j=0; j<m_numThreads; j++)
			m_rank_normal_threads[j].clear();
		delete[] m_rank_normal_threads;
		m_rank_normal_threads = NULL;
	}
	if(m_rank_threads !=NULL){
		for(j=0; j<m_numThreads; j++)
			m_rank_threads[j].clear();
		delete[] m_rank_threads;
		m_rank_threads = NULL;
	}
	if(m_rank_d!=NULL){
		CSeekTools::Free2DArray(m_rank_d);
		m_rank_d = NULL;
	}

	m_master_rank.clear();
	m_sum_weight.clear();
	m_counts.clear();
	m_weight.clear();
	m_final.clear();

	m_vecstrAllQuery.clear();
	m_Query.clear();

	m_vp.clear();

	m_mapstriPlatform.clear();
	m_vecstrPlatform.clear();

	if(m_DB!=NULL){
		if(!m_bSharedDB){
			delete m_DB;
		}
		m_DB = NULL;
	}
	m_iDatasets = 0;
	m_iGenes = 0;
	m_numThreads = 0;
	m_mapLoadTime.clear();
	m_output_dir = "";
	DEBUG = false;
	m_bSharedDB = false;
}

bool CSeekCentral::CalculateRestart(){
	set<string> ss;
	m_mapLoadTime.clear();
	ushort i, prev;
	prev = 0;
	for(i=0; i<m_vecstrAllQuery.size(); i++){
		if(m_vecstrAllQuery[i].size()>m_maxNumDB){
			fprintf(stderr, "Cannot fit query %d in buffer\n", i);
			return false;
		}
	}

	vector<vector<string> > *vv = NULL;
	m_mapLoadTime[prev] = vector< vector<string> >();
	for(i=0; i<m_vecstrAllQuery.size(); i++){
		ss.insert(m_vecstrAllQuery[i].begin(), m_vecstrAllQuery[i].end());
		vv = &(m_mapLoadTime[prev]);
		vv->push_back(m_vecstrAllQuery[i]);
		if(ss.size()>m_maxNumDB){
			vv->pop_back();
			ss.clear();
			ss.insert(m_vecstrAllQuery[i].begin(), m_vecstrAllQuery[i].end());
			prev = i;
			m_mapLoadTime[prev] = vector< vector<string> >();
			vv = &(m_mapLoadTime[prev]);
			vv->push_back(m_vecstrAllQuery[i]);
		}
	}
	ss.clear();
	//check
	vector<char> vc;
	ushort tot = 0;
	CSeekTools::InitVector(vc, m_vecstrAllQuery.size(), (char) 0);
	map<ushort, vector< vector<string> > >::const_iterator ci =
		m_mapLoadTime.begin();
	for(; ci!=m_mapLoadTime.end(); ci++) tot+=(ci->second).size();
	vc.clear();
	if(tot!=m_vecstrAllQuery.size()){
		fprintf(stderr, "ERROR, size do not match\n");
		return false;
	}
	return true;
}

//for SeekServer
//assume DB has been read (with gvar, sinfo information)
//assume datasets and genes have been read
bool CSeekCentral::Initialize(string &output_dir, string &query, string &search_dset,
	CSeekCentral *src, float &query_min_required, bool &bCorrelation,
	bool &bSubtractGeneAvg, bool &bSubtractPlatformAvg, bool &bDividePlatformStdev){

	//fprintf(stderr, "B0 %lu\n", CMeta::GetMemoryUsage());
	m_output_dir = output_dir; //LATER, TO BE DELETED
	m_maxNumDB = src->m_maxNumDB;
	m_bSharedDB = true;
	m_numThreads = src->m_numThreads;
	m_fScoreCutOff = src->m_fScoreCutOff;
	m_fPercentQueryAfterScoreCutOff = query_min_required;
	m_bSquareZ = src->m_bSquareZ;
	m_bOutputText = src->m_bOutputText;
	m_bSubtractGeneAvg = bSubtractGeneAvg;
	m_bSubtractPlatformAvg = bSubtractPlatformAvg;
	m_bDividePlatformStdev = bDividePlatformStdev;
	m_bLogit = src->m_bLogit;
	m_bCorrelation = bCorrelation;
	m_vecstrGenes.resize(src->m_vecstrGenes.size());

	m_bRandom = false;
	m_iNumRandom = 1;
	m_randRandom = NULL;	

	copy(src->m_vecstrGenes.begin(), src->m_vecstrGenes.end(), m_vecstrGenes.begin());

	m_vecstrDatasets.resize(src->m_vecstrDatasets.size());
	copy(src->m_vecstrDatasets.begin(), src->m_vecstrDatasets.end(), m_vecstrDatasets.begin());

	m_mapstrintDataset.insert(src->m_mapstrintDataset.begin(), 
		src->m_mapstrintDataset.end());

	m_mapstrintGene.insert(src->m_mapstrintGene.begin(), src->m_mapstrintGene.end());

	m_mapstrstrDatasetPlatform.insert(src->m_mapstrstrDatasetPlatform.begin(),
		src->m_mapstrstrDatasetPlatform.end());
	m_mapstriPlatform.insert(src->m_mapstriPlatform.begin(), src->m_mapstriPlatform.end());

	m_vecstrPlatform.resize(src->m_vecstrPlatform.size());
	copy(src->m_vecstrPlatform.begin(), src->m_vecstrPlatform.end(), m_vecstrPlatform.begin());

	m_vecstrDP.resize(src->m_vecstrDP.size());
	copy(src->m_vecstrDP.begin(), src->m_vecstrDP.end(), m_vecstrDP.begin());

	m_quant = src->m_quant;
	ushort i, j;
	omp_set_num_threads(m_numThreads);

	m_iDatasets = m_vecstrDatasets.size();
	m_iGenes = m_vecstrGenes.size();

	//fprintf(stderr, "%d %d\n", m_iDatasets, m_iGenes);
	//fprintf(stderr, "B1 %lu\n", CMeta::GetMemoryUsage());

	//read search datasets
	vector<string> sd;
	CMeta::Tokenize(search_dset.c_str(), sd, "|", false);
	m_vecstrSearchDatasets.resize(sd.size());
	for(i=0; i<sd.size(); i++){
		CMeta::Tokenize(sd[i].c_str(), m_vecstrSearchDatasets[i], " ", false);
		//fprintf(stderr, "%s\n", sd[i].c_str());
	}
	//read queries
	vector<string> sq;
	CMeta::Tokenize(query.c_str(), sq, "|", false);
	m_vecstrAllQuery.resize(sq.size());
	for(i=0; i<sq.size(); i++){
		CMeta::Tokenize(sq[i].c_str(), m_vecstrAllQuery[i], " ", false);
		//fprintf(stderr, "%s\n", sq[i].c_str());
	}
	//fprintf(stderr, "%s\n", output_dir.c_str());
	m_searchdsetMap.resize(m_vecstrAllQuery.size());
	for(i=0; i<m_vecstrAllQuery.size(); i++){
		m_searchdsetMap[i] = new CSeekIntIntMap(m_vecstrDatasets.size());
		for(j=0; j<m_vecstrSearchDatasets[i].size(); j++)
			m_searchdsetMap[i]->Add(
				m_mapstrintDataset[m_vecstrSearchDatasets[i][j]]);
	}

	//fprintf(stderr, "B2 %lu\n", CMeta::GetMemoryUsage());

	m_DB = src->m_DB; //shared DB

	CSeekTools::LoadDatabase(*m_DB, m_vc, src->m_vc, m_vp, src->m_vp,
		m_vecstrDatasets, m_mapstrstrDatasetPlatform, m_mapstriPlatform);

	//fprintf(stderr, "B3 %lu\n", CMeta::GetMemoryUsage());

	if(!CalculateRestart()) return false;
	return true;
}

//network mode, meant to be run after Initialize()
bool CSeekCentral::EnableNetwork(
	//network parameters
	const int &iClient, const bool &bNetworkSendData){
	m_bEnableNetwork = true;
	m_bNetworkSendStatus = true;
	m_bNetworkSendData = bNetworkSendData;
	m_iClient = iClient; //assume client connection is already open
	return true;
}

//optional step
//Checks how many datasets contain the query
//requires the queries and searchdatasets to be loaded!
bool CSeekCentral::CheckDatasets(const bool &replace){
	ushort dd, j;
	ushort l;
	stringstream ss; //search dataset (new!)
	stringstream sq; //query availability
	stringstream aq; //query (new!)

	for(l=0; l<m_searchdsetMap.size(); l++){
		ushort iUserDatasets = m_searchdsetMap[l]->GetNumSet();
		const vector<ushort> &allRDatasets = m_searchdsetMap[l]->GetAllReverse();	
		vector<int> count;
		CSeekTools::InitVector(count, m_vecstrAllQuery[l].size(), (int) 0);
		bool isFirst = true;

		for(dd=0; dd<iUserDatasets; dd++){
			ushort i = allRDatasets[dd];
			CSeekIntIntMap *si = m_vc[i]->GetGeneMap();
			ushort present = 0;
			for(j=0, present=0; j<m_vecstrAllQuery[l].size(); j++){
				if(m_mapstrintGene.find(m_vecstrAllQuery[l][j])==
					m_mapstrintGene.end()) continue;
				if(CSeekTools::IsNaN(si->GetForward(
					m_mapstrintGene[m_vecstrAllQuery[l][j]]))) continue;
				count[j]++;
				present++;
			}
			//datasets that contains all query genes (very stringent)
			//if(present==m_vecstrAllQuery[l].size()){

			int minRequired = (int) (m_fPercentQueryAfterScoreCutOff * 
				m_vecstrAllQuery[l].size());
			//datasets containing some query genes (relaxed) [ 0 ] 
			if(present>0 && present>=minRequired){
				if(isFirst){
					isFirst = false;
					ss << m_vecstrDatasets[i];
				}else{
					ss << " " << m_vecstrDatasets[i];
				}
			}
		}

		if(isFirst){
			string err = "Error: no dataset contains any of the query genes";
			fprintf(stderr, "%s\n", err.c_str());
			if(m_bEnableNetwork)
				CSeekNetwork::Send(m_iClient, err);
			return false;
		}

		if(l!=m_searchdsetMap.size()-1){
			ss << "|";
		}

		isFirst = true;		
		for(j=0; j<m_vecstrAllQuery[l].size(); j++){
			sq << m_vecstrAllQuery[l][j] << ":" << count[j];
			if(j!=m_vecstrAllQuery[l].size()-1){
				sq << ";";
			}
			if(count[j]==0) continue;
			if(isFirst){
				isFirst = false;
				aq << m_vecstrAllQuery[l][j];
			}else{
				aq << " " << m_vecstrAllQuery[l][j];
			}
		}

		if(isFirst){
			string err = "Error: no dataset contains any of the query genes";
			fprintf(stderr, "%s\n", err.c_str());
			if(m_bEnableNetwork)
				CSeekNetwork::Send(m_iClient, err);
			return false;
		}

		if(l!=m_searchdsetMap.size()-1){
			aq << "|";
			sq << "|";
		}

	}

	string refinedQuery = aq.str();
	string refinedSearchDataset = ss.str();
	string refinedGeneCount = sq.str();
	if(m_bEnableNetwork){
		CSeekNetwork::Send(m_iClient, refinedSearchDataset);
		CSeekNetwork::Send(m_iClient, refinedGeneCount);
	}

	if(replace){
		vector<string> qq;
		ushort i;
		CMeta::Tokenize(refinedQuery.c_str(), qq, "|", true);
		m_vecstrAllQuery.resize(qq.size());
		for(i=0; i<qq.size(); i++){
			m_vecstrAllQuery[i].clear();
			CMeta::Tokenize(qq[i].c_str(), m_vecstrAllQuery[i], " ", true);
		}

		//Change the search datasets
		vector<string> sd;
		CMeta::Tokenize(refinedSearchDataset.c_str(), sd, "|", false);
		m_vecstrSearchDatasets.resize(sd.size());
		for(i=0; i<sd.size(); i++){
			m_vecstrSearchDatasets[i].clear();
			CMeta::Tokenize(sd[i].c_str(), m_vecstrSearchDatasets[i], " ", false);
		}

		for(i=0; i<m_searchdsetMap.size(); i++){
			delete m_searchdsetMap[i];
		}
		m_searchdsetMap.clear();

		m_searchdsetMap.resize(m_vecstrAllQuery.size());
		for(i=0; i<m_vecstrAllQuery.size(); i++){
			m_searchdsetMap[i] = new CSeekIntIntMap(m_vecstrDatasets.size());
			for(j=0; j<m_vecstrSearchDatasets[i].size(); j++)
				m_searchdsetMap[i]->Add(
					m_mapstrintDataset[m_vecstrSearchDatasets[i][j]]);
		}

	}

	return true;	
}

//load everything except query, search datasets, output directory
bool CSeekCentral::Initialize(const char *gene, const char *quant,
	const char *dset, const char *platform, const char *db,
	const char *prep, const char *gvar, const char *sinfo,
	const bool &useNibble, const ushort &num_db,
	const ushort &buffer, const bool &to_output_text,
	const bool &bCorrelation, const bool &bSubtractAvg,
	const bool &bSubtractPlatformAvg, const bool &bDividePlatformStdev,
	const bool &bLogit, const float &fCutOff, const float &fPercentRequired, 
	const bool &bSquareZ, 
	//three new ones
	const bool &bRandom, const int &iNumRandom, gsl_rng *rand){

	m_maxNumDB = buffer;
	m_numThreads = 8; //changed from 8
	m_fScoreCutOff = fCutOff;
	m_fPercentQueryAfterScoreCutOff = fPercentRequired;
	m_bSquareZ = bSquareZ;

	//random retrieval==========================
	m_bRandom = bRandom;
	m_iNumRandom = iNumRandom;
	m_randRandom = rand; //random-case only
	//===

	if(!m_bRandom){
		m_randRandom = NULL;
	}

	ushort i, j;

	omp_set_num_threads(m_numThreads);

	m_bOutputText = to_output_text;
	m_bSubtractGeneAvg = bSubtractAvg;
	m_bSubtractPlatformAvg = bSubtractPlatformAvg;
	m_bDividePlatformStdev = bDividePlatformStdev;
	m_bLogit = bLogit;
	m_bCorrelation = bCorrelation;

	string strGvarDirectory = gvar;
	string strSinfoDirectory = sinfo;
	if(m_bCorrelation && sinfo=="NA"){
		fprintf(stderr, "Error: not specifying sinfo!\n");
		return false;
	}

	if(m_bCorrelation && (m_bSubtractGeneAvg || m_bSubtractPlatformAvg ||
		m_bDividePlatformStdev || m_bLogit)){
		fprintf(stderr, 
			"Warning: setting subtract_avg, subtract_platform to false\n");
		m_bSubtractGeneAvg = false;
		m_bSubtractPlatformAvg = false;
		m_bDividePlatformStdev = false;
		m_bLogit = false;
	}

	//read genes
	vector<string> vecstrGeneID;
	if(!CSeekTools::ReadListTwoColumns(gene, vecstrGeneID, m_vecstrGenes))
		return false;

	for(i=0; i<m_vecstrGenes.size(); i++)
		m_mapstrintGene[m_vecstrGenes[i]] = i;

	CSeekTools::ReadQuantFile(quant, m_quant);
	m_DB = new CDatabase(useNibble);

	//read datasets
	if(!CSeekTools::ReadListTwoColumns(dset, m_vecstrDatasets, m_vecstrDP))
		return false;

	for(i=0; i<m_vecstrDatasets.size(); i++){
		m_mapstrstrDatasetPlatform[m_vecstrDatasets[i]] = m_vecstrDP[i];
		m_mapstrintDataset[m_vecstrDatasets[i]] = i;
	}

	vector<string> vecstrPlatforms;
	CSeekTools::ReadPlatforms(platform, m_vp, vecstrPlatforms,
		m_mapstriPlatform);

	m_iDatasets = m_vecstrDatasets.size();
	m_iGenes = m_vecstrGenes.size();

	m_DB->Open(db, m_vecstrGenes, m_iDatasets, num_db);
	CSeekTools::LoadDatabase(*m_DB, prep, gvar, sinfo, m_vecstrDatasets,
		m_mapstrstrDatasetPlatform, m_mapstriPlatform, m_vp, m_vc);

	return true;
}


bool CSeekCentral::Initialize(const char *gene, const char *quant,
	const char *dset, const char *search_dset,
	const char *query, const char *platform, const char *db,
	const char *prep, const char *gvar, const char *sinfo,
	const bool &useNibble, const ushort &num_db,
	const ushort &buffer, const char *output_dir, const bool &to_output_text,
	const bool &bCorrelation, const bool &bSubtractAvg,
	const bool &bSubtractPlatformAvg, const bool &bDividePlatformStdev,
	const bool &bLogit, const float &fCutOff, const float &fPercentRequired, 
	const bool &bSquareZ,
	//three new ones
	const bool &bRandom, const int &iNumRandom, gsl_rng *rand){

	if(!CSeekCentral::Initialize(gene, quant, dset, platform, 
		db, prep, gvar, sinfo, useNibble, num_db, buffer, to_output_text, 
		bCorrelation, bSubtractAvg, bSubtractPlatformAvg, bDividePlatformStdev, 
		bLogit, fCutOff, fPercentRequired, bSquareZ, bRandom, iNumRandom, rand)){
		return false;
	}

	ushort i, j;
	omp_set_num_threads(m_numThreads);
	m_output_dir = output_dir;

	//read search datasets
	if(!CSeekTools::ReadMultipleQueries(search_dset, m_vecstrSearchDatasets))
		return false;
	//read queries
	if(!CSeekTools::ReadMultipleQueries(query, m_vecstrAllQuery))
		return false;

	m_searchdsetMap.resize(m_vecstrAllQuery.size());
	for(i=0; i<m_vecstrAllQuery.size(); i++){
		m_searchdsetMap[i] = new CSeekIntIntMap(m_vecstrDatasets.size());
		for(j=0; j<m_vecstrSearchDatasets[i].size(); j++)
			m_searchdsetMap[i]->Add(
				m_mapstrintDataset[m_vecstrSearchDatasets[i][j]]);
	}

	if(!CalculateRestart()) return false;


	return true;
}

bool CSeekCentral::PrepareQuery(const vector<string> &vecstrQuery,
	CSeekQuery &query){
	vector<ushort> queryGenes;
	ushort j;
	for(j=0; j<vecstrQuery.size(); j++){
		size_t m = m_DB->GetGene(vecstrQuery[j]);
		if(m==-1) continue;
		queryGenes.push_back(m);
	}
	queryGenes.resize(queryGenes.size());
	query.InitializeQuery(queryGenes, m_iGenes);

	return true;
}

bool CSeekCentral::PrepareOneQuery(CSeekQuery &query,
	CSeekIntIntMap &dMap, vector<float> &weight){

	assert(m_master_rank_threads==NULL && m_counts_threads==NULL &&
		m_sum_weight_threads==NULL);
	assert(m_rank_normal_threads==NULL && m_rank_threads==NULL);
	assert(m_rData==NULL);

	ushort j;
	const vector<ushort> &queryGenes = query.GetQuery();
	const vector<ushort> &allRDatasets = dMap.GetAllReverse();
	ushort iSearchDatasets = dMap.GetNumSet();
	ushort iQuery = queryGenes.size();

	for(j=0; j<iSearchDatasets; j++)
		m_vc[allRDatasets[j]]->InitializeQuery(queryGenes);

	m_rData = new ushort**[m_numThreads];
	for(j=0; j<m_numThreads; j++)
		m_rData[j] = CSeekTools::Init2DArray(m_iGenes, iQuery, (ushort)0);
	
	m_master_rank_threads =
		CSeekTools::Init2DArray(m_numThreads, m_iGenes, (float)0);
	m_sum_weight_threads =
		CSeekTools::Init2DArray(m_numThreads, m_iGenes, (float)0);
	m_counts_threads =
		CSeekTools::Init2DArray(m_numThreads, m_iGenes, (ushort)0);
	
	m_rank_normal_threads = new vector<ushort>[m_numThreads];
	m_rank_threads = new vector<ushort>[m_numThreads];

	for(j=0; j<m_numThreads; j++){
		m_rank_normal_threads[j].resize(m_iGenes);
		m_rank_threads[j].resize(m_iGenes);
		//CSeekTools::InitVector(m_rank_normal_threads[j], m_iGenes, (ushort) 255);
		//CSeekTools::InitVector(m_rank_threads[j], m_iGenes, (ushort) 255);
	}
	
	CSeekTools::InitVector(m_master_rank, m_iGenes, (float) 0);
	CSeekTools::InitVector(m_sum_weight, m_iGenes, (float) 0);
	CSeekTools::InitVector(m_counts, m_iGenes, (ushort) 0);
	CSeekTools::InitVector(weight, m_iDatasets, (float)0);
	
	return true;
}

bool CSeekCentral::AggregateThreads(){
	assert(m_master_rank_threads!=NULL && m_counts_threads!=NULL &&
		m_sum_weight_threads!=NULL);
	assert(m_rank_normal_threads!=NULL && m_rank_threads!=NULL);

	//Aggregate into three vectors: m_master_rank, m_counts, m_sum_weight
	ushort j, k;
	for(j=0; j<m_numThreads; j++){
		for(k=0; k<m_iGenes; k++){
			m_master_rank[k] += m_master_rank_threads[j][k];
			m_counts[k] += m_counts_threads[j][k];
			m_sum_weight[k]+=m_sum_weight_threads[j][k];
		}
	}

	CSeekTools::Free2DArray(m_master_rank_threads);
	CSeekTools::Free2DArray(m_counts_threads);
	CSeekTools::Free2DArray(m_sum_weight_threads);
	m_master_rank_threads=NULL;
	m_counts_threads=NULL;
	m_sum_weight_threads = NULL;

	for(j=0; j<m_numThreads; j++){
		m_rank_normal_threads[j].clear();
		m_rank_threads[j].clear();
	}

	delete[] m_rank_normal_threads;
	delete[] m_rank_threads;
	m_rank_normal_threads = NULL;
	m_rank_threads = NULL;

	return true;
}

bool CSeekCentral::FilterResults(const ushort &iSearchDatasets){
	ushort j, k;
	bool DEBUG = false;
	if(DEBUG) fprintf(stderr, "Aggregating genes\n");

	for(j=0; j<m_iGenes; j++){
		if(m_counts[j]<(int)(0.5*iSearchDatasets))
			m_master_rank[j] = -320;
		else if(m_sum_weight[j]==0)
			m_master_rank[j] = -320;
		else{
			m_master_rank[j] =
				(m_master_rank[j] / m_sum_weight[j] - 320) / 100.0;
			if(m_bCorrelation){
				m_master_rank[j] = m_master_rank[j] / 3.0;
			}
			//m_master_rank[j] = m_master_rank[j];
			//m_master_rank[j] = (m_master_rank[j] / iSearchDatasets - 320) / 100.0;
		}
		if(DEBUG) fprintf(stderr, "Gene %d %.5f\n", j, m_master_rank[j]);
	}
	return true;
}

bool CSeekCentral::Sort(vector<AResultFloat> &final){
	if(DEBUG) fprintf(stderr, "Sorting genes\n");
	final.resize(m_iGenes);
	ushort j;
	for(j=0; j<m_iGenes; j++){
		//fprintf(stderr, "%d %s\n", j, DB.GetGene((size_t) j).c_str());
		final[j].i = j;
		final[j].f = m_master_rank[j];
	}
	if(DEBUG) fprintf(stderr, "Begin Sorting genes\n");
	sort(final.begin(), final.end());
	return true;
}

bool CSeekCentral::Display(CSeekQuery &query, vector<AResultFloat> &final){
	if(DEBUG) fprintf(stderr, "Results:\n");
	ushort jj, ii;
	const vector<char> &cQuery = query.GetQueryPresence();
	for(ii=0, jj=0; jj<500; ii++){
		if(cQuery[final[ii].i]==1) continue;
		fprintf(stderr, "%s %.5f\n",
			m_DB->GetGene((size_t)final[ii].i).c_str(), final[ii].f);
		jj++;
	}
	return true;
}

bool CSeekCentral::Write(const ushort &i){
	//assume m_bRandom = false
	char acBuffer[1024];
	sprintf(acBuffer, "%s/%d.query", m_output_dir.c_str(), i);
	CSeekTools::WriteArrayText(acBuffer, m_vecstrAllQuery[i]);

	sprintf(acBuffer, "%s/%d.dweight", m_output_dir.c_str(), i);
	CSeekTools::WriteArray(acBuffer, m_weight[i]);
	
	sprintf(acBuffer, "%s/%d.gscore", m_output_dir.c_str(), i);
	CSeekTools::WriteArray(acBuffer, m_master_rank);

	//send data to client
	if(m_bEnableNetwork){
		if(CSeekNetwork::Send(m_iClient, m_weight[i])==-1){
			fprintf(stderr, "Error sending message to client\n");
			return false;
		}
		if(CSeekNetwork::Send(m_iClient, m_master_rank)==-1){
			fprintf(stderr, "Error sending message to client\n");
			return false;
		}
	}

	if(!m_bRandom && m_bOutputText){
		const vector<ushort> &allRDatasets =
			m_searchdsetMap[i]->GetAllReverse();
		ushort iSearchDatasets = m_searchdsetMap[i]->GetNumSet();
		vector<vector<string> > vecOutput;
		vecOutput.resize(2);
		vecOutput[0] = vector<string>();
		vecOutput[1] = vector<string>();
		sprintf(acBuffer, "%s/%d.results.txt", m_output_dir.c_str(), i);
		vector<AResultFloat> w;
		w.resize(m_iDatasets);
		ushort j;
		for(j=0; j<m_iDatasets; j++){
			w[j].i = j;
			w[j].f = m_weight[i][j];
		}
		sort(w.begin(), w.end());
		for(j=0; j<200 && j<iSearchDatasets; j++){
			if(w[j].f==0) break;
			vecOutput[0].push_back(m_vecstrDatasets[w[j].i]);
		}
		vector<AResultFloat> wd;
		Sort(wd);
		for(j=0; j<2000; j++){
			if(wd[j].f==-320) break;
			vecOutput[1].push_back(m_vecstrGenes[wd[j].i]);
		}
		CSeekTools::Write2DArrayText(acBuffer, vecOutput);
	}
	return true;
}


bool CSeekCentral::Common(CSeekCentral::SearchMode &sm,
	gsl_rng *rnd, const CSeekQuery::PartitionMode *PART_M,
	const ushort *FOLD, const float *RATE,
	const vector< vector<float> > *providedWeight,
	const vector< vector<string> > *newGoldStd){

	ushort i, j, d, dd;
	int k; //keeps track of genes (for random case)
	ushort l; //keeps track of random repetition (for random case)
	char acBuffer[1024];
	
	m_Query.resize(m_vecstrAllQuery.size());
	m_weight.resize(m_vecstrAllQuery.size());
	m_final.resize(m_vecstrAllQuery.size());
		
	//random-ranking case =========================
	vector<vector<float> > vecRandWeight, vecRandScore;
	vecRandWeight.resize(m_iNumRandom);
	vecRandScore.resize(m_iNumRandom);
	for(l=0; l<m_iNumRandom; l++){
		CSeekTools::InitVector(vecRandWeight[l], m_iDatasets, (float) 0);
		CSeekTools::InitVector(vecRandScore[l], m_iGenes, (float) 0);
	}

	//NEED TO MAKE THIS A PARAMETER
	bool simulateWeight = true;

	//output weight component (Mar 19)
	bool weightComponent = true;

	l = 0;
	//oct 20, 2012: whether to redo current query with equal weighting
	int redoWithEqual = 0; //tri-mode: 0, 1, 2
	CSeekCentral::SearchMode current_sm;
	CSeekQuery equalWeightGold;

	//backup of scores (Feb 3)
	vector<float> backupScore;
	CSeekTools::InitVector(backupScore, m_iGenes, (float) 0);

	//fprintf(stderr, "0 %lu\n", CMeta::GetMemoryUsage());
	current_sm = sm;

	for(i=0; i<m_vecstrAllQuery.size(); i++){


		//simulated weight case ======================
		/*if(simulateWeight && redoWithEqual>=1) //1 or 2 
			current_sm = EQUAL;
		else //0
			current_sm = sm;*/
		//============================================

		if(m_mapLoadTime.find(i)!=m_mapLoadTime.end()){
			if(!m_bRandom || l==0){ //l==0: first random repetition
				CSeekTools::ReadDatabaselets(*m_DB, m_mapLoadTime[i], m_vc, 
				m_iClient, m_bEnableNetwork);
			}
		}

		const vector<ushort> &allRDatasets =
			m_searchdsetMap[i]->GetAllReverse();
		ushort iSearchDatasets = m_searchdsetMap[i]->GetNumSet();

		//fprintf(stderr, "1 %lu\n", CMeta::GetMemoryUsage());

		CSeekQuery &query = m_Query[i];
		vector<float> &weight = m_weight[i];
		vector<AResultFloat> &final = m_final[i];

		PrepareQuery(m_vecstrAllQuery[i], query);
		PrepareOneQuery(query, *(m_searchdsetMap[i]), weight);
		ushort iQuery = query.GetQuery().size();

		//fprintf(stderr, "1b %lu\n", CMeta::GetMemoryUsage());
		
		//for CV_CUSTOM
		CSeekQuery customGoldStd;
		if(current_sm==CV_CUSTOM)
			PrepareQuery((*newGoldStd)[i], customGoldStd);

		if(current_sm==CV || current_sm==CV_CUSTOM)
			query.CreateCVPartitions(rnd, *PART_M, *FOLD);

		if(current_sm==ORDER_STATISTICS)
			m_rank_d = CSeekTools::Init2DArray(iSearchDatasets, m_iGenes,
				(ushort) 0);

		//For outputing component weights!
		vector<float> wc;
		if(weightComponent){
			if(current_sm==CV || current_sm==CV_CUSTOM){
				wc.resize((int)query.GetNumFold()*(int)m_iDatasets);
			}else{
				wc.resize((int)query.GetQuery().size()*(int)m_iDatasets);
			}
			fill(wc.begin(), wc.end(), (float)0);
		}


		//fprintf(stderr, "2 %lu\n", CMeta::GetMemoryUsage());

		#pragma omp parallel for \
		shared(customGoldStd) \
		private(dd, d, j) \
		firstprivate(iSearchDatasets, iQuery) \
		schedule(dynamic)
		for(dd=0; dd<iSearchDatasets; dd++){
			d = allRDatasets[dd];
			ushort tid = omp_get_thread_num();
			if(DEBUG) fprintf(stderr, "Dataset %d, %s\n",
				d, m_vecstrDatasets[d].c_str());

			CSeekIntIntMap *mapG = m_vc[d]->GetGeneMap();
			CSeekIntIntMap *mapQ = m_vc[d]->GetQueryMap();

			//if(mapG->GetNumSet()<10000){
			//	continue;
			//}

			if(mapQ==NULL ||mapQ->GetNumSet()==0){
				if(DEBUG) fprintf(stderr, "This dataset is skipped\n");
				continue;
			}

			vector<ushort> this_q;
			for(j=0; j<mapQ->GetNumSet(); j++)
				this_q.push_back(mapQ->GetReverse(j));

			if(DEBUG) fprintf(stderr, "Initializing %d\n",
				(int) this_q.size());
			m_vc[d]->InitializeDataMatrix(m_rData[tid], m_quant, m_iGenes,
				iQuery, m_bSubtractGeneAvg, m_bSubtractPlatformAvg, m_bLogit,
				m_bCorrelation, m_fScoreCutOff, m_bRandom, m_randRandom);
			//m_bSubtractPlatformStdev is not used, it's assumed

			float w = -1;
			float report_w = -1; //for showing weight of dataset

			if(current_sm==CV || current_sm==CV_CUSTOM){
				if(DEBUG) fprintf(stderr, "Weighting dataset\n");
				if(current_sm==CV)
					CSeekWeighter::CVWeighting(query, *m_vc[d], *RATE,
						m_fPercentQueryAfterScoreCutOff, m_bSquareZ,
						&m_rank_threads[tid]);
				else
					CSeekWeighter::CVWeighting(query, *m_vc[d], *RATE,
						m_fPercentQueryAfterScoreCutOff, m_bSquareZ,
						&m_rank_threads[tid], &customGoldStd);

				if( (w = m_vc[d]->GetDatasetSumWeight())==-1){
					if(DEBUG) fprintf(stderr, "Bad weight\n");
					continue;
				}

				if(weightComponent && current_sm==CV){
					ushort numFold = query.GetNumFold();
					float ww;
					for(j=0; j<numFold; j++){
						//if((ww=m_vc[d]->GetCVWeight(j))==-1) continue;
						wc[(int)d*(int)numFold+(int)j] = m_vc[d]->GetCVWeight(j);
					}
				}
			}
			else if(current_sm==EQUAL && redoWithEqual==0){
				w = 1.0;
			}
			else if(current_sm==USE_WEIGHT){
				w = (*providedWeight)[i][d];
			}
			//simulated weight case ======================
			else if(current_sm==EQUAL && redoWithEqual==1){
				w = 1.0;
				if(DEBUG) fprintf(stderr, "Before doing one gene weighting\n");
				//calculate reported weight here!
				CSeekWeighter::OneGeneWeighting(query, *m_vc[d], 0.95,
					m_fPercentQueryAfterScoreCutOff, m_bSquareZ,
					&m_rank_threads[tid], &equalWeightGold);
				report_w = m_vc[d]->GetDatasetSumWeight();
			}
			//============================================

			if(DEBUG) fprintf(stderr, "Doing linear combination\n");

			const ushort MIN_REQUIRED = max((ushort) 1, (ushort) (
				m_fPercentQueryAfterScoreCutOff * this_q.size()));
			CSeekWeighter::LinearCombine(m_rank_normal_threads[tid], this_q,
				*m_vc[d], MIN_REQUIRED, m_bSquareZ);

			if(DEBUG) fprintf(stderr,
				"Adding contribution of dataset %d to master ranking: %.5f\n", d, w);

			ushort iGeneSet = mapG->GetNumSet();
			const vector<ushort> &allRGenes = mapG->GetAllReverse();
			vector<ushort>::const_iterator iterR = allRGenes.begin();
			vector<ushort>::const_iterator endR = allRGenes.begin() + iGeneSet;
			vector<ushort> &Rank_Normal = m_rank_normal_threads[tid];
			float* Master_Rank = &m_master_rank_threads[tid][0];
			float* Sum_Weight = &m_sum_weight_threads[tid][0];
			ushort* Counts = &m_counts_threads[tid][0];

			if(current_sm==ORDER_STATISTICS)
				for(; iterR!=endR; iterR++){
					//if(Rank_Normal[*iterR]==0) continue;
					m_rank_d[dd][*iterR] = Rank_Normal[*iterR];
					Counts[*iterR]++;
				}
			else
				for(; iterR!=endR; iterR++){
					//if(Rank_Normal[*iterR]==0) continue;
					Master_Rank[*iterR] += (float) Rank_Normal[*iterR] * w;
					Sum_Weight[*iterR] += w;
					Counts[*iterR]++;
				}

			//simulated weight case ======================
			if(current_sm==EQUAL && redoWithEqual==1){
				weight[d] = report_w;
			}
			//============================================
			else if((current_sm==EQUAL || current_sm==ORDER_STATISTICS) 
			&& redoWithEqual==0){
				weight[d] = 0;
			}
			else{
				weight[d] = w;
			}

		}
		//omp finishes
		//fprintf(stderr, "3 %lu\n", CMeta::GetMemoryUsage());
		for(j=0; j<iSearchDatasets; j++)
			m_vc[allRDatasets[j]]->DeleteQuery();

		assert(m_rData!=NULL);
		for(j=0; j<m_numThreads; j++)
			CSeekTools::Free2DArray(m_rData[j]);
		delete[] m_rData;
		m_rData = NULL;

		AggregateThreads();

		if(current_sm!=ORDER_STATISTICS){
			FilterResults(iSearchDatasets);
		}else{
			CSeekWeighter::OrderStatisticsRankAggregation(iSearchDatasets,
				m_iGenes, m_rank_d, m_counts, m_master_rank, m_numThreads);
			CSeekTools::Free2DArray(m_rank_d);
			m_rank_d = NULL;
		}

		//Display(query, final);
		//fprintf(stderr, "4 %lu\n", CMeta::GetMemoryUsage());
		SetQueryScoreNull(query);
		Sort(final);

		if(m_bRandom){
			/*ushort z, cz;
			for(cz=0, z=0; z<m_iDatasets; z++)
				if(m_weight[i][z]!=0) 
					cz++;
			fprintf(stderr, "Number of weighted dataset: %d\n", cz);
			*/
		}
		else if(simulateWeight){
			if((current_sm==EQUAL || current_sm==ORDER_STATISTICS) && !CheckWeight(i)){
				fprintf(stderr, "Calculate dataset ordering\n"); system("date +%s%N 1>&2");
				if(m_bEnableNetwork && CSeekNetwork::Send(m_iClient, 
					"Calculate dataset ordering")==-1){
					fprintf(stderr, "Error sending message to client\n");
				}
				CopyTopGenes(equalWeightGold, final, 100);
				redoWithEqual = 1;
				if(current_sm==ORDER_STATISTICS){
					current_sm = EQUAL;
					//backup genes to old
					copy(m_master_rank.begin(), m_master_rank.end(), backupScore.begin());	
				}
				i--;
				continue;
			}
			else if(current_sm==CV && !CheckWeight(i)){
				fprintf(stderr, "Redo with equal weighting\n"); system("date +%s%N 1>&2");
				if(m_bEnableNetwork && CSeekNetwork::Send(m_iClient, 
					"Redo with equal weighting")==-1){
					fprintf(stderr, "Error sending message to client\n");
				}
				current_sm = EQUAL;
				i--;
				continue;
			}
			else if(current_sm==EQUAL && CheckWeight(i)){
				redoWithEqual = 0;
				current_sm = sm;
				//copy genes from old to new
				if(sm==ORDER_STATISTICS){
					copy(backupScore.begin(), backupScore.end(), m_master_rank.begin());	
				}
			}
		}

		fprintf(stderr, "Done search\n"); system("date +%s%N 1>&2");

		if(m_bEnableNetwork && CSeekNetwork::Send(m_iClient, "Done Search")==-1){
			fprintf(stderr, "Error sending message to client\n");
		}
	
		//random-ranking case =========================
		if(m_bRandom){
			sort(m_master_rank.begin(), m_master_rank.end(), greater<float>());
			sort(weight.begin(), weight.end(), greater<float>());
			copy(m_master_rank.begin(), m_master_rank.end(), vecRandScore[l].begin());
			copy(weight.begin(), weight.end(), vecRandWeight[l].begin());
			l++;

			if(l==m_iNumRandom){ //last repetition
				float confidence = 0.50;
				int n = (int) (confidence * (float) m_iNumRandom); 
				vector<float> new_weight, new_score;
				CSeekTools::InitVector(new_weight, m_iDatasets, (float) 0);
				CSeekTools::InitVector(new_score, m_iGenes, (float) 0);
				for(k=0; k<m_iGenes; k++){
					vector<float> all_s;
					all_s.resize(m_iNumRandom);
					for(l=0; l<m_iNumRandom; l++)
						all_s[l] = vecRandScore[l][k];
					std::nth_element(all_s.begin(), all_s.begin()+n, all_s.end());
					new_score[k] = all_s[n];
				}
				/*for(k=0; k<m_iGenes; k++){
					fprintf(stderr, "%d %.5f\n", k, new_score[k]);
				}
				getchar();*/
				for(k=0; k<m_iDatasets; k++){
					vector<float> all_w;
					all_w.resize(m_iNumRandom);
					for(l=0; l<m_iNumRandom; l++)
						all_w[l] = vecRandWeight[l][k];	
					std::nth_element(all_w.begin(), all_w.begin()+n, all_w.end());
					new_weight[k] = all_w[n];
				}
				sprintf(acBuffer, "%s/%d.rand.dweight", m_output_dir.c_str(), i);
				CSeekTools::WriteArray(acBuffer, new_weight);
				sprintf(acBuffer, "%s/%d.rand.gscore", m_output_dir.c_str(), i);
				CSeekTools::WriteArray(acBuffer, new_score);				
				l = 0;
			}else{
				i--;
			}
		}

		if(!m_bRandom){
			//if m_bRandom, write at the very end when all repetitions are done
			Write(i);
			if(weightComponent){
				sprintf(acBuffer, "%s/%d.dweight_comp", m_output_dir.c_str(), i);
				CSeekTools::WriteArray(acBuffer, wc);
				vector<vector<string> > vecParts;
				vecParts.resize(query.GetNumFold());
				ushort kk;
				string strParts = "";
				for(j=0; j<query.GetNumFold(); j++){
					vecParts[j] = vector<string>();
					const vector<ushort> &vu = query.GetCVQuery(j);
					string s = "";
					for(kk=0; kk<vu.size(); kk++){
						vecParts[j].push_back(m_vecstrGenes[vu[kk]]);
						s+=m_vecstrGenes[vu[kk]];
						if(kk!=vu.size()-1)
							s+=" ";
					}
					strParts+=s;
					if(j!=query.GetNumFold()-1)
						strParts+="|";
				}
				sprintf(acBuffer, "%s/%d.query_cvpart", m_output_dir.c_str(), i);
				CSeekTools::Write2DArrayText(acBuffer, vecParts);
				//send data to client
				if(m_bEnableNetwork){
					if(CSeekNetwork::Send(m_iClient, wc)==-1){
						fprintf(stderr, "Error sending message to client\n");
						return false;
					}
					if(CSeekNetwork::Send(m_iClient, strParts)==-1){
						fprintf(stderr, "Error sending message to client\n");
						return false;
					}
				}
			}
		}
	}

	//fprintf(stderr, "4b %lu\n", CMeta::GetMemoryUsage());

	return true;
}

bool CSeekCentral::CheckWeight(const ushort &i){
	ushort j = 0;
	bool valid = false;
	for(j=0; j<m_iDatasets; j++){
		if(m_weight[i][j]!=0){
			valid = true;
			break;
		}
	}
	return valid;
}

bool CSeekCentral::CopyTopGenes(CSeekQuery &csq, 
	const vector<AResultFloat> &src, const ushort top){
	ushort i, j;
	vector<ushort> topGenes;
	for(i=0; i<top; i++){
		if(src[i].f==-320) continue;
		topGenes.push_back(src[i].i);
	}
	if(topGenes.size()==0){
		fprintf(stderr, "Error in CopyTopGenes!\n");
		return false;
	}
	csq.InitializeQuery(topGenes, m_iGenes);
	return true;
}

bool CSeekCentral::SetQueryScoreNull(const CSeekQuery &csq){
	ushort j;
	const vector<ushort> &query = csq.GetQuery();
	for(j=0; j<query.size(); j++){
		m_master_rank[query[j]] = -320;
	}
	return true;	
}

bool CSeekCentral::EqualWeightSearch(){
	CSeekCentral::SearchMode sm = EQUAL;
	CSeekCentral::Common(sm);
}

bool CSeekCentral::CVSearch(gsl_rng *rnd, const CSeekQuery::PartitionMode &PART_M,
	const ushort &FOLD, const float &RATE){
	CSeekCentral::SearchMode sm = CV;
	CSeekCentral::Common(sm, rnd, &PART_M, &FOLD, &RATE);
}

/*	perform CVSearch, except that the weighting is based not on co-expression
	of query genes, but based on similarity of query genes to some custom gold
	standard gene-set */
bool CSeekCentral::CVCustomSearch(const vector< vector<string> > &newGoldStd,
	gsl_rng *rnd, const CSeekQuery::PartitionMode &PART_M,
	const ushort &FOLD, const float &RATE){
	CSeekCentral::SearchMode sm = CV_CUSTOM;
	CSeekCentral::Common(sm, rnd, &PART_M, &FOLD, &RATE,
		NULL, &newGoldStd);
}

bool CSeekCentral::WeightSearch(const vector< vector<float> > &weights){
	CSeekCentral::SearchMode sm = USE_WEIGHT;
	CSeekCentral::Common(sm, NULL, NULL, NULL, NULL, &weights);
}

bool CSeekCentral::OrderStatistics(){
	CSeekCentral::SearchMode sm = ORDER_STATISTICS;
	CSeekCentral::Common(sm);
}

/* to be implemented */
bool CSeekCentral::SingleGeneMetaCorrelation(){
	CSeekCentral::SearchMode sm = SINGLE_GENE_META;
	return false;
}

bool CSeekCentral::VarianceWeightSearch(){
	vector<vector<float> > weights;
	weights.resize(m_vecstrAllQuery.size());
	ushort i, j, k;
	for(i=0; i<m_vecstrAllQuery.size(); i++){
		weights[i] = vector<float>();
		CSeekTools::InitVector(weights[i], m_iDatasets, (float)0);
		const vector<ushort> &allRDatasets =
			m_searchdsetMap[i]->GetAllReverse();
		ushort iSearchDatasets = m_searchdsetMap[i]->GetNumSet();

		CSeekQuery q;
		PrepareQuery(m_vecstrAllQuery[i], q);
		const vector<ushort> &qGenes = q.GetQuery();

		for(j=0; j<iSearchDatasets; j++){
			ushort d = allRDatasets[j];
			for(k=0; k<qGenes.size(); k++){
				float gv = m_vc[d]->GetGeneVariance(qGenes[k]);
				if(CMeta::IsNaN(gv)) continue;
				weights[i][d] += gv;
			}
		}
	}
	return WeightSearch(weights);
}

bool CSeekCentral::Destruct(){
	ushort j;
	//for(j=0; j<m_iDatasets; j++) m_vc[j]->DeleteQueryBlock();
	for(j=0; j<m_iDatasets; j++){
		if(m_vc[j]!=NULL) continue;
		delete m_vc[j];
		m_vc[j] = NULL;
	}
	return true;
}

const vector< vector<AResultFloat> >& CSeekCentral::GetAllResult() const{
	return m_final;
}

const vector<CSeekQuery>& CSeekCentral::GetAllQuery() const{
	return m_Query;
}

ushort CSeekCentral::GetGene(const string &strGene) const{
	return (ushort) m_DB->GetGene(strGene);
}
string CSeekCentral::GetGene(const ushort &geneID) const{
	return m_DB->GetGene((size_t) geneID);
}

const vector<vector<float> >& CSeekCentral::GetAllWeight() const{
	return m_weight;
}
}

