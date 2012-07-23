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
#include "stdafx.h"
#include "seekcentral.h"
#include "seekevaluate.h"
#include "database.h"
#include "datapair.h"
#include "seekweight.h"


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
	m_dsetMap = NULL;
	m_DB = NULL;
	m_rData = NULL;
	m_maxNumDB = 50;

	m_master_rank_threads = NULL;
	m_sum_weight_threads = NULL;
	m_counts_threads = NULL;
	m_rank_normal_threads = NULL;
	m_rank_threads = NULL;

	m_master_rank.clear();
	m_sum_weight.clear();
	m_weight.clear();
	m_counts.clear();
	m_mapLoadTime.clear();

	m_final.clear();

	m_iDatasets = 0;
	m_iGenes = 0;
	m_numThreads = 0;
	m_bSubtractGeneAvg = false;
	m_bSubtractPlatformAvg = false;
	m_bDividePlatformStdev = false;
	m_bLogit = false;
	DEBUG = false;
}

CSeekCentral::~CSeekCentral(){
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
	m_mapLoadTime.clear();
	if(m_dsetMap!=NULL){
		delete m_dsetMap;
		m_dsetMap = NULL;
	}
	if(m_DB!=NULL){
		delete m_DB;
		m_DB = NULL;
	}
	m_iDatasets = 0;
	m_iGenes = 0;
	m_numThreads = 0;
	DEBUG = false;
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
	//check
	vector<char> vc;
	ushort tot = 0;
	CSeekTools::InitVector(vc, m_vecstrAllQuery.size(), (char) 0);
	map<ushort, vector< vector<string> > >::const_iterator ci =
		m_mapLoadTime.begin();
	for(; ci!=m_mapLoadTime.end(); ci++) tot+=(ci->second).size();
	if(tot!=m_vecstrAllQuery.size()){
		fprintf(stderr, "ERROR, size do not match\n");
		return false;
	}
	return true;
}

bool CSeekCentral::Initialize(const char *gene, const char *quant,
	const char *dset, const char *search_dset,
	const char *query, const char *platform, const char *db,
	const char *prep, const bool &useNibble, const ushort &num_db,
	const ushort &buffer, const bool &bSubtractAvg,
	const bool &bSubtractPlatformAvg, const bool &bDividePlatformStdev,
	const bool &bLogit){

	m_maxNumDB = buffer;
	m_numThreads = 8;
	omp_set_num_threads(m_numThreads);

	m_bSubtractGeneAvg = bSubtractAvg;
	m_bSubtractPlatformAvg = bSubtractPlatformAvg;
	m_bDividePlatformStdev = bDividePlatformStdev;
	m_bLogit = bLogit;

	vector<string> vecstrGeneID;
	if(!CSeekTools::ReadListTwoColumns(gene, vecstrGeneID, m_vecstrGenes))
		return false;

	CSeekTools::ReadQuantFile(quant, m_quant);
	m_DB = new CDatabase(useNibble);

	vector<string> vecstrDP;
	if(!CSeekTools::ReadListTwoColumns(dset, m_vecstrDatasets, vecstrDP))
		return false;

	vector<string> vecstrSDP;
	if(!CSeekTools::ReadListTwoColumns(search_dset,
		m_vecstrSearchDatasets, vecstrSDP))
		return false;

	map<string, ushort> mapstrintDataset;
	ushort i;
	for(i=0; i<m_vecstrDatasets.size(); i++){
		m_mapstrstrDatasetPlatform[m_vecstrDatasets[i]] = vecstrDP[i];
		mapstrintDataset[m_vecstrDatasets[i]] = i;
	}

	m_dsetMap = new CSeekIntIntMap(m_vecstrDatasets.size());
	for(i=0; i<m_vecstrSearchDatasets.size(); i++)
		m_dsetMap->Add(mapstrintDataset[m_vecstrSearchDatasets[i]]);

	if(!CSeekTools::ReadMultipleQueries(query, m_vecstrAllQuery))
		return false;

	vector<string> vecstrPlatforms;
	CSeekTools::ReadPlatforms(platform, m_vp, vecstrPlatforms,
		m_mapstriPlatform);

	m_iDatasets = m_vecstrDatasets.size();
	m_iGenes = m_vecstrGenes.size();

	m_DB->Open(db, m_vecstrGenes, m_iDatasets, num_db);

	CSeekTools::LoadDatabase(*m_DB, prep, m_vecstrDatasets,
		m_mapstrstrDatasetPlatform, m_mapstriPlatform, m_vp, m_vc);

	if(!CalculateRestart()) return false;
	//CSeekTools::ReadDatabaselets(*m_DB, m_vecstrAllQuery, m_vc);

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

bool CSeekCentral::PrepareOneQuery(CSeekQuery &query){
	vector<ushort> &queryGenes = query.GetQuery();
	ushort j;
	const vector<ushort> &allRDatasets = m_dsetMap->GetAllReverse();
	ushort iSearchDatasets = m_dsetMap->GetNumSet();
	ushort iQuery = queryGenes.size();

	for(j=0; j<iSearchDatasets; j++)
		m_vc[allRDatasets[j]]->InitializeQuery(queryGenes);

	m_rData = new ushort**[m_numThreads];
	for(j=0; j<m_numThreads; j++)
		m_rData[j] = CSeekTools::Init2DArray(m_iGenes, iQuery, (ushort)0);

	m_master_rank_threads =
		CSeekTools::Init2DArray(m_numThreads, m_iGenes, (float) 0);
	m_sum_weight_threads =
		CSeekTools::Init2DArray(m_numThreads, m_iGenes, (float) 0);
	m_counts_threads =
		CSeekTools::Init2DArray(m_numThreads, m_iGenes, (ushort) 0);
	m_rank_normal_threads = new vector<ushort>[m_numThreads];
	m_rank_threads = new vector<ushort>[m_numThreads];

	for(j=0; j<m_numThreads; j++){
		m_rank_normal_threads[j].resize(m_iGenes);
		m_rank_threads[j].resize(m_iGenes);
	}

	CSeekTools::InitVector(m_master_rank, m_iGenes, (float) 0);
	CSeekTools::InitVector(m_sum_weight, m_iGenes, (float) 0);
	CSeekTools::InitVector(m_counts, m_iGenes, (ushort) 0);
	CSeekTools::InitVector(m_weight, m_iDatasets, (float) 0);

	return true;
}

bool CSeekCentral::PostSearch(){
	ushort j, k;
	bool DEBUG = false;
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

	for(j=0; j<m_numThreads; j++) CSeekTools::Free2DArray(m_rData[j]);
	delete[] m_rData;
	delete[] m_rank_normal_threads;

	if(DEBUG) fprintf(stderr, "Aggregating genes\n");
	const vector<ushort> &allRDatasets = m_dsetMap->GetAllReverse();
	ushort iSearchDatasets = m_dsetMap->GetNumSet();

	for(j=0; j<m_iGenes; j++){
		if(m_counts[j]<(int)(0.5*iSearchDatasets))
			m_master_rank[j] = -320;
		else if(m_sum_weight[j]==0)
			m_master_rank[j] = -320;
		else
			m_master_rank[j] =
				(m_master_rank[j] / m_sum_weight[j] - 320) / 100.0;
		if(DEBUG) fprintf(stderr, "Gene %d %.5f\n", j, m_master_rank[j]);
	}
	return true;
}

bool CSeekCentral::Sort(){
	if(DEBUG) fprintf(stderr, "Sorting genes\n");
	m_final.clear();
	m_final.resize(m_iGenes);
	ushort j;
	for(j=0; j<m_iGenes; j++){
		//fprintf(stderr, "%d %s\n", j, DB.GetGene((size_t) j).c_str());
		m_final[j].i = j;
		m_final[j].f = m_master_rank[j];
	}
	if(DEBUG) fprintf(stderr, "Begin Sorting genes\n");
	sort(m_final.begin(), m_final.end());
	return true;
}

bool CSeekCentral::Display(CSeekQuery &query){
	if(DEBUG) fprintf(stderr, "Results:\n");
	ushort jj, ii;
	vector<char> &cQuery = query.GetQueryPresence();
	for(ii=0, jj=0; jj<500; ii++){
		if(cQuery[m_final[ii].i]==1) continue;
		fprintf(stderr, "%s %.5f\n",
			m_DB->GetGene((size_t)m_final[ii].i).c_str(), m_final[ii].f);
		jj++;
	}
	return true;
}

bool CSeekCentral::Write(const ushort &i){
	char acBuffer[1024];
	sprintf(acBuffer, "results/%d.query", i);
	CSeekTools::WriteArrayText(acBuffer, m_vecstrAllQuery[i]);
	sprintf(acBuffer, "results/%d.dweight", i);
	CSeekTools::WriteArray(acBuffer, m_weight);
	sprintf(acBuffer, "results/%d.gscore", i);
	CSeekTools::WriteArray(acBuffer, m_master_rank);
	return true;
}

bool CSeekCentral::Common(enum SearchMode &sm,
		gsl_rng *rnd, const enum PartitionMode *PART_M,
		const ushort *FOLD, const float *RATE){

	ushort i, j, d, dd;
	const vector<ushort> &allRDatasets = m_dsetMap->GetAllReverse();
	ushort iSearchDatasets = m_dsetMap->GetNumSet();

	for(i=0; i<m_vecstrAllQuery.size(); i++){
		if(m_mapLoadTime.find(i)!=m_mapLoadTime.end()){
			CSeekTools::ReadDatabaselets(*m_DB, m_mapLoadTime[i], m_vc);
		}

		CSeekQuery query;
		PrepareQuery(m_vecstrAllQuery[i], query);
		PrepareOneQuery(query);
		ushort iQuery = query.GetQuery().size();

		if(sm==CV){
			query.CreateCVPartitions(rnd, *PART_M, *FOLD);
		}

		#pragma omp parallel for \
		shared(allRDatasets, query) \
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

			if(mapQ==NULL ||mapQ->GetNumSet()==0){
				if(DEBUG) fprintf(stderr, "This dataset is skipped\n");
				continue;
			}

			vector<ushort> this_q;
			for(j=0; j<mapQ->GetNumSet(); j++)
				this_q.push_back(mapQ->GetReverse(j));


			if(DEBUG) fprintf(stderr, "Initializing %d\n", this_q.size());
			m_vc[d]->InitializeDataMatrix(m_rData[tid], m_quant, m_iGenes,
				iQuery, m_bSubtractGeneAvg, m_bSubtractPlatformAvg, m_bLogit);
			//m_bSubtractPlatformStdev is not used, it's assumed

			float w = -1;
			if(sm==CV){
				if(DEBUG) fprintf(stderr, "Weighting dataset\n");
				CSeekWeighter::CVWeighting(query, *m_vc[d], *RATE,
					&m_rank_threads[tid], false);
				if( (w = m_vc[d]->GetDatasetSumWeight())==-1){
					if(DEBUG) fprintf(stderr, "Bad weight\n");
					continue;
				}
			}else if(sm==EQUAL){
				w = 1.0;
			}

			if(DEBUG) fprintf(stderr, "Doing linear combination\n");

			CSeekWeighter::LinearCombine(m_rank_normal_threads[tid], this_q,
				*m_vc[d], false);

			m_vc[d]->DeleteQuery();

			if(DEBUG) fprintf(stderr,
				"Adding contribution of dataset to master ranking: %.5f\n", w);

			ushort iGeneSet = mapG->GetNumSet();
			const vector<ushort> &allRGenes = mapG->GetAllReverse();
			vector<ushort>::const_iterator iterR = allRGenes.begin();
			vector<ushort>::const_iterator endR = allRGenes.begin() + iGeneSet;
			vector<ushort> &Rank_Normal = m_rank_normal_threads[tid];
			float* Master_Rank = &m_master_rank_threads[tid][0];
			float* Sum_Weight = &m_sum_weight_threads[tid][0];
			ushort* Counts = &m_counts_threads[tid][0];

			for(; iterR!=endR; iterR++){
				if(Rank_Normal[*iterR]==0) continue;
				Master_Rank[*iterR] += (float) Rank_Normal[*iterR] * w;
				Sum_Weight[*iterR] += w;
				Counts[*iterR]++;
			}

			m_weight[d] = w;
		}
		//omp finishes
		PostSearch();
		Sort();
		//Display(query);
		fprintf(stderr, "Done search\n"); system("date +%s%N 1>&2");
		Write(i);
	}

	return true;
}


bool CSeekCentral::EqualWeightSearch(){
	enum SearchMode sm = EQUAL;
	return CSeekCentral::Common(sm);
}

bool CSeekCentral::CVSearch(gsl_rng *rnd, const enum PartitionMode &PART_M,
	const ushort &FOLD, const float &RATE){
	enum SearchMode sm = CV;
	return CSeekCentral::Common(sm, rnd, &PART_M, &FOLD, &RATE);
}

bool CSeekCentral::Destruct(){
	ushort j;
	for(j=0; j<m_iDatasets; j++) m_vc[j]->DeleteQueryBlock();
	for(j=0; j<m_iDatasets; j++) delete m_vc[j];
	return true;
}

}

