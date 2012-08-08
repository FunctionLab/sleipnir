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
	m_bOutputText = false;
	DEBUG = false;
	m_output_dir = "";

}

CSeekCentral::~CSeekCentral(){
	m_vecstrGenes.clear();
	m_vecstrDatasets.clear();
	m_vecstrSearchDatasets.clear();
	m_mapstrstrDatasetPlatform.clear();

	ushort i;
	for(i=0; i<m_vc.size(); i++){
		if(m_vc[i]==NULL) continue;
		delete m_vc[i];
	}
	m_vc.clear();
	m_quant.clear();

	for(i=0; i<m_searchdsetMap.size(); i++){
		if(m_searchdsetMap[i]==NULL) continue;
		delete m_searchdsetMap[i];
	}
	m_searchdsetMap.clear();

	//m_rData, m_master_rank_threads,
	//m_sum_weight_threads, m_counts_threads
	//m_rank_normal_threads, and m_rank_threads
	//should be already freed

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
		delete m_DB;
		m_DB = NULL;
	}
	m_iDatasets = 0;
	m_iGenes = 0;
	m_numThreads = 0;
	m_mapLoadTime.clear();
	m_output_dir = "";
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

bool CSeekCentral::Initialize(const char *gene, const char *quant,
	const char *dset, const char *search_dset,
	const char *query, const char *platform, const char *db,
	const char *prep, const bool &useNibble, const ushort &num_db,
	const ushort &buffer, const char *output_dir, const bool &to_output_text,
	const bool &bSubtractAvg,
	const bool &bSubtractPlatformAvg, const bool &bDividePlatformStdev,
	const bool &bLogit, const float &fCutOff, const float &fPercentRequired){

	m_output_dir = output_dir;
	m_maxNumDB = buffer;
	m_numThreads = 8;
	m_fScoreCutOff = fCutOff;
	m_fPercentQueryAfterScoreCutOff = fPercentRequired;
	ushort i, j;

	omp_set_num_threads(m_numThreads);

	m_bOutputText = to_output_text;
	m_bSubtractGeneAvg = bSubtractAvg;
	m_bSubtractPlatformAvg = bSubtractPlatformAvg;
	m_bDividePlatformStdev = bDividePlatformStdev;
	m_bLogit = bLogit;

	//read genes
	vector<string> vecstrGeneID;
	if(!CSeekTools::ReadListTwoColumns(gene, vecstrGeneID, m_vecstrGenes))
		return false;

	CSeekTools::ReadQuantFile(quant, m_quant);
	m_DB = new CDatabase(useNibble);

	//read datasets
	vector<string> vecstrDP;
	if(!CSeekTools::ReadListTwoColumns(dset, m_vecstrDatasets, vecstrDP))
		return false;

	map<string, ushort> mapstrintDataset;
	for(i=0; i<m_vecstrDatasets.size(); i++){
		m_mapstrstrDatasetPlatform[m_vecstrDatasets[i]] = vecstrDP[i];
		mapstrintDataset[m_vecstrDatasets[i]] = i;
	}

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
				mapstrintDataset[m_vecstrSearchDatasets[i][j]]);
	}

	vector<string> vecstrPlatforms;
	CSeekTools::ReadPlatforms(platform, m_vp, vecstrPlatforms,
		m_mapstriPlatform);

	m_iDatasets = m_vecstrDatasets.size();
	m_iGenes = m_vecstrGenes.size();

	m_DB->Open(db, m_vecstrGenes, m_iDatasets, num_db);
	CSeekTools::LoadDatabase(*m_DB, prep, m_vecstrDatasets,
		m_mapstrstrDatasetPlatform, m_mapstriPlatform, m_vp, m_vc);

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
	}

	CSeekTools::InitVector(m_master_rank, m_iGenes, (float) 0);
	CSeekTools::InitVector(m_sum_weight, m_iGenes, (float) 0);
	CSeekTools::InitVector(m_counts, m_iGenes, (ushort) 0);
	CSeekTools::InitVector(weight, m_iDatasets, (float)0);

	return true;
}

bool CSeekCentral::AggregateThreads(){
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
	for(j=0; j<m_numThreads; j++){
		m_rank_normal_threads[j].clear();
		m_rank_threads[j].clear();
	}
	delete[] m_rank_normal_threads;
	delete[] m_rank_threads;
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
		else
			m_master_rank[j] =
				(m_master_rank[j] / m_sum_weight[j] - 320) / 100.0;
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
	char acBuffer[1024];
	sprintf(acBuffer, "%s/%d.query", m_output_dir.c_str(), i);
	CSeekTools::WriteArrayText(acBuffer, m_vecstrAllQuery[i]);
	sprintf(acBuffer, "%s/%d.dweight", m_output_dir.c_str(), i);
	CSeekTools::WriteArray(acBuffer, m_weight[i]);
	sprintf(acBuffer, "%s/%d.gscore", m_output_dir.c_str(), i);
	CSeekTools::WriteArray(acBuffer, m_master_rank);
	if(m_bOutputText){
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
		for(j=0; j<200; j++){
			if(w[j].f==0) continue;
			vecOutput[0].push_back(m_vecstrDatasets[w[j].i]);
		}
		vector<AResultFloat> wd;
		Sort(wd);
		for(j=0; j<2000; j++){
			if(wd[j].f==-320) continue;
			vecOutput[1].push_back(m_vecstrGenes[wd[j].i]);
		}
		CSeekTools::Write2DArrayText(acBuffer, vecOutput);
	}
	return true;
}



bool CSeekCentral::Common(enum SearchMode &sm,
	gsl_rng *rnd, const enum PartitionMode *PART_M,
	const ushort *FOLD, const float *RATE,
	const vector< vector<float> > *providedWeight,
	const vector< vector<string> > *newGoldStd){

	ushort i, j, d, dd;

	m_Query.resize(m_vecstrAllQuery.size());
	m_weight.resize(m_vecstrAllQuery.size());
	m_final.resize(m_vecstrAllQuery.size());

	for(i=0; i<m_vecstrAllQuery.size(); i++){
		if(m_mapLoadTime.find(i)!=m_mapLoadTime.end())
			CSeekTools::ReadDatabaselets(*m_DB, m_mapLoadTime[i], m_vc);

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

		//for CV_CUSTOM
		CSeekQuery customGoldStd;
		if(sm==CV_CUSTOM)
			PrepareQuery((*newGoldStd)[i], customGoldStd);

		if(sm==CV || sm==CV_CUSTOM)
			query.CreateCVPartitions(rnd, *PART_M, *FOLD);

		if(sm==ORDER_STATISTICS)
			m_rank_d = CSeekTools::Init2DArray(iSearchDatasets, m_iGenes,
				(ushort) 0);

		//fprintf(stderr, "2 %lu\n", CMeta::GetMemoryUsage());

		#pragma omp parallel for \
		shared(allRDatasets, query, customGoldStd) \
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

			if(DEBUG) fprintf(stderr, "Initializing %d\n",
				(int) this_q.size());
			m_vc[d]->InitializeDataMatrix(m_rData[tid], m_quant, m_iGenes,
				iQuery, m_bSubtractGeneAvg, m_bSubtractPlatformAvg, m_bLogit,
				m_fScoreCutOff);
			//m_bSubtractPlatformStdev is not used, it's assumed

			float w = -1;
			if(sm==CV || sm==CV_CUSTOM){
				if(DEBUG) fprintf(stderr, "Weighting dataset\n");
				if(sm==CV)
					CSeekWeighter::CVWeighting(query, *m_vc[d], *RATE,
						m_fPercentQueryAfterScoreCutOff,
						&m_rank_threads[tid]);
				else
					CSeekWeighter::CVWeighting(query, *m_vc[d], *RATE,
						m_fPercentQueryAfterScoreCutOff,
						&m_rank_threads[tid], &customGoldStd);

				if( (w = m_vc[d]->GetDatasetSumWeight())==-1){
					if(DEBUG) fprintf(stderr, "Bad weight\n");
					continue;
				}
			}
			else if(sm==EQUAL || sm==ORDER_STATISTICS) w = 1.0;
			else if(sm==USE_WEIGHT) w = (*providedWeight)[i][d];

			if(DEBUG) fprintf(stderr, "Doing linear combination\n");

			const ushort MIN_REQUIRED = max((ushort) 1, (ushort) (
				m_fPercentQueryAfterScoreCutOff * this_q.size()));
			CSeekWeighter::LinearCombine(m_rank_normal_threads[tid], this_q,
				*m_vc[d], MIN_REQUIRED);

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

			if(sm==ORDER_STATISTICS)
				for(; iterR!=endR; iterR++){
					if(Rank_Normal[*iterR]==0) continue;
					m_rank_d[dd][*iterR] = Rank_Normal[*iterR];
					Counts[*iterR]++;
				}
			else
				for(; iterR!=endR; iterR++){
					if(Rank_Normal[*iterR]==0) continue;
					Master_Rank[*iterR] += (float) Rank_Normal[*iterR] * w;
					Sum_Weight[*iterR] += w;
					Counts[*iterR]++;
				}
			weight[d] = w;
		}
		//omp finishes
		//fprintf(stderr, "3 %lu\n", CMeta::GetMemoryUsage());
		for(j=0; j<iSearchDatasets; j++)
			m_vc[allRDatasets[j]]->DeleteQuery();

		for(j=0; j<m_numThreads; j++)
			CSeekTools::Free2DArray(m_rData[j]);
		delete[] m_rData;

		AggregateThreads();

		if(sm!=ORDER_STATISTICS){
			FilterResults(iSearchDatasets);
		}else{
			CSeekWeighter::OrderStatisticsRankAggregation(iSearchDatasets,
				m_iGenes, m_rank_d, m_counts, m_master_rank);
			CSeekTools::Free2DArray(m_rank_d);
			m_rank_d = NULL;
		}

		Sort(final);
		//Display(query, final);
		//fprintf(stderr, "4 %lu\n", CMeta::GetMemoryUsage());

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

/*	perform CVSearch, except that the weighting is based not on co-expression
	of query genes, but based on similarity of query genes to some custom gold
	standard gene-set */
bool CSeekCentral::CVCustomSearch(const vector< vector<string> > &newGoldStd,
	gsl_rng *rnd, const enum PartitionMode &PART_M,
	const ushort &FOLD, const float &RATE){
	enum SearchMode sm = CV_CUSTOM;
	return CSeekCentral::Common(sm, rnd, &PART_M, &FOLD, &RATE,
		NULL, &newGoldStd);
}

bool CSeekCentral::WeightSearch(const vector< vector<float> > &weights){
	enum SearchMode sm = USE_WEIGHT;
	return CSeekCentral::Common(sm, NULL, NULL, NULL, NULL, &weights);
}

bool CSeekCentral::OrderStatistics(){
	enum SearchMode sm = ORDER_STATISTICS;
	return CSeekCentral::Common(sm);
}

bool CSeekCentral::SingleGeneMetaCorrelation(){
	enum SearchMode sm = SINGLE_GENE_META;
	return ;
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
	for(j=0; j<m_iDatasets; j++) m_vc[j]->DeleteQueryBlock();
	for(j=0; j<m_iDatasets; j++){
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

