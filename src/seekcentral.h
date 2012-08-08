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
#ifndef SEEKCENTRAL_H
#define SEEKCENTRAL_H

#include "seekbasic.h"
#include "seekdataset.h"
#include "seekplatform.h"
#include "seekmap.h"
#include "seekreader.h"
#include "seekquery.h"
#include "seekevaluate.h"
#include "database.h"
#include "datapair.h"
#include "seekweight.h"

namespace Sleipnir {

enum SearchMode{
	CV=0, EQUAL=1, USE_WEIGHT=2, CV_CUSTOM=3, ORDER_STATISTICS=4,
	SINGLE_GENE_META=5
};

class CSeekCentral{
public:
	CSeekCentral();
	~CSeekCentral();

	bool Initialize(const char *gene, const char *quant,
		const char *dset, const char *search_dset,
		const char *query, const char *platform, const char* db,
		const char *prep, const bool &useNibble, const ushort &num_db,
		const ushort &, const char*, const bool&,
		const bool&, const bool&, const bool&, const bool&,
		const float&, const float&);

	bool CVSearch(gsl_rng*, const enum PartitionMode&, const ushort&, const float&);
	bool CVCustomSearch(const vector< vector<string> > &, gsl_rng*,
		const enum PartitionMode&, const ushort&, const float&);
	bool EqualWeightSearch();
	bool WeightSearch(const vector<vector<float> >&);
	bool SingleGeneMetaCorrelation();
	bool VarianceWeightSearch();

	bool OrderStatistics();

	bool Common(enum SearchMode&, gsl_rng* = NULL, const enum PartitionMode* = NULL,
		const ushort* = NULL, const float* = NULL,
		const vector< vector<float> >* = NULL,
		const vector< vector<string> >* = NULL);

	bool Destruct();
	bool PrepareQuery(const vector<string>&, CSeekQuery&);
	bool CalculateRestart();

	bool PrepareOneQuery(CSeekQuery &, CSeekIntIntMap &, vector<float>&);
	bool AggregateThreads();
	bool FilterResults(const ushort &);
	bool Sort(vector<AResultFloat> &);
	bool Write(const ushort &);
	bool Display(CSeekQuery &, vector<AResultFloat>&);
	const vector< vector<AResultFloat> >& GetAllResult()const;
	const vector<CSeekQuery>& GetAllQuery() const;
	ushort GetGene(const string &strGene) const;
	string GetGene(const ushort &geneID) const;

	const vector<vector<float> > &GetAllWeight() const;

private:
	vector<string> m_vecstrGenes;

	/* Dataset */
	vector<string> m_vecstrDatasets;
	map<string, string> m_mapstrstrDatasetPlatform;
	vector<vector<string> > m_vecstrSearchDatasets;
	vector<CSeekIntIntMap*> m_searchdsetMap;
	vector<CSeekDataset*> m_vc;
	vector<float> m_quant;
	bool m_bSubtractGeneAvg;
	bool m_bSubtractPlatformAvg;
	bool m_bDividePlatformStdev;
	bool m_bLogit;
	bool m_bOutputText;

	ushort ***m_rData;

	/* multi-threaded programming */
	float **m_master_rank_threads;
	float **m_sum_weight_threads;
	ushort **m_counts_threads;
	vector<ushort> *m_rank_normal_threads;
	vector<ushort> *m_rank_threads;

	/* Essential search results */
	vector<float> m_master_rank;
	vector<float> m_sum_weight;
	vector<ushort> m_counts;

	/* Holds results for all queries */
	vector< vector<float> > m_weight;
	vector< vector<AResultFloat> > m_final;

	/* Query */
	vector< vector<string> > m_vecstrAllQuery;
	vector<CSeekQuery> m_Query;

	/* Platform */
	vector<CSeekPlatform> m_vp;
	map<string, ushort> m_mapstriPlatform;
	vector<string> m_vecstrPlatform;

	CDatabase *m_DB;

	size_t m_iDatasets;
	size_t m_iGenes;
	ushort m_numThreads;

	ushort m_maxNumDB;
	map<ushort, vector< vector<string> > > m_mapLoadTime;
	bool DEBUG;

	string m_output_dir;
	float m_fScoreCutOff;
	float m_fPercentQueryAfterScoreCutOff;

	/* for order statistics, a datasets-by-genes matrix */
	ushort **m_rank_d;

};

}
#endif
