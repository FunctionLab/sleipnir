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

#include "seekbasic.h"
#include "seekmap.h"
#include "datapair.h"
#include "seekplatform.h"

namespace Sleipnir {

class CSeekDataset{
public:
	/*!
	 * \brief
	 * Default constructor
	 *
	 * Initialize the data-set weight, the gene-centric average, variance,
	 * and presence vectors, and prepare the data-matrix that will store
	 * the pairwise gene correlation (for group of genes in the DB)
	 * in this data-set
	 */
	CSeekDataset();

	/*!
	 * \brief
	 * Destruct the data-set instance
	 *
	 * Reset the data-matrix, and all gene-centric vectors
	 */
	~CSeekDataset();

	bool ReadDatasetAverageStdev(const string &);

	/*!
	 * \brief
	 * Given the gene-average file-name, read the gene averages in this
	 * data-set.
	 *
	 * \param strFileName
	 * File-name of the gene-average file (.gavg extension)
	 *
	 * The .gavg format lists the average of a gene (to all other
	 * genes in this data-set), in float binary format. The gene-averages
	 * are ordered by the gene-ID.
	 */
	bool ReadGeneAverage(const string &);

	/*!
	 * \brief
	 * Given the gene-variance file-name, read the gene variances in this
	 * data-set
	 *
	 * \param strFileName
	 * File-name of the gene-variance file (.gvar extension)
	 *
	 * The .gvar format lists the variance of a gene (calculated from
	 * expressions across all conditions in this data-set), in float binary
	 * format. The gene-variances are ordered by the gene-ID.
	 */
	bool ReadGeneVariance(const string &);

	/*!
	 * \brief
	 * Given the gene-presence file-name, read the gene-presence vector in
	 * this data-set
	 *
	 * \param strFileName
	 * File-name of the gene-presence file (.gpres extension)
	 *
	 * The .gpres format lists the presence/absence value (0/1) of the
	 * gene in the data-set (based on what genes are present in the PCL
	 * file), in char binary format. The gene-presences are ordered by
	 * gene-ID.
	 */
	bool ReadGenePresence(const string &);

	/*!
	 * \brief
	 * Initialize the gene-map for the data-set
	 *
	 * Based on the gene-presence vector already loaded, initialize a
	 * CSeekIntIntMap instance for the genes in this data-set.
	 * This map will map the gene-ID to the position of the gene in the
	 * data-set (forward mapping), and vice-versa (reverse mapping).
	 * For example, let's suppose gene-ID 10246 appears in position 3
	 * in the data-set. Then: the forward mapping maps 10246 to 3, and
	 * reverse mapping maps 3 to 10246. Because the data for only the
	 * gene-ID's available in the data-set are stored, the map tells us
	 * how to find a gene-ID in the data-set.
	 */
	bool InitializeGeneMap();

	/*!
	 * \brief
	 * Initialize a map for a specific query
	 *
	 * \param query
	 * Genes in the query
	 *
	 * Suppose the query is a strict subset of genes in the queryBlock.
	 * Create a CSeekIntIntMap based on the overlap between query and
	 * queryBlock.
	 */
	bool InitializeQuery(const vector<ushort> &);

	/*!
	 * \brief
	 * Initialize a map for the query-block genes
	 *
	 * \param queryBlock
	 * Vector of genes in the query-block
	 *
	 * Create a CSeekIntIntMap for the query-block genes, given also what
	 * genes are present in the data-set. If a gene is present in the
	 * data-set and is also a gene in queryBlock, then add it to the map.
	 * Assumes presence-vector is loaded InitializeGeneMap() has been called.
	 */
	bool InitializeQueryBlock(const vector<ushort> &);

	/*!
	 * \brief
	 * Delete query
	 *
	 * Reset all query-related data structures, such as data-set weight,
	 * iQuerySize, queryMap, etc.
	 */
	bool DeleteQuery();

	/*!
	 * \brief
	 * Delete query block
	 *
	 * Reset all query-block related data.
	 */
	bool DeleteQueryBlock();
	bool DeleteR();

	/*!
	 * \brief
	 * Initialize the matrix storing the data
	 */
	bool InitializeDataMatrix(ushort**, const vector<float> &,
		const ushort&, const ushort&, const bool=true, const bool=true,
		const bool=false, const bool=false,
		const float cutoff=-1.0*CMeta::GetNaN(), 
		const bool=false, gsl_rng *rand=NULL);

	//Copy CSeekDataset object, mainly for SeekServer
	bool Copy(CSeekDataset *);

	ushort** GetDataMatrix();

	unsigned char** GetMatrix();
	CSeekIntIntMap* GetGeneMap();
	CSeekIntIntMap* GetDBMap();
	CSeekIntIntMap* GetQueryMap();

	const vector<ushort>& GetQuery() const;
	const vector<ushort>& GetQueryIndex() const;

	float GetGeneVariance(const ushort&) const;
	float GetGeneAverage(const ushort&) const;
	float GetDatasetAverage() const;
	float GetDatasetStdev() const;
	ushort GetNumGenes() const;
	bool InitializeCVWeight(const ushort&);
	bool SetCVWeight(const ushort&, const float&);
	float GetDatasetSumWeight();
	void SetPlatform(CSeekPlatform &);
	CSeekPlatform& GetPlatform() const;

private:
	CSeekPlatform *platform;
	vector<float> geneAverage;
	vector<float> geneVariance;
	vector<char> genePresence;
	CSeekIntIntMap *geneMap;

	/* previously known as sinfo file */
	float m_fDsetAverage;
	float m_fDsetStdev;

	CSeekIntIntMap *dbMap;
	CSeekIntIntMap *queryMap;
	vector<ushort> query;
	vector<ushort> queryIndex;

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
