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
/*!
 * \brief A microarray dataset structure
 *
 * A \c CSeekDataset encapsulates the following information about the dataset:
 *
 * \li The gene-gene correlation matrix
 * \li Each gene's expression variance
 * \li Each gene's average correlation
 * \li The genes in the query that are present in the dataset
 * \li The genes in the genome that are present in the dataset
 * \li The platform which the dataset belongs to
 * \li The weight of the dataset that is assigned by the search algorithm
 *
 * This dataset structure is designed to be used by Seek.
 */

class CSeekDataset{
public:
	/*!
	 * \brief
	 * Constructor
	 */
	CSeekDataset();

	/*!
	 * \brief
	 * Destructor
	 */
	~CSeekDataset();

	/*!
	 * \brief
	 * Read the \c *.sinfo file
	 *
	 * \param strFileName The file name
	 *
	 * The \c *.sinfo file contains the mean and the standard deviation
	 * of the global gene-gene correlation distribution for this dataset.
	 */
	bool ReadDatasetAverageStdev(const string &);

	/*!
	 * \brief
	 * Read the gene average correlation file \c *.gavg
	 *
	 * \param strFileName The file name
	 *
	 * The \c *.gavg is an array that stores the average correlation of each
	 * gene.
	 */
	bool ReadGeneAverage(const string &);

	/*!
	 * \brief
	 * Read the gene variance file \c *.gvar
	 *
	 * \param strFileName The file name
	 *
	 * The \c *.gvar file is an array that stores the expression variance of each
	 * gene.
	 */
	bool ReadGeneVariance(const string &);

	/*!
	 * \brief
	 * Read the gene presence file \c *.gpres
	 *
	 * \param strFileName The file name
	 *
	 * The \c *.gpres is a 2-value array that contains the presence (1), absence (0)
	 * status of genes.
	 */
	bool ReadGenePresence(const string &);

	/*!
	 * \brief
	 * Initialize the genome presence map
	 *
	 * Indicates which genes of the genome are present in the dataset.
	 */
	bool InitializeGeneMap();

	/*!
	 * \brief
	 * Initialize the query presence map
	 *
	 * \param query The query genes
	 *
	 * Indicates which query genes are present in the dataset.
	 */
	bool InitializeQuery(const vector<ushort> &);

	/*!
	 * \brief
	 * Initialize a presence map for a block of queries
	 *
	 * \param queryBlock A vector of queries
	 *
	 * Flattens all the queries into one vector that contains only the unique query genes, then
	 * constructs a presence map based on this vector.
	 */
	bool InitializeQueryBlock(const vector<ushort> &);

	/*!
	 * \brief
	 * Delete the query
	 *
	 * Resets all query-related data, such as dataset weight, query presence map, etc.
	 */
	bool DeleteQuery();

	/*!
	 * \brief
	 * Delete query block
	 *
	 * Resets all query-block related data.
	 */
	bool DeleteQueryBlock();

	/*!
	 * \brief
	 * Initialize the gene-gene correlation matrix
	 */
	bool InitializeDataMatrix(ushort**, const vector<float> &,
		const ushort&, const ushort&, const bool=true, const bool=true,
		const bool=false, const bool=false,
		const float cutoff=-1.0*CMeta::GetNaN(), 
		const bool=false, gsl_rng *rand=NULL);

	/*!
	 * \brief
	 * Copy constructor
	 */
	bool Copy(CSeekDataset *);

	/*!
	 * \brief
	 * Get the gene-gene correlation matrix
	 *
	 * A two-dimensional array of type \c ushort is returned. Note that the
	 * correlation has been scaled to a integer range from 0 to 640.
	 * See CSeekDataset::InitializeDataMatrix.
	 *
	 */
	ushort** GetDataMatrix();

	/*!
	 * \brief
	 * Get the gene-gene correlation matrix
	 *
	 * A two-dimensional array of type \c unsigned \c char** is returned.
	 */
	unsigned char** GetMatrix();

	/*!
	 * \brief Get the genome presence map
	 */
	CSeekIntIntMap* GetGeneMap();

	/*!
	 * \brief Get the query-block presence map
	 */
	CSeekIntIntMap* GetDBMap();

	/*!
	 * \brief Get the query presence map
	 */
	CSeekIntIntMap* GetQueryMap();

	/*!
	 * \brief Get the query genes
	 */
	const vector<ushort>& GetQuery() const;

	/*!
	 * \brief Get the query gene indices
	 */
	const vector<ushort>& GetQueryIndex() const;

	/*!
	 * \brief Get the gene expression variance vector
	 */
	float GetGeneVariance(const ushort&) const;
	/*!
	 * \brief Get the gene average correlation vector
	 */
	float GetGeneAverage(const ushort&) const;
	/*!
	 * \brief Get the mean of the global gene-gene correlation distribution
	 */
	float GetDatasetAverage() const;
	/*!
	 * \brief Get the standard deviation of the global gene-gene correlation distribution
	 */
	float GetDatasetStdev() const;
	/*!
	 * \brief Get the genome size
	 */
	ushort GetNumGenes() const;

	/*!
	 * \brief Initialize the weight of the dataset
	 *
	 * \param i The number of cross-validations
	 *
	 * Initializes the total dataset weight, and the score of the
	 * individual cross-validation (CV) runs. 
	 */
	bool InitializeCVWeight(const ushort&);

	/*!
	 * \brief Set the score for a particular cross-validation
	 * \param i The index
	 * \param f The validation score
	 */
	bool SetCVWeight(const ushort&, const float&);

	/*!
	 * \brief Get the score for a particular cross-validation
	 * \param i The index
	 */
	float GetCVWeight(const ushort&);

	/*!
	 * \brief Get all the cross-validation scores
	 */
	const vector<float>& GetCVWeight() const;

	/*!
	 * \brief Get the dataset weight
	 */
	float GetDatasetSumWeight();

	/*!
	 * \brief Set the platform
	 */
	void SetPlatform(CSeekPlatform &);
	/*!
	 * \brief Get the platform
	 */
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
