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
 * \brief Representation of a microarray dataset that is used by Seek
 *
 * A \c CSeekDataset encapsulates the following information about the dataset:
 *
 * \li The gene-gene \a correlation matrix
 * \li Each gene's expression variance
 * \li Each gene's average \a correlation
 * \li The genes in the query that are present in the dataset
 * \li The genes in the genome that are present in the dataset
 * \li The platform which the dataset belongs to
 * \li The weight of the dataset that is assigned by the search algorithm
 *
 * This dataset structure is designed to be used by Seek.
 * \remarks
 * The word \a correlation refers to the standardized z-scores of Pearson correlations, which
 * is derived from a 2-step process:
 * \f[f(x,y)=\frac{1}{2}ln\frac{1+p(x,y)}{1-p(x,y)}\f]
 * where \f$p(x,y)\f$ is the Pearson correlation, \f$f(x,y)\f$ is the Fisher's transformed score.
 * \f[z(x,y)=\frac{f(x,y) - \bar{f}}{\sigma_{f}}\f]
 * where \f$z(x,y)\f$ is the z-score, \f$\bar{f}\f$ is the mean, and \f$\sigma_{f}\f$ is
 * the standard deviation.
 * 
 * From here on, \a correlation always refers to the above z-score definition.
 */

class CSeekDataset{
public:

	/*!
	 * \brief
	 * Distance measure
	 */
	enum DistanceMeasure{
		CORRELATION = 0,
		Z_SCORE = CORRELATION + 1
	};

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
	 * of the global gene-gene Pearson distribution for this dataset.
	 */
	bool ReadDatasetAverageStdev(const string &);

	/*!
	 * \brief
	 * Read the gene average \a correlation file \c *.gavg
	 *
	 * \param strFileName The file name
	 *
	 * The \c *.gavg is an array that stores the average \a correlation of each
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
	 * Initialize the gene-gene \a correlation matrix
	 *
	 * \param rD A two-dimensional array storing the discretized gene-gene \a correlations
	 * \param quant The discretization function
	 * \param iRows The number of rows for the \a correlation matrix
	 * \param iColumns The number of columns for the \a correlation matrix
	 * \param bSubtractAvg If true, subtract the \a correlation by the dataset average
	 * \param bNormPlatform If true, subtract the \a correlation by the platform average and divide by standard deviation
	 * \param logit If true, apply the logit transform on \a correlations
	 * \param dist_measure Distance measure: z-score or correlations
	 * \param cutoff Apply a hard cutoff on \a correlations
	 * \param bRandom If true, shuffle the \a correlation vector
	 * \param rand The random generator for the shuffling operation above
	 * \remarks
	 * The discretized \a correlation in the matrix \c rD is bounded by 0 to 255 (the limit of
	 * \c unsigned \c char). The parameter \c quant specifies how a \a correlation is
	 * discretized. For example, if the \c quant has 5 bins:
	 * \code 
	 * [0, 1, 2, 3, 4]
	 * \endcode
	 * Then if a \a correlation is 2.5, the discretized value would be 2.
	 */
	bool InitializeDataMatrix(ushort**, const vector<float> &,
		const ushort&, const ushort&, const bool=true, 
		const bool=false, const bool=false,
		const enum DistanceMeasure=Z_SCORE,
		const float cutoff=-1.0*CMeta::GetNaN(), 
		const bool=false, gsl_rng *rand=NULL);

	/*!
	 * \brief
	 * Copy constructor
	 * \param src A given dataset
	 */
	bool Copy(CSeekDataset *);

	/*!
	 * \brief
	 * Get the gene-gene \a correlation matrix
	 *
	 * \return A two-dimensional array of type \c ushort. Note that the
	 * \a correlation has been scaled to a integer range from 0 to 640.
	 * See CSeekDataset::InitializeDataMatrix.
	 *
	 */
	ushort** GetDataMatrix();

	/*!
	 * \brief
	 * Get the gene-gene \a correlation matrix
	 *
	 * \return A two-dimensional array of type \c unsigned \c char**.
	 */
	unsigned char** GetMatrix();

	/*!
	 * \brief Get the genome presence map
	 * \return The genome presence map
	 */
	CSeekIntIntMap* GetGeneMap();

	/*!
	 * \brief Get the query-block presence map
	 * \return The query-block presence map
	 */
	CSeekIntIntMap* GetDBMap();

	/*!
	 * \brief Get the query presence map
	 * \return The query presence map
	 */
	CSeekIntIntMap* GetQueryMap();

	/*!
	 * \brief Get the query genes
	 * \return A vector of queries
	 */
	const vector<ushort>& GetQuery() const;

	/*!
	 * \brief Get the query gene indices
	 * \return A vector of query gene indices
	 */
	const vector<ushort>& GetQueryIndex() const;

	/*!
	 * \brief Get the gene expression variance vector
	 * \return The variance vector
	 */
	float GetGeneVariance(const ushort&) const;
	/*!
	 * \brief Get the gene average \a correlation vector
	 * \return The average \a correlation vector
	 */
	float GetGeneAverage(const ushort&) const;
	/*!
	 * \brief Get the mean of the global gene-gene Pearson distribution
	 * \return The mean Pearson for the dataset
	 */
	float GetDatasetAverage() const;
	/*!
	 * \brief Get the standard deviation of the global gene-gene Pearson distribution
	 * \return The standard deviation of the Pearson distribution
	 */
	float GetDatasetStdev() const;
	/*!
	 * \brief Get the genome size
	 * \return The genome size
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
	 * \return A vector of cross-validation scores
	 */
	const vector<float>& GetCVWeight() const;

	/*!
	 * \brief Get the dataset weight
	 * \return The dataset weight
	 */
	float GetDatasetSumWeight();

	/*!
	 * \brief Set the platform
	 * \param cp The platform
	 */
	void SetPlatform(CSeekPlatform &);
	/*!
	 * \brief Get the platform
	 * \return The platform of this dataset
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
