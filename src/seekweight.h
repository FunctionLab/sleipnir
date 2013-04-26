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
#ifndef SEEKWEIGHT_H
#define SEEKWEIGHT_H

#include "seekbasic.h"
#include "seekreader.h"
#include "seekquery.h"
#include "seekevaluate.h"

namespace Sleipnir {

/*!
 * \brief
 * Provide functions to assign dataset weight using the query gene.
 *
 * For dataset weighting, one way is to use CSeekWeighter::CVWeighting. The CSeekWeighter::CVWeighting
 * uses a cross-validation (CV) framework, where it partitions the query and performs a search
 * instance on one sub-query, using the remainder of the queries as the evaluation of the search instance.
 *
 * The CSeekWeighter::OrderStatisticsRankAggregation is a rank-based technique described by Adler et al (2009). This combines
 * dataset weighting and dataset gene-ranking aggregation all into one step.
 *
 */
class CSeekWeighter{
public:

	/*!
	 * \brief
	 * Calculates for each gene the average \a correlation to all of the query genes in a dataset.
	 *
	 * \param rank
	 * A vector that stores the \a correlation of each gene to all of the query genes
	 *
	 * \param cv_query
	 * A vector that stores the query genes
	 *
	 * \param sDataset
	 * A dataset
	 *
	 * \param MIN_REQUIRED
	 * A ushort that specifies how many query genes are required to be present in a dataset.
	 * If not enough query genes are present, then the averaging is not performed.
	 *
	 * \param bSquareZ
	 * If true, square the \a correlation values before adding \a correlations.
	 *
	 * \remark
	 * The word \a correlations refer to z-scored, standardized Pearson correlations.
	 * The result is returned in the parameter \c rank.
	 *
	 */
	/*cv_query must be present in sDataset */
	static bool LinearCombine(vector<ushort> &rank,
		const vector<ushort> &cv_query, CSeekDataset &sDataset,
		const ushort &, const bool &);

	/*!
	 * \brief
	 * Cross-validates query-genes in a dataset
	 *
	 * \param sQuery
	 * The query and its partitions
	 *
	 * \param sDataset
	 * A dataset
	 *
	 * \param rate
	 * RBP parameter \a p
	 *
	 * \param percent_required
	 * Percentage of query genes required to be present in the dataset
	 *
	 * \param bSquareZ
	 * Whether or not to square \a correlations
	 *
	 * \param rrank
	 * Temporary vector storing intermediary \a correlations
	 *
	 * \param goldStd
	 * If a gold-standard gene-set is provided, use this to evaluate the retrieval of a cross-validation
	 *
	 * This performs multiple cross-validation runs to validate
	 * the query genes in retrieving themselves in a dataset. The sum of the evaluation of all the
	 * runs then becomes the dataset weight. For evaluation, we use <em>rank-biased precision (RBP)</em>.
	 * The parameter \a p needs to be provided for RBP; the default value is 0.99.
	 *
	 */
	static bool CVWeighting(CSeekQuery &sQuery, CSeekDataset &sDataset,
		const float &rate, const float &percent_required, const bool &bsquareZ,
		vector<ushort> *rrank, const CSeekQuery *goldStd = NULL);

	/*!
	 * \brief
	 * Performs OrderStatisticsAggregation, also known as the MEM algorithm
	 * 
	 * \param iDatasets
	 * The number of datasets
	 * 
	 * \param iGenes
	 * The number of genes
	 * 
	 * \param rank_d
	 * Two-dimensional vectors storing correlation-ranks to the query genes. 
	 * First dimension: datasets. Second dimension: genes.
	 * 
	 * \param counts
	 * A vector storing the count of datasets for each gene
	 * 
	 * \param master_rank
	 * A vector storing the integrated gene-score
	 * 
	 * \param numThreads
	 * The number of threads to be used (in a parallel setup)
	 *
	 * \c rank_d needs to be prepared as follows: a <tt>correlation rank</tt> vector is obtained from sorting Pearson correlations 
	 * in a dataset, and then it is normalized by <tt>(rank of correlation) / (number of genes)</tt>. The result is stored
	 * in \c rank_d.
	 *
	 * Afterward, for each gene \a g, the algorithm compares this gene's \c rank_d distribution across datasets with
	 * that derived from a set of datasets with randomly ordered correlation vectors (ie a null distribution). 
	 * A significance p-value is calculated for this gene, and \a -log(p) values are stored in master_rank. 
	 */
	static bool OrderStatisticsRankAggregation(const ushort&, const ushort&,
		ushort**, const vector<ushort> &, vector<float>&, const ushort&);
	static bool OrderStatisticsPreCompute();


	/*!
	 * \brief
	 * Simulates a dataset weight for one-gene query
	 *
	 * \param sQuery
	 * The query
	 * 
	 * \param sDataset
	 * The dataset
	 * 
	 * \param rate
	 * RBP parameter \a p
	 * 
	 * \param percent_required
	 * Percentage of query genes required to be present in a dataset (assumed to be 1 in this case)
	 * 
	 * \param bSquareZ
	 * Whether or not to square \a correlations
	 * 
	 * \param rrank
	 * Final gene-score
	 * 
	 * \param goldStd
	 * Gold-standard gene-set for weighting a dataset
	 * 
	 * This function is mainly used for <em>equal weighting</em>. 
	 * Although <em>equal weighting</em> integrates all datasets with weight = 1, 
	 * for the purpose of displaying datasets, the datasets need to be ranked according to the distance to the
	 * average gene-ranking. 
	 *
	 * This <em>average gene-ranking</em> is produced by summing gene-rankings from all datasets and divided by the number of datasets.
	 * To score a dataset, we calculate the RBP precision of this dataset in retrieving the top 100 genes of the <em>average ranking</em>.
	 * 
	 */
	static bool OneGeneWeighting(CSeekQuery&, 
		CSeekDataset&, const float&, const float&, const bool&, 
		vector<ushort>*, const CSeekQuery*);
};


}
#endif
