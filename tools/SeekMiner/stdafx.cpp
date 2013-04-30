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

/*!
 * \page SeekMiner SeekMiner
 * 
 * SeekMiner returns a gene-ranking based on the coexpressions to the user-specified
 * query genes. It finds relevant datasets by using one of the many dataset weighting
 * algorithms, including the query-coexpression weighting, the order statistics 
 * weighting, etc. Afterward, it performs a weighted integration of coexpressions
 * using the computed dataset weights.
 * The search algorithms employed by Seek are designed to be quick and efficient, and
 * they support the on-the-fly weight calculations for thousands of microarray 
 * datasets.
 *
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * SeekMiner -x <dset_platform_map> -D <search_dset> -i <gene_map> -q <query> -P <platform_dir> -p <prep_dir> -n <num_db> -d <db_dir> -Q <quant> -o <output_dir> -m
 * \endcode
 * This performs the coexpression mining for a list of queries in the file \c query, and outputs
 * the gene-ranking and the dataset weights in the \c output_dir.
 * The outputs are divided by queries. Starting with the first query (named 0), its final results
 * will consist of three files: \c 0.query, 0.dweight, 0.gscore.
 * \li The file base name (0) indicates the query index in the list.
 * \li The \c 0.query stores the space-delimited query gene-set in text.
 * \li The \c 0.dweight stores the weightings of datasets as a binary one-dimensional float vector
 * (see SeekEvaluator for displaying a DWEIGHT extension file).
 * \li The \c 0.gscore stores the gene scores as a binary one-dimensional float vector
 * (see SeekEvaluator for displaying a GSCORE extension file).
 *
 * Users can choose among several dataset weighting schemes:
 * \li Query cross-validated weighting (shortened CV, default), where we iteratively use a subset of the
 * query to construct a search instance to retrieve the remaining query genes. The retrieval precision
 * score becomes the dataset weight.
 * \li Equal weighting (shortened EQ), where all datasets are weighted equally.
 * \li Order statistics integration (shortened ORDER_STAT), which is outlined in Adler et al (2009).
 * This method computes P-value statistics by comparing the rank of correlation across datasets to the
 * ranks that would have been generated a null distribution (where correlations are randomly scattered
 * and all ranks appear equally likely).
 *
 * Users can also choose the type of transformations to perform on the distance matrix of <b>each dataset</b>.
 * Among these are:
 * \li No transformation: use only Pearson correlations. (\c -T)
 * \li Default z-scores: Fisher's transform followed by standardization. No arguments required.
 * \li On top of the default z-scores, subtracts the average z-score of the result genes. (\c -m).
 * It is highly recommended to turn \c -m on.
 *
 * Users may also define the datasets that they wish to use for integrations in a query-specific way, using \c -D  argument.
 * If this argument is absent, all datasets in the compendium will be integrated. If \c -D is used, the structure of \c search_dset should be:
 * \code
 * GSE15913.GPL570.pcl GSE16122.GPL2005.pcl GSE16836.GPL570.pcl ...
 * GSE14933.GPL570.pcl GSE15162.GPL2005.pcl GSE15566.GPL570.pcl ...
 * ...
 * \encode
 * where each line, corresponding to a query, is a space-separated dataset list.
 *
 * The remainder of arguments are query-independent search setting files and environment files.
 * 
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include SeekMiner/SeekMiner.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>The dataset-platform mapping file. 
 * Tab-delimited text file containing two columns, the dataset name,
 * and the corresponding platform name. For example, below is a few sample lines
 * of what datasets and platforms can be:
 * \code
 * GSE15913.GPL570.pcl  GPL570
 * GSE16122.GPL2005.pcl GPL2005
 * GSE16797.GPL570.pcl  GPL570
 * GSE16836.GPL570.pcl  GPL570
 * GSE17351.GPL570.pcl  GPL570
 * GSE17537.GPL570.pcl  GPL570
 * \endcode
 * Note that although the dataset name looks like a file name, it does not 
 * need to be a valid file name, as long as it properly and uniquely describes
 * the dataset. Here, the dataset is uniquely identified by a GSE ID and a GPL ID 
 * combination. In addition, the ordering of the datasets in this file is 
 * non-random: it must match the order of the datasets in the CDatabaselet 
 * (\c *.db).
 * Typically, the \c dataset_map.txt file that generated the
 * CDatabaselet collection can be used for this field.
 *	</td>
 * </tr><tr>
 *	<td>-D</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>The search dataset file. This multi-line file contains the list of 
 * space-separated dataset names to be used for the coexpression search. 
 * Each query can use its own list, so each query requires its own definition in 
 * a separate line. Below are some sample lines:
 * \code
 * GSE12386.GPL7149.pcl GSE20934.GPL10200.pcl GSE17941.GPL9188.pcl GSE13901.GPL7763.pcl GSE10616.GPL5760.pcl ...
 * GSE12386.GPL7149.pcl GSE20934.GPL10200.pcl GSE17941.GPL9188.pcl GSE13901.GPL7763.pcl ...
 * \endcode
 * This specifies the search datasets for two queries. The dataset 
 * must be one of the dataset names from the dataset-platform map file (-x).
 * </td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>The gene-map file. Maps the genes to an ID between 0 to N where N is
 * the genome size. Example:
 * \code
 * 1    1
 * 2    10
 * 3    100
 * 4    1000
 * 5    10000
 * 6    100008589
 * 7    100009676
 * 8    10001
 * 9    10002
 * 10   10003
 * 11   100033413
 * 12   100033414
 * \endcode
 * The ordering of the genes in this file must match the order of genes 
 * in the CDatabaselets (\c *.db). Typically, the \c gene_map.txt file 
 * that generated the CDatabaselet collection can be used for this field.
 * </td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>.</td>
 *	<td>Text file</td>
 *	<td>Input file containing list of CDatabaselets to combine</td>
 * </tr></table>
 */
