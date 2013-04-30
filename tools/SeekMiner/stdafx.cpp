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
 * SeekMiner ranks the genome based on the coexpressions to the user-specified
 * query genes. It finds relevant datasets using one of the many dataset weighting
 * algorithms, including the query-coexpression weighting, the order statistics 
 * weighting, etc. These algorithms are designed to be quick and efficient, and 
 * they support the on-the-fly weight calculations for thousands of microarray 
 * datasets. Finally, SeekMiner integrates the individual dataset's query-relevant
 * genes to return a master gene-ranking.
 *
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * SeekMiner -x <dset_platform_map> -D <search_dset> -i <gene_map> -q <query> -P <platform_dir> -p <prep_dir> -n <num_db> -d <db_dir> -Q <quant> -b <num_db_in_memory> [-m] [-M] [-r] -o <output_dir> 
 * \endcode
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
