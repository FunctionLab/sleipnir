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
 * \page KNNImputer KNNImputer
 * 
 * KNNImputer imputs missing values in microarray data as described in Troyanskaya et al 2001.  Given an
 * input PCL, for each gene with missing values, some number of nearest neighbors (by a configurable
 * similarity measure) are found, and the missing value is replaced with a weighted average of the equivalent
 * value in those neighbors.  KNNImputer can optionally remove genes with too many missing values to
 * impute.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * KNNImputer -i <data.pcl> -o <imputed.pcl>
 * \endcode
 * 
 * Replace missing values in the microarray data \c data.pcl based on their nearest neighbors, remove any
 * genes with too many missing values, and save the result in \c imputed.pcl.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include KNNImputer/KNNImputer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>PCL text file</td>
 *	<td>Input PCL file in which missing values are to be imputed.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>PCL text file</td>
 *	<td>Output PCL file in which missing values have been replaced and genes with too many missing values
 *		have been removed.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>10</td>
 *	<td>Integer</td>
 *	<td>Number of neighbors to use for each missing value imputation.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>euclidean</td>
 *	<td>euclidean, pearson, kendalls, kolm-smir, spearman, pearnorm, or hypergeom</td>
 *	<td>Similarity measure to use for finding nearest neighbors.  The default (Euclidean distance) is
 *		highly recommended.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>0.7</td>
 *	<td>Double</td>
 *	<td>Fraction of a gene's expression vector that must be present; genes with less than this many
 *		non-missing values are removed from the output.  For example, in a PCL with 10 columns, genes with
 *		more than three missing values would be removed by default.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Gene text file</td>
 *	<td>If given, only genes in the given gene set are included in the output.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>PCL text file</td>
 *	<td>If given, a PCL file with dimensions equal to the data given with \c -i.  However, the values in the
 *		cells of the weights PCL represent the relative weight given to each gene/experiment pair.  If no
 *		weights file is given, all weights default to 1.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, autocorrelate similarity scores (find the maximum similarity score over all possible
 *		lags of the two vectors; see Sleipnir::CMeasureAutocorrelate).</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>-1</td>
 *	<td>Integer</td>
 *	<td>Maximum number of genes in input file before in-memory score caching is disabled.  If -1, caching is
 *		never performed.  Caching greatly speeds up processing, but can consume large amounts of memory
 *		for inputs with many genes (rows).</td>
 * </tr></table>
 */
