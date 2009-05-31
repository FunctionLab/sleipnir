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
 * \page Matcher Matcher
 * 
 * Matcher calculates the similarity between pairs of input datasets using a variety of similarity/distance
 * measures.  It is similar to \ref MIer, but makes no assumption that the input datasets cover the same
 * gene sets (or even organisms), and it optimized to quickly assess subsets of the input datasets.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Matcher -d <data1_dir> <data2.dab>*
 * \endcode
 * 
 * Compute pairwise similarities between each pair of datasets from the directory of DAT/DAB files
 * \c data1_dir and the individual DAT/DAB files \c data2.dab, outputting a table of scores to standard output.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Matcher/Matcher.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * 	<td>None</td>
 * 	<td>None</td>
 * 	<td>DAT/DAB files</td>
 * 	<td>Datasets for which pairwise similarities will be calculated relative to the datasets in \c -i.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Directory of DAT/DAB files for which pairwise similarities will be calculated relative to the
 *		datasets given on the command line.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>kolm-smir</td>
 *	<td>pearson, quickpear, euclidean, kendalls, kolm-smir, hypergeom, innerprod, bininnerprod, mi</td>
 *	<td>Similarity measure to be used for dataset comparisons.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>Minimum number of data points to subsample from each dataset.</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>1000000000</td>
 *	<td>Integer</td>
 *	<td>Maximum number of data points to subsample from each dataset.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, format output as a two-dimensional table; otherwise, format as a list of pairs.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.</td>
 * </tr></table>
 */
