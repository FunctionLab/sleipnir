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
 * \page MIer MIer
 * 
 * MIer calculates the mutual information (or other similarity measure) between pairs of input datasets.
 * This can be used to approximate how much information is shared between two experimental datasets or
 * how similar two predicted functional relationship networks are.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * MIer <data.dab>*
 * \endcode
 * 
 * Compute pairwise mutual information scores for each pair of datasets in \c data.dab and output them to
 * standard output.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * Note that to use \ref Counter with the output from \c MIer, you should first convert the table of raw
 * (bit) mutual information scores to exponentially scaled sums of relative shared information.  This can
 * be done using the \c half2relative.rb and \c half2weights.rb scripts included with Sleipnir.  The
 * combination of these two files' outputs creates a weights file appropriate for use with \ref Counter 's
 * alphas parameters.
 * 
 * \include MIer/MIer.ggo
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
 * 	<td>Datasets for which pairwise mutual information (or other similarity measure) will be calculated.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>mi</td>
 *	<td>mi, pearson, quickpear, euclidean, kendalls, kolm-smir, hypergeom, innerprod, bininnerprod</td>
 *	<td>Similarity measure to be used for dataset comparisons.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, assume that all missing gene pairs in all datasets have a value of 0 (i.e. the first bin).</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>None</td>
 *	<td>Tab-delimited text file</td>
 *	<td>If given, argument must be a tab-delimited text file containing two columns, the first node
 *		IDs (see \ref BNCreator) and the second bin numbers (zero indexed).  For each node ID present in
 *		this file, missing values will be substituted with the given bin number.</td>
 * </tr><tr>
 *	<td>-R</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, assign missing values randomly; this generally results in much better approximations of
 *		mutual information.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, format output as a tab-delimited table; otherwise, format as one pair per line.</td>
 * </tr><tr>
 *	<td>-y</td>
 *	<td>-1</td>
 *	<td>Integer</td>
 *	<td>If nonnegative, process only pairs of datasets containing (and beginning with) the given dataset
 *		index.  This can be used to parallelize many mutual information calculations by running
 *		processes with different \c -y values.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
