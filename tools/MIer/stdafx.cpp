#include "stdafx.h"

/*!
 * \page MIer MIer
 * 
 * MIer calculates the mutual information (or other similarity measure) between pairs of input datasets.
 * This can be used to approximate how much information is shared between two experimental datasets or
 * how similar two predicted functional relationship networks are.
 * 
 * \section sec_overview Overview
 * 
 * \subsection ssec_usage_detailed Detailed Usage
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
