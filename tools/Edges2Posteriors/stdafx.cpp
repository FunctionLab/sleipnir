#include "stdafx.h"

/*!
 * \page Edges2Posteriors Edges2Posteriors
 * 
 * Edges2Posteriors calculates the individual contribution of every dataset to every gene pair in a predicted
 * functional relationship network.  This is similar to \ref BNTruster, but the contributions are broken
 * down per gene pair (rather than per dataset value, and the output can thus be \e much larger) and given
 * as the actual difference in posterior probability incurred by each gene pair's individual data values.
 * 
 * \section sec_overview Overview
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Edges2Posteriors -i <predictions.dab> -n <learned.xdsl> -d <data_dir>
 * \endcode
 * 
 * Output (to standard output) the contribution of each dataset DAT/DAB file in \c data_dir to the
 * predicted functional relationship probabilities in \c predictions.dab based on the Bayesian network
 * \c learned.xdsl; this network must contain node IDs corresponding to the data file names.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Edges2Posteriors/Edges2Posteriors.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>DAT/DAB file</td>
 *	<td>Input DAT/DAB file containing predicted functional relationships.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>Input Bayesian network for which posterior probabilities are to be analyzed.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory from which DAT/DAB data files are read; must have file names corresponding to the input
 *		(X)DSL's node IDs and be accompanied by appropriate QUANT files agreeing with that (X)DSL's
 *		probability tables.</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>None</td>
 *	<td>Tab-delimited text file</td>
 *	<td>If given, argument must be a tab-delimited text file containing two columns, the first node
 *		IDs and the second bin numbers (zero indexed).  For each node ID present in this file, missing values
 *		will be substituted with the given bin number.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>None</td>
 *	<td>Gene pair text file</td>
 *	<td>Tab-delimited text file containing one pair of gene IDs per row.  If given, lookup values for
 *		all requested pairs.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
