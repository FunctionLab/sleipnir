#include "stdafx.h"

/*!
 * \page Funcographer Funcographer
 * 
 * Funcographer mines heavy subgraphs (clusters) from three inputs: gene set functional associations
 * (e.g. from \ref Funcifier), functional activities between datasets and functions (e.g. from
 * \ref BNTruster), and correlations between datasets' functional activities.  This provides a way of
 * mining large data collections for pathways that are interacting and the datasets in which those
 * pathways are most active.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Funcographer -f <functions.dab> -t <trusts.pcl> -d <datasets.dab> -o <cograph.dab>
 *		-r <initial_specificity> -w <final_specificity_ratio> -Z <datasets>
 * \endcode
 * 
 * Using the process/process associations in \c functions.dab, the dataset/process associations in
 * \c trusts.pcl (from \ref BNTruster), and the dataset/dataset similarities in \c datasets.dab, construct
 * a single graph \c cograph.dab containing both node types and all edges.  Then mine the \c functions.dab
 * portion of this graph, outputting (to standard output) all heavy subgraphs with an initial specificity
 * ratio of at least \c initial_specificity, a final specificity ratio that's at least
 * \c final_specificity_ratio fraction of the initial value, and the \c datasets datasets most strongly
 * associated with those functions.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Funcographer/Funcographer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Input DAT/DAB file containing pairwise functional associations between gene sets (e.g. from
 *		\ref Funcifier).</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>None</td>
 *	<td>PCL text file</td>
 *	<td>Input PCL file containing functional activities for each dataset/gene set pair (e.g. from
 *		\ref BNTruster).</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Input DAT/DAB file containing shared functional activity scores between datasets (e.g. from
 *		running \ref Distancer on \ref BNTruster output).</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdou</td>
 *	<td>DAT/DAB file</td>
 *	<td>Output DAT/DAB file containing nodes for both gene sets and datasets; gene set/gene set scores
 *		are taken from \c -f, gene set/dataset scores from \c -t, and dataset/dataset scores from \c -d.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>0</td>
 *	<td>Double</td>
 *	<td>Adjustment value added to each dataset/dataset score; useful for ensuring that all three score types
 *		fall within more or less the same range.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>-1</td>
 *	<td>Integer</td>
 *	<td>Number of dense subgraphs to output; -1 will run until no appropriate subgraphs are available.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>0.5</td>
 *	<td>Double</td>
 *	<td>Ratio of initial to final specificity scores for heavy subgraphs.  For example, if searching for
 *		dense subgraphs with a seed specificity of 25, a ratio of 0.5 will stop building the subgraph
 *		when a specificity of 12.5 is reached.  Value should not exceed \c -r or fall below 1 / \c -r.</td>
 * </tr><tr>
 *	<td>-r</td>
 *	<td>25</td>
 *	<td>Double</td>
 *	<td>Initial specificity score for heavy subgraphs.  Guarantees that the ratio of in- to out-connectivity
 *		for a cluster seed is at least the given value.  A value of 0 will find full cliques of edges above
 *		\c -c instead of dense subgraphs.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>Minimum number of gene sets a heavy subgraph must contain to be output.</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>10</td>
 *	<td>Integer</td>
 *	<td>Number of datasets to be output in association with each heavy subgraph.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
