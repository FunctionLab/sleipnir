#include "stdafx.h"

/*!
 * \page Cliquer Cliquer
 * 
 * Cliquer mines an interaction network for dense subgraphs, i.e. clusters.  These can be (by default)
 * heavily weighted subgraphs (clusters) based on a greedy algorithm due to Charikar 2000 or full
 * cliques above some cutoff (computationally expensive!)
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Cliquer -i <network.dab> -r <initial_specificity> -w <final_specificity_ratio>
 * \endcode
 * 
 * Output (to standard output) all heavy subgraphs (clusters) in \c network.dab with initial specificity
 * ratio at least \c initial_specificity and final specificity ratio at least \c final_specificity_ratio
 * fraction of the initial value.
 * 
 * \code
 * Cliquer -i <network.dab> -r 0 -s <subgraphs> -S <size>
 * \endcode
 * 
 * Output at most \c subgraphs cliques of non-missing, non-zero edges in \c network.dab of size \c size.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Cliquer/Cliquer.ggo
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
 *	<td>Interaction network which will be mined for dense subgraphs (clusters/cliques).</td>
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
 *	<td>-s</td>
 *	<td>100</td>
 *	<td>Integer</td>
 *	<td>Number of cliques (not dense subgraphs) to output.</td>
 * </tr><tr>
 *	<td>-S</td>
 *	<td>3</td>
 *	<td>Integer</td>
 *	<td>Size of cliques (not dense subgraphs) to output.  Use with caution; exponential growth is not your
 *		friend.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>If given, ignore all edges present in the given DAT/DAB file during clique finding.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>None</td>
 *	<td>Double</td>
 *	<td>If given, remove all input edges below the given cutoff (after optional normalization).</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
