#include "stdafx.h"

/*!
 * \page Dab2Dad Dab2Dad
 * 
 * Dab2Dad combines multiple DAT/DAB files, possibly including a gold standard answer file in addition to
 * multiple datasets, into a DAD file as described by Sleipnir::CDatasetCompact.
 * 
 * \section sec_overview Overview
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Dab2Dad/Dab2Dad.ggo
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
 * 	<td>Input DAT/DAB files to be combined into a DAD.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>DAD file</td>
 *	<td>If given, load this binary DAD file and save it as text to standard output or to \c -o.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>None</td>
 *	<td>DAD file</td>
 *	<td>If given, load this binary DAD file and retain a shared memory map (can theoretically be used to
 *		force a DAD into virtual memory for rapid access by other processes).</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>If given, combine DAT/DABs from \c -w and \c -d into a DAD file appropriate for use with the
 *		given Bayes net.  Discretization and node order information will be read from the given network.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>DAD file</td>
 *	<td>If given, output DAD is written as binary to the requested file.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Must be given with \c -n; indicates an answer file DAT/DAB to be included in the DAD along with
 *		datasets from the command line.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Must be given with \c -n; indicates the directory in which DAT/DAB files corresponding to the
 *		input network's node IDs will be found.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output DAD will include data for all gene pairs; if off, it will only contain data for
 *		gene pairs with a value present in the answer file.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which both genes are in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which neither gene is in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>None</td>
 *	<td>String</td>
 *	<td>If given, lookup all values for pairs involving the requested gene.</td>
 * </tr><tr>
 *	<td>-L</td>
 *	<td>None</td>
 *	<td>String</td>
 *	<td>If given with \c -l, lookup all values for the requested gene pair.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>None</td>
 *	<td>Gene text file</td>
 *	<td>If given with \c -l, lookup all pairs between \c -l and the given gene set.  If given alone,
 *		lookup all pairs between genes in the given set.</td>
 * </tr><tr>
 *	<td>-T</td>
 *	<td>None</td>
 *	<td>Gene pair text file</td>
 *	<td>Tab-delimited text file containing one pair of gene IDs per row.  If given, lookup values for
 *		all requested pairs.</td>
 * </tr><tr>
 *	<td>-q</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, quantize data from input DAT/DABs before outputting lookup results.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>If given, process only gene pairs present in the given DAT/DAB with a score greater than zero.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
