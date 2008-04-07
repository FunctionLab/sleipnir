#include "stdafx.h"

/*!
 * \page Overlapper Overlapper
 * 
 * Overlapper outputs a confusion matrix for two discretized input DAT/DAB files, usually gold standards.
 * This summarizes the degree to which the two inputs overlap and agree (or disagree) for all pairwise
 * scores.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Overlapper -i <data1.dab> -I <data2.dab>
 * \endcode
 * 
 * Outputs a confusion matrix comparing values for gene pairs in \c data1.dab and \c data2.dab, both of
 * which must have associated QUANT files.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Overlapper/Overlapper.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>First input DAT/DAB file to inspect.</td>
 * </tr><tr>
 *	<td>-I</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Second input DAT/DAB file to inspect.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
