#include "stdafx.h"

/*!
 * \page Overlapper Overlapper
 * 
 * Overlapper outputs a confusion matrix for two discretized input DAT/DAB files, usually gold standards.
 * This summarizes the degree to which the two inputs overlap and agree (or disagree) for all pairwise
 * scores.
 * 
 * \section sec_overview Overview
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
