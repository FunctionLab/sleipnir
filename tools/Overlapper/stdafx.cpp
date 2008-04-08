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
