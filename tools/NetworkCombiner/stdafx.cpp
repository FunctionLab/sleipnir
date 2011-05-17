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
 * \page NetworkCombiner NetworkCombiner
 * 
 * NetworkCombiner can replace the DAT mode of \ref Combiner.  It is faster, but achieves this by making some
 * important assumptions.  NetworkCombiner uses the mean method of combining and assumes that all gene pairs 
 * be present in all DABs being combined. Therfore, NetworkCombiner should only be used if these assumptions 
 * are met.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * Create a new DAB file from a directory of existing DAB files by calculating the mean for each gene pair
 * across all of the existing DAB files.
 * 
 * \code
 * NetworkCombiner -d <directory of dabs> -o <combined.dab>
 * \endcode
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include NetworkCombiner/NetworkCombiner.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * 	<td>-d</td>
 * 	<td>None</td>
 * 	<td>DAT/DAB Directory</td>
 * 	<td>Input directory (must only contain input files as DAT/DAB).</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>Map gene index among the network dabs to combine (should be used when the gene indices are not identical among network dabs).</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>DAB file</td>
 *	<td>Output file for combined network.</td>
 * </tr></table>
 */
