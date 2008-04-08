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
 * \page Randomizer Randomizer
 * 
 * Randomizer removes all values from the given DAT/DAB and randomly inserts the same number (or a requested
 * number) of ones.  This can be used to quickly create a random gold standard with a known number of positive
 * examples.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Randomizer -i <data.dab>
 * \endcode
 * 
 * Removes all non-missing values in \c data.dab and assigns an equivalent number of 1s randomly across all
 * gene pairs.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Randomizer/Randomizer.ggo
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
 *	<td>Input DAT/DAB file to be randomized.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>If nonzero, desired number of non-missing values in the output file.</td>
 * </tr></table>
 */
