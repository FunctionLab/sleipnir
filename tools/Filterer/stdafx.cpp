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
 * \page Filterer Filterer
 * 
 * Filterer removes values within a given range or ranges from a DAT/DAB file.  These ranges can be
 * specified using a rudimentary query language allowing inclusions, exclusions, intersections,
 * and unions.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Filterer -i <input.dab> -o <output.dab> i0=1
 * \endcode
 * 
 * Include in \c output.dab only the gene pair values from \c input.dab that are between 0 and 1.
 * 
 * \code
 * Filterer -i <input.dab> -o <output.dab> x0=1
 * \endcode
 * 
 * Include in \c output.dab only the gene pair values from \c input.dab that are not between 0 and 1.
 * 
 * \code
 * Filterer -i <input.dab> -o <output.dab> i0=1 i3=
 * \endcode
 * 
 * Include in \c output.dab only the gene pair values from \c input.dab that are between 0 and 1 or
 * at least 3.
 * 
 * \code
 * Filterer -i <input.dab> -o <output.dab> x0=1 i=-0.5
 * \endcode
 * 
 * Include in \c output.dab only the gene pair values from \c input.dab that are at least -0.5 and
 * not between 0 and 1.
 * 
 * \code
 * Filterer -i <input.dab> -o <output.dab> i=-2 i-1= x-2.5=-0.5
 * \endcode
 * 
 * Include in \c output.dab only the gene pair values from \c input.dab that are at most -2 or at
 * least -1 and are not between -2.5 and -0.5.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Dat2Dab/Dat2Dab.ggo
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
 *	<td>Input DAT/DAB file.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>DAT/DAB file</td>
 *	<td>Output (filtered) DAT/DAB file.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
