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
 * \page Mat2Txt Mat2Txt
 * 
 * Mat2Txt converts a binary full matrix (Sleipnir::CFullMatrix) into a human-readable text format.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Mat2Txt -i <matrix.bin> -o <matrix.txt>
 * \endcode
 * 
 * Output to \c matrix.txt a tab-delimited text representation of the full floating point matrix stored in
 * binary format in \c matrix.bin.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Mat2Txt/Mat2Txt.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Binary matrix file</td>
 *	<td>Binary full matrix (Sleipnir::CFullMatrix) stored by another tool (e.g. \ref Hubber).</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>Text matrix file</td>
 *	<td>Tab-delimited textual representation of the input matrix.</td>
 * </tr></table>
 */
