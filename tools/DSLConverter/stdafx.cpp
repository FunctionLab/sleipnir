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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"

/*!
 * \page DSLConverter DSLConverter
 * 
 * DSLConverter converts between XDSL and DSL files and vice versa.  This is useful in that certain versions
 * of SMILE won't parse XDSL files on Windows, and other versions won't parse DSL files on Linux.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * DSLConverter -i <input.xdsl> -o <output.dsl>
 * \endcode
 * 
 * Output a DSL file \c output.dsl equivalent to the input XDSL file \c input.xdsl.  Of course, you can
 * input a DSL file and output an XDSL file instead.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include DSLConverter/DSLConverter.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>Input (X)DSL file.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>Output (X)DSL file.</td>
 * </tr></table>
 */
