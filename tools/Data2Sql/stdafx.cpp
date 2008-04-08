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
 * \page Data2Sql Data2Sql
 * 
 * Data2Sql converts a collection of DAT/DAB files into relational tables appropriate for insertion into a
 * SQL database.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Data2Sql -i <genes.txt> -t <table> <data.dab>*
 * \endcode
 * 
 * Output (to standard output) SQL commands to construct a table named \c table containing the pairwise data
 * in the DAT/DAB files \c data.dab, with gene names mapped to numerical IDs using \c genes.txt
 * (tab-delimited file with two columns, one-based integer indices and gene names).
 * 
 * \code
 * Data2Sql -i <genes.txt> -d <data.dab>*
 * \endcode
 * 
 * Output (to standard output) a table containing a numerical index for each input dataset \c data.dab.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Data2Sql/Data2Sql.ggo
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
 * 	<td>Input DAT/DAB files from which data is drawn to be converted into the output SQL file.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Text file</td>
 *	<td>Tab-delimited text file containing two columns, numerical gene IDs (one-based) and unique gene names
 *		(matching those in the input DAT/DAB files).</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>datapairs</td>
 *	<td>String</td>
 *	<td>Database table name.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output a table listing dataset ID/name relations; if off, output a table listing individual
 *		gene pair values from the input datasets.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>1000</td>
 *	<td>Integer</td>
 *	<td>Initiate a new INSERT command after each block of this many data values.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
