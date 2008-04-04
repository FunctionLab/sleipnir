#include "stdafx.h"

/*!
 * \page Data2Sql Data2Sql
 * 
 * Data2Sql converts a collection of DAT/DAB files into relational tables appropriate for insertion into a
 * SQL database.
 * 
 * \section sec_overview Overview
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
