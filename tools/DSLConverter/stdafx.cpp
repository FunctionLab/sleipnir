#include "stdafx.h"

/*!
 * \page DSLConverter DSLConverter
 * 
 * DSLConverter converts between XDSL and DSL files and vice versa.  This is useful in that certain versions
 * of SMILE won't parse XDSL files on Windows, and other versions won't parse DSL files on Linux.
 * 
 * \section sec_overview Overview
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
