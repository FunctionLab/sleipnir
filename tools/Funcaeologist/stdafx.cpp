#include "stdafx.h"

/*!
 * \page Funcaeologist Funcaeologist
 * 
 * Funcaeologist evaluates the cohesiveness of a given gene set within a collection of functional
 * relationship networks or experimental result datasets.  The cohesiveness is a normalized ratio of the
 * average in-connectivity (weight of an edge connecting two genes in the set) to out-connectivity (weight
 * of an edge including at least one gene in the set).  Gene sets with higher cohesiveness a overall more
 * functionally related and generally behave more similarly in experimental data.
 * 
 * \section sec_overview Overview
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Funcaeologist/Funcaeologist.ggo
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
 * 	<td>Functional relationship networks or datasets in which the cohesiveness of the given gene set is
 *		evaluated.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>stdin</td>
 *	<td>Gene text file</td>
 *	<td>Gene set (one gene ID per line) to be investigated for cohesiveness in the given networks.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory containing functional relationship or dataset DAT/DAB files.  DAT/DAB files in this
 *		directory must have names corresponding to the contexts given on the command line (e.g. if the
 *		context file "mitotic_cell_cycle.txt" is given on the command line, "mitotic_cell_cycle.dab" must
 *		exist in the directory).</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
