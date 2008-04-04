#include "stdafx.h"

/*!
 * \page Hubber Hubber
 * 
 * Hubber calculates the cohesiveness of one or more input gene sets in a given functional relationship (or
 * experimental data) network; it also calculates the association of each gene with those gene sets.  This
 * is primarily useful for predicting function for genes based on "guilt by association" with genes already
 * known to be in those functions (i.e. gene sets).
 * 
 * \section sec_overview Overview
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Hubber/Hubber.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * 	<td>None</td>
 * 	<td>None</td>
 * 	<td>Gene text files</td>
 * 	<td>Gene sets which will be tested for cohesiveness and gene associations in the given dataset.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>DAT/DAB file</td>
 *	<td>Interaction network which will be analyzed for gene/function association.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>Number of genes to output in association with each gene set.  0 will display no genes, only gene
 *		set cohesiveness scores, and -1 will produce a table of every gene's association with every
 *		function.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>If given, for each gene set, produce a DAB file in the given directory containing only the genes
 *		(and edges) in that gene set.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>Binary file</td>
 *	<td>If given, binary file into which each gene's average connection weight to every other gene in
 *		each context is placed.  Used with \ref BNServer.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>If \c -o is given, tab-delimited text file containing one header row and at least two columns, the
 *		second of which must be a context name corresponding to that context's functional relationship DAT/DAB
 *		file.  Format nearly identical to the \ref BNServer context list.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>If \c -o is given, tab-delimited text file containing one header row and at least two columns, the
 *		second of which must be a gene name; should list all genes for which backgrounds are to be calculated.
 *		Format nearly identical to the \ref BNServer gene list.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>If \c -o is given, directory from which DAT/DAB files with names corresponding to contexts from \c -x
 *		are loaded.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
