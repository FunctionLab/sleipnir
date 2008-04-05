#include "stdafx.h"

/*!
 * \page Contexter Contexter
 * 
 * Contexter calculates gene/context association weights, either from a database and (X)DSL files on an
 * as-needed basis (like \ref BNServer) or from a given pre-calculated global functional relationship
 * network (DAT/DAB file) and context files.
 * 
 * \section sec_overview Overview
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Contexter -i <contexts.txt> -d <database_dir> -c <contexts_genes.txt> -e <context>
 *		-g <genes.txt> -n <contexts_dir> -b <global.xdsl>
 * \endcode
 * 
 * Output (to standard output) the specific association score for each gene in the genome to context number
 * \c context, using (as with \ref BNServer) the context ID/name mappings in \c contexts.txt, the context/gene
 * mappings in \c contexts_genes.txt, the gene ID/name mappings in \c genes.txt, the Sleipnir::CDatabase
 * in \c database_dir, the context-specific Bayesian classifiers in \c contexts_dir, and the global
 * Bayesian classifier \c global.xdsl.
 * 
 * \code
 * Contexter -t -i <network.dab> -g <genes.txt> -c <contexts_genes.txt>
 * \endcode
 * 
 * Output the specific association score for each gene in the genome to context number \c context, using the
 * context/gene mappings in \c contexts_genes.txt, the gene ID/name mappings in \c genes.txt, and the
 * pre-computed interaction network \c network.dab.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Contexter/Contexter.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>DAT/DAB file or text file</td>
 *	<td>When \c -t is on, an input DAT/DAB file containing a global predicted functional relationship
 *		network.  When \c -t is off, a tab-delimited text file from which context indices, names, and
 *		string IDs are read.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, read a global functional relationship network from \c -i; if off, read a context ID
 *		file (along with additional context/gene information from other arguments).</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>When \c -t is off, directory from which Sleipnir::CDatabase data files are read.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>Tab-delimited text file containing two columns, context IDs and gene IDs (both one-based).</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>None</td>
 *	<td>Integer</td>
 *	<td>Context ID (one-based) for which gene associations should be calculated.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>When \c -t is on, tab-delimited text file containing two columns, gene IDs (one-based) and
 *		unique gene names (matching those in the DAT/DAB file).</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory from which context-specific Bayesian classifiers ((X)DSL files) are read.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>Bayesian classifier for the default (global) context.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, assume XDSL files will be used instead of DSL files.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, read binary stored Bayesian classifiers from a file specified by \c -n; if off, assume
 *		\c -n specifies a directory of (X)DSL files.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>None</td>
 *	<td>Binary file</td>
 *	<td>If given, store Bayesian classifiers in a custom binary format in the given filename.  Cannot be
 *		used with \c -m on.  Loading the classifiers from a binary file can be faster than loading several
 *		hundred separate (X)DSLs.  If your classifiers aren't changing, you can load them from (X)DSL
 *		files once, leaving \c -m off and saving them with \c -M, then on subsequent runs turn \c -M off and
 *		load the binary file with \c -m on.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>None</td>
 *	<td>String</td>
 *	<td>If given, calculate only the association weight for the given gene name rather than for the entire
 *		genome.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
