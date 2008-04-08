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
 * \page Hubber Hubber
 * 
 * Hubber calculates the cohesiveness of one or more input gene sets in a given functional relationship (or
 * experimental data) network; it also calculates the association of each gene with those gene sets.  This
 * is primarily useful for predicting function for genes based on "guilt by association" with genes already
 * known to be in those functions (i.e. gene sets).
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Hubber -i <network.dab> -g <genes> <contexts.txt>*
 * \endcode
 * 
 * Output (to standard output) the average in- and out-connectivity of each gene set \c contexts.txt in the
 * network \c network.dab; for each context, display the \c genes genes most specifically associated with
 * that context.  In-connectivity is the weight of an edge between two genes in the context; out-connectivity
 * is the weight of an edge incident to a gene in the context.
 * 
 * \code
 * Hubber -i <network.dab> -g -1 <contexts.txt>*
 * \endcode
 * 
 * Output (to standard output) the specificity with which every gene in the genome is associated with
 * each context \c contexts.txt in the network \c network.dab.  Scores for genes already annotated in the
 * context are negative.
 * 
 * \code
 * Hubber -i <global.dab> -o <backgrounds.bin> -x <context_ids.txt> -s <gene_ids.txt>
 *		-d <contexts_dir> <contexts.txt>*
 * \endcode
 * 
 * Save in \c backgrounds.bin the average weight of edges incident to each gene in the global network
 * \c global.dab and in every context-specific network in \c contexts_dir, each corresponding to a
 * context gene list \c contexts.txt.  \c context_ids.txt and \c gene_ids.txt are tab-delimited text
 * files as described below (and used in \ref BNServer) with one-based indices in the first column and
 * context or gene names, respectively, in the second column.
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
 *	<td>If \c -o is given, tab-delimited text file containing at least two columns, the second of which must
 *		be a context name corresponding to that context's functional relationship DAT/DAB file.  Identical to
 *		the \ref BNServer context list.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>If \c -o is given, tab-delimited text file containing at least two columns, the second of which must
 *		be a gene name; should list all genes for which backgrounds are to be calculated.  Identical to \c -x
 *		in format.</td>
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
