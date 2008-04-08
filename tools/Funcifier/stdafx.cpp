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
 * \page Funcifier Funcifier
 * 
 * Funcifier (pronounced funk-if-eyer, not fun-siffier, with flagrant disregard for rules of standard
 * English) calculates the strength of functional association between processes or pathways in a given
 * functional relationship (or data) network.  A functional association between two gene sets (e.g. pathways)
 * measures the amount of cross-talk between those processes, which is calculated from the given individual
 * gene pair functional relationships.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Funcifier -i <network.dab> -o <functions.dab> -l <colors.txt> <contexts.txt>*
 * \endcode
 * 
 * Find the strength of functional association between every pair of contexts \c contexts.txt in the
 * functional relationship network \c network.dab and store the resulting function pair network in
 * \c functions.dab; if requested, output the cohesiveness score for each function in \c colors.txt for
 * use with \ref Dat2Graph.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Funcifier/Funcifier.ggo
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
 * 	<td>Gene sets which will be tested for functional association in the given dataset.  Each set
 *		generally represents a biological context, e.g. process or pathway.  Elements in the output
 *		DAT/DAB file will be named according to these input filenames.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>DAT/DAB file</td>
 *	<td>Input DAT/DAB file containing the functional relationship or data network to be mined for
 *		gene set functional associations.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>DAT/DAB file</td>
 *	<td>Output DAT/DAB file containing functional associations between gene sets.  While the input
 *		network contains individual gene pairs, the output network contains pairs of gene sets (e.g. each
 *		element is a pair of labels from gene sets given on the command line).</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>discard</td>
 *	<td>discard, ignore, or oneonly</td>
 *	<td>Way in which genes present in multiple gene sets are handled.  \c discard prevents shared genes
 *		from contributing to any gene set's associations, \c ignore allows them to contribute to every
 *		gene set in which they're included, and \c oneonly allows genes to contribute only if they are
 *		shared between exactly two processes (i.e. they have only one overlap).</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>If given, output file into which cohesiveness scores for each gene set are placed.  Output is
 *		compatible with color input for \ref Dat2Graph.</td>
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
