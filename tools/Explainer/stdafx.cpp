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
 * \page Explainer Explainer
 * 
 * Explainer display information on the best- and/or worst-predicted gene pairs in a functional relationship
 * network or, equivalently, the most and least similar gene pairs in two experimental datasets.  This can be
 * used to analyze predictions (to understand where errors are coming from) or to analyze experimental data
 * (to see where and why laboratory results disagree with each other or a gold standard).
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Explainer -i <predictions.dab> -w <answers.dab> [-f SGD_features.tab]
 * \endcode
 * 
 * Output (to standard output) all gene pairs in \c predictions.dab as compared to \c answers.dab, sorted
 * from greatest to least difference; the optional \c SGD_features.tab file can make this output more
 * informative for networks containing yeast genes.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Explainer/Explainer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>DAT/DAB file</td>
 *	<td>Input DAT/DAB file to be compared against a gold standard answer file.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Gold standard answer DAT/DAB file against which input predictions/data are compared.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>-1</td>
 *	<td>Integer</td>
 *	<td>Number of gene pair comparisons to output; -1 displays all pairs.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output only gene pairs marked as positive (1) in the gold standard.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output all gene pairs with any data; if off, only output gene pairs present in the gold
 *		standard.</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>exclude</td>
 *	<td>exclude, include, or only</td>
 *	<td>If exclude, do not output any gene pairs containing an uncharacterized gene.  If include, allow
 *		gene pairs containing uncharacterized genes.  If only, output only gene pairs containing
 *		uncharacterized genes.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Randomly subsample the requested fraction of the possible output data.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which both genes are in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which neither gene is in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs passing a "term" filter against the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output one minus the input's values.</td>
 * </tr><tr>
 *	<td>-r</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, sort output from least to greatest error; otherwise, sort from greatest to least error.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>OBO text file</td>
 *	<td>OBO file containing the structure of the Gene Ontology.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>None</td>
 *	<td>Annotation text file</td>
 *	<td>Gene Ontology annotation file for the desired organism.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>None</td>
 *	<td>SGD features text file</td>
 *	<td>If given, use gene names from the given SGD_features.tab file to label graph nodes.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
