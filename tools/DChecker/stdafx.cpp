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
 * \page DChecker DChecker
 * 
 * DChecker inputs a gold standard answer file and a DAT/DAB of predicted functional relationships (or
 * other interactions) and outputs the information necessary to perform a performance analysis (ROC curve,
 * precision/recall curve, or AUC score) for the given predictions.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * DChecker -w <answers.dab> -i <predictions.dab>
 * \endcode
 * 
 * Output (to standard output) positive gene counts, true and false positive pair counts, true and false
 * negative pair counts, and an overall AUC score using the default binning of continuous data in
 * \c predictions.dab and the binary gold standard answers in \c answers.dab.
 * 
 * \code
 * DChecker -w <answers.dab> -i <predictions.dab> -f
 * \endcode
 * 
 * Output true/false positive and negative pair counts assuming that \c predictions.dab contains only a
 * finite number of different values and using these as bins.  This is appropriate for inherently
 * discrete data, e.g. cocluster counts.
 * 
 * \code
 * DChecker -w <answers.dab> -i <predictions.dab> -b 0 -n -m 0 -M 1 -e 0.01
 * \endcode
 * 
 * Output true/false positive and negative pair counts by normalizing the given data \c predictions.dab to
 * the range [0,1] and creating bin cutoffs at 0.01 increments between 0 and 1.  This can provide a finer
 * grained binning than the default \c -b setting for some prediction/data sets (and a less fine grained
 * binning for others; when in doubt, try both).
 * 
 * \code
 * DChecker -w <answers.dab> -i <predictions.dab> -c <context.txt>
 * \endcode
 * 
 * Output true/false positive and negative pair counts for the given \c predictions.dab using only the gene
 * pairs relevant to the given biological function \c context.txt.  This is appropriate for evaluating
 * context-specific functional relationship predictions.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include DChecker/DChecker.ggo
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
 * 	<td>If given, contexts in which multiple context-specific evaluations are performed.  Each gene set is
 *		read, treated as a "term" filter (see Sleipnir::CDat::FilterGenes) on the given answer file, and a
 *		context-specific evaluation is saved in the directory \c -d.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>DAT/DAB file</td>
 *	<td>Input DAT, DAB, DAS, or PCL file.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Functional gold standard for learning.  Should consist of gene pairs with scores of 0 (unrelated),
 *		1 (related), or missing (NaN).</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>If multiple contexts are being checked, output directory in which individual contexts' score files
 *		are placed.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>1000</td>
 *	<td>Integer</td>
 *	<td>If nonzero, number of quantile bins into which input scores are sorted.  Each bin is then used as a
 *		cutoff for predicted positives and negatives.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, assume the input predictions contain a small, finite number of distinct values and bin
 *		quantiles appropriate.  Bad things will happen if \c -f is on and there are actually a large number
 *		of distinct input values.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>0</td>
 *	<td>Float</td>
 *	<td>If \c -b is zero and \c -f is off, minimum input score to treat as a positive/negative cutoff.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>1</td>
 *	<td>Float</td>
 *	<td>If \c -b is zero and \c -f is off, maximum input score to treat as a positive/negative cutoff.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>0.01</td>
 *	<td>Double</td>
 *	<td>If \c -b is zero and \c -f is off, size of step to take for cutoffs between \c -m and \c -M.</td>
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
 *	<td>-C</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs passing an "edge" filter against the list.  For details, see
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
 *	<td>-s</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output sum of squared error between input predictions and answer file (assumes a continuous
 *		rather than discrete answer file).</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
