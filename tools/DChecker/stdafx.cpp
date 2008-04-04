#include "stdafx.h"

/*!
 * \page DChecker DChecker
 * 
 * DChecker inputs a gold standard answer file and a DAT/DAB of predicted functional relationships (or
 * other interactions) and outputs the information necessary to perform a performance analysis (ROC curve,
 * precision/recall curve, or AUC score) for the given predictions.
 * 
 * \section sec_overview Overview
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
