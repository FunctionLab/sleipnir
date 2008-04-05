#include "stdafx.h"

/*!
 * \page Data2Bnt Data2Bnt
 * 
 * Data2Bnt converts a collection of DAT/DAB files into a features file appropriate for use with either
 * the MATLAB Bayes Net Toolkit (BNT) or the Weka machine learning package.  Data2Bnt provides one way to
 * bridge the single gene function/gene pair functional relationship gap by producing one example per
 * gene pair, but labeling these examples based on the inclusion of one of that pair's genes in a positive
 * gene set (e.g. pathway/process/complex/etc.)  Similar to \ref Data2Features.
 * 
 * \section sec_overview Overview
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Data2Bnt -i <positives.txt> -f <features.txt> -d <data.txt> -q <bins.quant> <data.pcl/dab>*
 * \endcode
 * 
 * Generate a feature matrix compatible with BNT in which each row represents a gene pair, labeled positive
 * if one of the genes is in \c positives.txt, with features as described in \c features.txt, default
 * values drawn from \c data.txt, and feature values from the microarray PCL files \c data.pcl or
 * DAT/DAB files \c data.dab.  If PCL files are used, they are discretized using the bin edges in the
 * QUANT file \c bins.quant.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Data2Bnt/Data2Bnt.ggo
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
 * 	<td>Input DAT/DAB files from which data is drawn for features in the output BNT/Weka file.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Gene text file</td>
 *	<td>Set of genes such that pairs including them will be labeled as positive examples.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>Tab-delimited text file containing three columns: feature name, |-delimited feature values, and an
 *		optional default value.  Lines starting with # are ignored as comments.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>Tab-delimited text file containing one dataset per line.  The first tab-delimited token of each
 *		line should be a dataset name, with all subsequent tokens of the form
 *		&lt;feature name>|&lt;feature value>.</td>
 * </tr><tr>
 *	<td>-q</td>
 *	<td>None</td>
 *	<td>QUANT text file</td>
 *	<td>Tab-delimited QUANT file containing exactly one line of bin edges; these are used to discretize
 *		all continuous input data.  For details, see Sleipnir::CDataPair.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>SGD features text file</td>
 *	<td>SGD_features.tab file; if given, process only genes appearing in this file.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Randomly subsample the requested fraction of the available data.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, display only non-default values in the output matrix.  This is incompatible with BNT/Weka and
 *		should only be used for informational purposes.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, include #-prefixed comments in the output indicating which examples represent which gene
 *		pairs.  This is incompatible with BNT/Weka and should only be used for informational purposes.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, generate Weka-compatible XRFF output; if off, generate BNT-compatible textual output.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, generated weighted XRFF output by combining and upweighting identical examples.</td>
 * </tr></table>
 */
