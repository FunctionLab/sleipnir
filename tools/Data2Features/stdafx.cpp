#include "stdafx.h"

/*!
 * \page Data2Features Data2Features
 * 
 * Data2Features converts a collection of DAT/DAB or PCL files into a features file appropriate for use with
 * the Weka machine learning package.  Data2Features focuses on machine learning about entire datasets,
 * collapsing a collection of per-gene or per-gene-pair values into a single summary feature.  Similar to
 * \ref Data2Bnt.
 * 
 * \section sec_overview Overview
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Data2Features -p <positives.txt> -e <features.txt> -d <data.txt> <data.pcl/dab>*
 * \endcode
 * 
 * Output (to standard output) an XRFF feature file appropriate for use with Weka; each example represents
 * a dataset from the microarray PCL files \c data.pcl or DAT/DAB files \c data.dab, with features as
 * specified in \c features.txt, defaults from \c data.txt, and values calculated as the average across
 * gene pairs in \c positives.txt.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Data2Features/Data2Features.ggo
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
 * 	<td>Input DAT/DAB files from which data is drawn for features in the output Weka file.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>stdin</td>
 *	<td>Gene text file</td>
 *	<td>List of genes labeled as positives for machine learning; can be drawn from the same
 *		pathway/process/complex/GO term/etc.</td>
 * </tr><tr>
 *	<td>-e</td>
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
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>SGD features text file</td>
 *	<td>SGD_features.tab file; if given, process only genes appearing in this file.</td>
 * </tr><tr>
 *	<td>-D</td>
 *	<td>pearnorm</td>
 *	<td>pearson, euclidean, kendalls, kolm-smir, spearman, pearnorm, hypergeom, innerprod, bininnerprod,
 *		quickpear, or mi</td>
 *	<td>Similarity measure to be used for converting PCL inputs into pairwise scores.</td>
 * </tr><tr>
 *	<td>-N</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to z-scores (subtract mean, divide by standard deviation) before
 *		processing.</td>
 * </tr><tr>
 *	<td>-S</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
