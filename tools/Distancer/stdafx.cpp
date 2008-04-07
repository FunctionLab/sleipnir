#include "stdafx.h"

/*!
 * \page Distancer Distancer
 * 
 * Distancer converts a microarray PCL file into a collection of pairwise similarity scores (DAT/DAB file)
 * using one of a large variety of similarity/distance measures (Pearson correlation, Euclidean distance,
 * etc.)
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Distancer -i <data.pcl> -o <data.dab>
 * \endcode
 * 
 * Output a DAT/DAB file \c data.dab containing pairwise similarity scores calculated from the microarray
 * data in \c data.pcl.  The default settings will produce z-score, z-transformed Pearson correlations;
 * this behavior can be extensively configured using the detailed options.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Distancer/Distancer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>PCL text file</td>
 *	<td>Input PCL file of microarray data to be scored.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>DAT/DAB file</td>
 *	<td>Output DAT/DAB file to contain pairwise scores.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>pearnorm</td>
 *	<td>Lots!</td>
 *	<td>Similarity measure to be used for scoring.  See Sleipnir::IMeasure.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>PCL text file</td>
 *	<td>If given, a PCL file with dimensions equal to the data given with \c -i.  However, the values in the
 *		cells of the weights PCL represent the relative weight given to each gene/experiment pair.  If no
 *		weights file is given, all weights default to 1.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, autocorrelate similarity scores (find the maximum similarity score over all possible
 *		lags of the two vectors; see Sleipnir::CMeasureAutocorrelate).</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to z-scores (subtract mean, divide by standard deviation) before
 *		processing.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output one minus the input's values.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which both genes are in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>None</td>
 *	<td>Double</td>
 *	<td>If given, remove all input edges below the given cutoff (after optional normalization).</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>-1</td>
 *	<td>Integer</td>
 *	<td>Maximum number of genes in input file before in-memory score caching is disabled.  If -1, caching is
 *		never performed.  Caching greatly speeds up processing, but can consume large amounts of memory
 *		for inputs with many genes (rows).</td>
 * </tr></table>
 */
