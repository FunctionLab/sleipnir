#include "stdafx.h"

/*!
 * \page Clusterer Clusterer
 * 
 * Clusterer performs non-hierarchical clustering, k-means or quality threshhold clustering (QTC), on an
 * input microarray dataset (PCL) using one of Sleipnir's many similarity measures.
 * 
 * \section sec_overview Overview
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Clusterer/Clusterer.ggo
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
 *	<td>Input PCL file of microarray data to be clustered.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>kmeans</td>
 *	<td>kmeans or qtc</td>
 *	<td>Clustering algorithm to be used.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>PCL text file</td>
 *	<td>If given, a PCL file with dimensions equal to the data given with \c -i.  However, the values in the
 *		cells of the weights PCL represent the relative weight given to each gene/experiment pair.  If no
 *		weights file is given, all weights default to 1.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>pearson</td>
 *	<td>pearson, euclidean, kendalls, kolm-smir, spearman, or quickpear</td>
 *	<td>Similarity measure to be used for clustering.  "quickpear" is a simplified Pearson correlation
 *		that cannot deal with missing values or weights not equal to 1.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>10</td>
 *	<td>Integer</td>
 *	<td>For k-means clustering, the desired number of clusters k.  For QTC, the minimum cluster size (i.e.
 *		minimum number of genes in a cluster).</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>0.5</td>
 *	<td>Double</td>
 *	<td>For QTC, the maximum cluster diameter.  Note that this is similarity measure dependent.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>DAT/DAB file</td>
 *	<td>For QTC cocluster threshholding, the output DAT/DAB file to contain the minimum diameter at which
 *		each gene pair coclusters.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>0</td>
 *	<td>Double</td>
 *	<td>For QTC cocluster threshholding, the smallest maximum cluster diameter to consider.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>0</td>
 *	<td>Double</td>
 *	<td>For QTC cocluster threshholding, the size of steps to take between \c -M and \c -m.  Coclustering
 *		is only performed when the given value is nonzero.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, autocorrelate similarity scores (find the maximum similarity score over all possible
 *		lags of the two vectors; see Sleipnir::CMeasureAutocorrelate).</td>
 * </tr></table>
 */
