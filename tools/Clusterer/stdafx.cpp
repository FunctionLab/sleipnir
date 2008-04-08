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
 * \page Clusterer Clusterer
 * 
 * Clusterer performs non-hierarchical clustering, k-means or quality threshhold clustering (QTC), on an
 * input microarray dataset (PCL) using one of Sleipnir's many similarity measures.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Clusterer -i <data.pcl> -k <clusters>
 * \endcode
 * 
 * Output (to standard output) a list of k-means clusters formed from the microarray PCL data \c data.pcl
 * with k equal to \c clusters using the default similarity measure (which can be modified using \c -d).
 * 
 * \code
 * Clusterer -i <data.pcl> -a qtc -k <min_size> -m <max_diameter)
 * \endcode
 * 
 * Output a list of quality threshhold clusters formed from the microarray PCL data \c data.pcl with a
 * minimum cluster size of \c min_size and a maximum cluster diameter of \c max_diameter.
 * 
 * \code
 * Clusterer -i <data.pcl> -o <cocluster.dab> -k <min_size> -M <min_diameter> -m <max_diameter>
 *		-e <delta_diameter>
 * \endcode
 * 
 * Create a DAT/DAB file \c cocluster.dab with gene pair scores indicating the minimum diameter size at
 * which each gene pair from \c data.pcl coclustered, using QTC with a minimum cluster size of \c min_size
 * and testing cluster diameters from \c min_diameter to \c max_diameter by steps fo \c delta_diameter.
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
