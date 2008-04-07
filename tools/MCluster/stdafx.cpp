#include "stdafx.h"

/*!
 * \page MCluster MCluster
 * 
 * MCluster performs hierarchical clustering of a given microarray dataset based on one of a variety of
 * similarity measures or on a precomputed gene pair similarity matrix (e.g. predicted functional
 * relationships).  Output files are compatible with visualization tools such as Java TreeView.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * MCluster -i <data.pcl> -o <clustered.gtr> > <clustered.cdt>
 * \endcode
 * 
 * Cluster the microarray data in \c data.pcl using the default similarity measure (Pearson correlation);
 * save the clustered expression data in \c clustered.cdt and the gene tree in \c clustered.gtr.
 * 
 * \code
 * MCluster -i <data.dab> -o <clustered.gtr> < <data.pcl> > <clustered.cdt>
 * \endcode
 * 
 * Cluster the microarray data in \c data.pcl using the pre-calculated similarity scores (e.g. probabilities
 * of functional relationship) in \c data.dab; save the clustered expression data in \c clustered.cdt and
 * the gene tree in \c clustered.gtr.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include MCluster/MCluster.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>PCL or DAT/DAB file</td>
 *	<td>Input PCL microarray data or DAT/DAB pairwise similarity data.  If a DAT/DAB is given here, a
 *		PCL file should be provided on standard input.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>GTR text file</td>
 *	<td>Output GTR file for use with visualization tools.  Output CDT is printed to standard output.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>pearson</td>
 *	<td>pearson, euclidean, kendalls, kolm-smir, or spearman</td>
 *	<td>Similarity measure to be used for clustering.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>PCL text file</td>
 *	<td>If given, a PCL file with dimensions equal to the data given with \c -i.  However, the values in the
 *		cells of the weights PCL represent the relative weight given to each gene/experiment pair.  If no
 *		weights file is given, all weights default to 1.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output one minus the input's values.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>None</td>
 *	<td>Double</td>
 *	<td>If given, remove all input edges below the given cutoff (after optional normalization).</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Raise all input similarity scores to the given power.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr></table>
 */
