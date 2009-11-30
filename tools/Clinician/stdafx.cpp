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
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"

/*!
 * \page PCLPlotter PCLPlotter
 * 
 * PCLPlotter produces summary statistics from a PCL file and is optimized to show the mean expression values
 * for subsets (biclusters) of genes and conditions.  It can also provide accompanying bicluster statistics
 * for associated FASTA files containing gene sequences.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * PCLPlotter -i <data.pcl>
 * \endcode
 * 
 * Produce a summary of the mean and standard deviation expression for each condition in \c data.pcl.  If a
 * bicluster is present, member genes should be marked with an initial \c * in their NAME (not ID) column, and
 * member conditions should be marked with an initial \c *.
 * 
 * \code
 * PCLPlotter -i <cluster.pcl> -b <genome.pcl>
 * \endcode
 * 
 * Produce a summary of the mean and standard deviation expression for each condition in \c cluster.pcl and
 * \c genome.pcl, which must contain the same conditions.  Genes in \c cluster.pcl are considered to be
 * members of the bicluster, and member conditions should be marked with an initial \c *.
 * 
 * \code
 * PCLPlotter -i <data.pcl> -g <genes.txt>
 * \endcode
 * 
 * Produce a summary of the mean and standard deviation expression for each condition in \c data.pcl.
 * Genes in \c genes.txt are considered to be members of the bicluster, and member conditions should be marked
 * with an initial \c *.
 * 
 * \code
 * PCLPlotter -i <data.pcl> -f <data.fasta>
 * \endcode
 * 
 * Produce a summary of the mean and standard deviation expression for each condition in \c data.pcl, as well
 * as an HMM summarizing sequence characteristics.  If a bicluster is present, member genes should be marked
 * with an initial \c * in their NAME (not ID) column, and member conditions should be marked with an initial
 * \c *.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include PCLPlotter/PCLPlotter.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>PCL file</td>
 *	<td>Input PCL file from which bicluster summary information is extracted.  In the absence of \c -b or \c -g
 *		options, genes and conditions in the bicluster should be marked with a \c * at the beginning of their
 *		NAME and label, respectively.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>None</td>
 *	<td>FASTA file</td>
 *	<td>If given, input FASTA sequence file from which cluster sequence summary information is extracted.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>None</td>
 *	<td>PCL file</td>
 *	<td>If given, input PCL file from which non-bicluster summary information is extracted; all genes in \c -i
 *		are considered to be in the bicluster, and conditions should be marked with a \c *.  PCL files for \c -i
 *		and \c -b should contain exactly the same conditions.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, input text file from which biclustered genes are read; other genes in \c -i are considered to
 *		be out of the bicluster.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>None</td>
 *	<td>Text motif list</td>
 *	<td>If given, input text file from which known motifs are read.  In conjunction with \c -f, frequencies
 *		of each motif in gene sequences in and out of the bicluster will be provided.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>7</td>
 *	<td>Integer (base pairs)</td>
 *	<td>Default number of base pairs per motif; largely unrelated to the contents of \c -m.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>Degree of HMM used to provide summary statistics of sequences given in \c -f.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr></table>
 */
