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
 * \page SVDer SVDer
 * 
 * SVDer performs singular value decomposition of the matrix of expression values in a given PCL file
 * (which should contain no missing values; see \ref KNNImputer).  The PCL file can then be transformed in
 * SVD space using a number of algorithms, e.g. up- or down-weighting low-variance singular values or
 * continuously transforming the weights of each value.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * SVDer -i <input.pcl> -o <output.pcl> -r <fraction>
 * \endcode
 * 
 * Decompose the input matrix \c input.pcl, retain only the minimum number of singular values necessary to
 * account for \c fraction of the total variance, set the rest to zero, and reproject the resulting
 * matrix in \c output.pcl.
 * 
 * \code
 * SVDer -i <input.pcl> -o <output.pcl> -b
 * \endcode
 * 
 * Decompose the input matrix \c input.pcl, set all singular values to the mean (thus preserving the total
 * variance), and reproject the resulting matrix in \c output.pcl.
 * 
 * \code
 * SVDer -i <input.pcl> -o <output.pcl> -t <parameter>
 * \endcode
 * 
 * Decompose the input matrix \c input.pcl, transform the singular values proportionally to the \c parameter
 * power, and reproject the resulting matrix in \c output.pcl.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include SVDer/SVDer.ggo
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
 *	<td>Input PCL file to be transformed; must not contain any missing values.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>PCL file</td>
 *	<td>Output PCL file containing the transformed version of \c -i.</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>None</td>
 *	<td>PCL file</td>
 *	<td>If given, output PCL file containing the U matrix (basis vectors) of the decomposition of \c -i.</td>
 * </tr><tr>
 *	<td>-r</td>
 *	<td>1</td>
 *	<td>Double (fraction)</td>
 *	<td>Fraction of variance to be reprojected; when less than 1, proportionally many singular values are
 *		set to zero before reprojecting the output matrix.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Power by which singular values are transformed before reprojection.  For singular values {s1, ..., sn}
 *		S = sum(si), ti = (si/S)^\c -t, and T = sum(ti), each si is replaced by ti/sum(ti).</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>Off</td>
 *	<td>Flag</td>
 *	<td>If given, set all singular values to their average value before reprojection; otherwise, transform
 *		them according to \c -t.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr></table>
 */
