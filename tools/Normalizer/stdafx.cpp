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
 * \page Normalizer Normalizer
 * 
 * Normalizer will perform simple normalization of PCL files (Sleipnir::CPCL::Normalize), transforming each
 * gene's expression vector to have mean zero and standard deviation one, or of DAT/DAB files
 * (Sleipnir::CDat::Normalize), transforming all values to either z-scores or the range [0, 1].
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Normalizer -i <data.dab> -o <normalized.dab> -t dat
 * \endcode
 * 
 * Normalize the given \c data.dab file by z-scoring all values (subtracting the mean and dividing by the
 * standard deviation) and store the result in \c normalized.dab.
 * 
 * \code
 * Normalizer -i <data.pcl> -o <normalized.pcl> -t pcl
 * \endcode
 * 
 * Normalize the given \c data.pcl file by z-scoring each row independently (guaranteeing each row has
 * mean zero, standard deviation one) and store the result in \c normalized.pcl.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Normalizer/Normalizer.ggo
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
 *	<td>Input data file to be normalized, PCLs by row (mean zero/stdev one) and DAT/DABs to [0,1] or
 *		z-scores.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>PCL or DAT/DAB file</td>
 *	<td>Output normalized data file.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>dat</td>
 *	<td>dat or pcl</td>
 *	<td>Type of data file to be normalized.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output one minus the input's values.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to z-scores (subtract mean, divide by standard deviation) before
 *		processing; otherwise, normalize to the range [0,1].</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr></table>
 */
