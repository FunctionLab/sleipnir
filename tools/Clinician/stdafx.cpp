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
 * \page Clinician Clinician
 * 
 * Clinician performs multiple correlation tests of a clinical (or any non-expression) variable against
 * genomewide expression values.  This can be used to determine transcript correlates of a molecular or
 * clinical phenotype, and HEFalMp/bioPIXIE queries can be used as a pre-screen to mitigate the effects
 * of multiple hypothesis testing.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * PCLPlotter -i <data.pcl>
 * \endcode
 * 
 * Produce a list of clinical correlates in the specially formatted \c data.pcl, which should have exactly
 * one non-data column after the initial ID column (in place of the standard NAME/GWEIGHT columns).  This
 * column should contain 0 for standard expression values and 1 for clinical variables with which they are
 * to be correlated.
 * 
 * \code
 * PCLPlotter -i <data.pcl> -I <network.dab>
 * \endcode
 * 
 * Produce a list of clinical correlates in the specially formatted \c data.pcl, formatted as described
 * above, but pre-screen the correlates using the interaction network \c network.dab.  For each clinical
 * variable, a small number of top correlates will be pre-selected and used as a HEFalMp/bioPIXIE query into
 * the given interaction network.  Only the nearest neighbors from this query will be tested for significant
 * clinical correlation, reducing the number of necessary multiple hypothesis tests.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Clinician/Clinician.ggo
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
 *	<td>Input PCL file from which expression and clinical variables are read.  Must be formatted with exactly
 *		one, rather than the standard two, non-ID columns.  Rows containing a 0 in this column will be treated
 *		as gene expression, rows containing a 1 will be treated as clinical correlates.</td>
 * </tr><tr>
 *	<td>-I</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>If given, input DAT/DAB file used to pre-screen potential clinical correlates.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>100</td>
 *	<td>Integer</td>
 *	<td>If given, number of top correlates used during pre-screening as a HEFalMp/bioPIXIE query.</td>
 * </tr><tr>
 *	<td>-N</td>
 *	<td>1000</td>
 *	<td>Integer</td>
 *	<td>If given, number of neighbors retrieved from a HEFalMp/bioPIXIE query to reduce multiple hypothesis
 *		testing.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>On</td>
 *	<td>Flag</td>
 *	<td>If given, perform a HEFalMp rather than bioPIXIE query.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>1</td>
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
