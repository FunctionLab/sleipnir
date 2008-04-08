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
 * \page DataDumper DataDumper
 * 
 * DataDumper opens an answer file, multiple datasets, and a Bayesian network and displays the data exactly
 * as it would be provided to the network for Bayesian learning or evaluation.  This provides information on
 * exactly what genes/gene pairs and data are being used for learning under different circumstances.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * DataDumper -w <answers.dab> <data.dab>*
 * \endcode
 * 
 * Output (to standard output) the discretized answers for all gene pairs in \c answers.dab with data from
 * the DAT/DAB files \c data.dab (and associated QUANT files) exactly as it would be used in Bayesian
 * learning.  In combination with the other command line arguments, this allows a user to see exactly what
 * data is being used for learning/evaluation under specific circumstances, e.g. context-specific learning,
 * a holdout test set, etc.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include DataDumper/DataDumper.ggo
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
 * 	<td>Input DAT/DAB files from which data is drawn for display in the output.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Functional gold standard for learning.  Should consist of gene pairs with scores of 0 (unrelated),
 *		1 (related), or missing (NaN).</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which both genes are in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which neither gene is in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs passing a "term" filter against the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, assume that all missing gene pairs in all datasets have a value of 0 (i.e. the first bin).</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>None</td>
 *	<td>Tab-delimited text file</td>
 *	<td>If given, argument must be a tab-delimited text file containing two columns, the first node
 *		IDs and the second bin numbers (zero indexed).  For each node ID present in this file, missing values
 *		will be substituted with the given bin number.</td>
 * </tr></table>
 */
