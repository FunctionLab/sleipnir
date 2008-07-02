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
 * \page BNs2Txt BNs2Txt
 * 
 * BNs2Txt converts a binary representation of a collection of naive Bayesian classifiers
 * (Sleipnir::CBayesNetMinimal) into a set of (X)DSL files.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNs2Txt -i <networks.bin> -o <network_dir> -d <datasets.txt>
 * \endcode
 * 
 * Reads the binary naive Bayesian classifiers stored in \c networks.bin and the text file containing dataset
 * (classifier node) names \c datasets.txt and creates one (X)DSL file per classifier in the output directory
 * \c network_dir.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNs2Txt/BNs2Txt.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Binary network file</td>
 *	<td>Binary naive Bayesian classifiers (Sleipnir::CBayesNetMinimal) stored by another tool
 *		(e.g. \ref Counter).</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory in which (X)DSL files are created.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>Text file containing one dataset name per line.  These correspond to the nodes of each classifier in
 *		the input file and are used to name the equivalent nodes in the output (X)DSL files.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, create XDSL files instead of DSL files.</td>
 * </tr></table>
 */
