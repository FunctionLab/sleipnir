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
 * \page Orthologer Orthologer
 * 
 * Orthologer will read in a text-based orthology file containing orthologous clusters of functionally
 * equivalent genes across multiple organisms.  It will then "align" gold standards for those organisms in
 * two ways.  First, gene pairs from the same organism within an orthologous cluster will be assigned one
 * weight (high) in the aligned standard.  Next, gene pairs in from the same organism spanning different
 * orthologous clusters that share a gene pair related in some organism will be assigned a second weight
 * (low).  Existing positive functional relationships in the input files are never overwritten.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Orthologer -i <orthology.txt> -w <weight1> -W <weight2> <answers.dab>*
 * \endcode
 * 
 * Aligns (modifies) the organism-specific answer files \c answers.dab based on the orthologous clusters in
 * \c orthology.txt, assigning score \c weight1 to previously unrelated genes sharing an orthologous cluster
 * and \c weight2 to previously unrelated genes in different orthologous clusters spanned by a gene pair
 * related in some organism.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Orthologer/Orthologer.ggo
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
 * 	<td>Input gold standard answer files for each organism within the orthology to be aligned.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Orthology text file</td>
 *	<td>Tab-delimited text file in which each line represents an orthologous cluster.  Tab-separated tokens
 *		within each line should be of the form &lt;organism id>|&lt;gene id>, e.g. SCE|YAL001C or
 *		MMU|MGI:88064.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Weight inserted into aligned gold standards for genes in the same orthologous cluster.</td>
 * </tr><tr>
 *	<td>-W</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Weight inserted into aligned gold standards for genes in different orthologous clusters spanned by
 *		a gene pair functionally related in at least one organism.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
