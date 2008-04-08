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
 * \page MEFIT MEFIT
 * 
 * MEFIT performs all steps necessary to produce context-specific predicted functional relationship
 * networks (DAT/DAB files) from input microarray PCL files as described in Huttenhower et al 2006.  This
 * is essentially a summarization of work performed by \ref Answerer, \ref Distancer, \ref BNCreator, and
 * \ref BNTruster.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * MEFIT -r <related_dir> -u <unrelated.txt> -b <bins.quant> -o <learned_dir> -O <learned.xdsl>
 *		-p <predictions_dir> -t <trusts.txt> <data.pcl>*
 * \endcode
 * 
 * First, construct a gold standard by reading related gene set text files from \c related_dir and unrelated
 * gene pairs from \c unrelated.txt; any two genes coannotated to some set in \c related_dir are considered
 * functionally related, and any gene pair listed in \c unrelated.txt is considered to be unrelated.  Next,
 * compute normalized pairwise correlations for all genes in the given \c data.pcl microarray data files and
 * discretize these scores based on the QUANT file \c bins.quant.  Learn a global Bayesian classifier
 * \c learned.xdsl and context-specific classifiers for each positive gene set, stored in \c learned_dir.
 * Finally, save a table of functional activity scores for each dataset in each context in \c trusts.txt,
 * and save context-specific predicted functional relationship networks in \c predictions_dir.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include MEFIT/MEFIT.ggo
 * 
 * <table><tr>
 * 	<th>Flag</th>
 * 	<th>Default</th>
 * 	<th>Type</th>
 * 	<th>Description</th>
 * </tr><tr>
 * 	<td>None</td>
 * 	<td>None</td>
 * 	<td>PCL text files</td>
 * 	<td>Microarray datasets which will be integrated by MEFIT.  Each dataset will correspond to one node
 *		in each of the learned Bayesian classifiers and assigned a trust score in each biological context.
 *		All input PCLs must have the same number of skip columns \c -s.</td>
 * </tr><tr>
 *	<td>-r</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Input directory containing related (positive) gene lists.  Each gene list is a text file containing
 *		one systematic gene ID per line (see \ref Answerer).</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>None</td>
 *	<td>Gene pair text file</td>
 *	<td>Input tab-delimited text file containing two columns; each line is a gene pair which is known to be
 *		functionally unrelated (e.g. annotated to two different Gene Ontology terms; see \ref Answerer).</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>pearnorm</td>
 *	<td>pearnorm, pearson, euclidean, kendalls, kolm-smir, or spearman</td>
 *	<td>Similarity measure to be used for converting microarray data into pairwise similarity scores.
 *		\c pearnorm is the recommended Fisher's z-transformed Pearson correlation.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>None</td>
 *	<td>QUANT text file</td>
 *	<td>Input tab-delimited QUANT file containing exactly one line of bin edges; these are used to discretize
 *		pairwise similarity scores.  For details, see Sleipnir::CDataPair.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Output directory in which learned context-specific Bayesian classifiers are saved as (X)DSL
 *		files (see \ref BNCreator).</td>
 * </tr><tr>
 *	<td>-O</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>Output file in which the learned global (non-context-specific) Bayesian classifier is saved
 *		(see \ref BNCreator).</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Directory in which predicted context-specific functional relationships (DAT/DAB files) are saved
 *		(see \ref BNCreator).</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>None</td>
 *	<td>PCL text file</td>
 *	<td>Output PCL file in which dataset/context functional activity scores are saved (see
 *		\ref BNTruster).</td>
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
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, assume that all missing gene pairs in all datasets have a value of 0 (i.e. the first bin).</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>None</td>
 *	<td>Double</td>
 *	<td>If given, remove all input edges below the given cutoff (after optional normalization).</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, assume XDSL files will be used instead of DSL files.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, output DAB files instead of DAT files.</td>
 * </tr></table>
 */
