#include "stdafx.h"

/*!
 * \page MEFIT MEFIT
 * 
 * MEFIT performs all steps necessary to produce context-specific predicted functional relationship
 * networks (DAT/DAB files) from input microarray PCL files as described in Huttenhower et al 2006.  This
 * is essentially a summarization of work performed by \ref Answerer, \ref Distancer, \ref BNCreator, and
 * \ref BNTruster.
 * 
 * \section sec_overview Overview
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
