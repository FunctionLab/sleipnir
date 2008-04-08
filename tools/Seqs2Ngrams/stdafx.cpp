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
 * \page Seqs2Ngrams Seqs2Ngrams
 * 
 * Seqs2Ngrams reads an input FASTA file and uses it to compute pairwise sequence similarity scores between
 * genes.  This is done simplistically by computing the similarity of each gene pair's vectors of k-mer
 * counts, but could easily be modified to do something more complex.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Seqs2Ngrams -i <sequences.fasta> -o <similarities.dab> -n <ngram>
 * \endcode
 * 
 * Breaks the sequences in \c sequences.fasta into k-mers of size \c ngram, computes k-mer frequency
 * counts for each gene, and outputs gene pair similarities based on comparisons of these frequency vectors
 * in \c similarities.dab.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Seqs2Ngrams/Seqs2Ngrams.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>FASTA text file</td>
 *	<td>Input FASTA file containing gene IDs and sequences.  Identifiers after each > record header should
 *		consist solely of the unique gene ID without any additiona information, e.g. > YAL001C.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>DAT/DAB file</td>
 *	<td>Output DAT/DAB file containing pairwise sequence similarity scores calculated from the given
 *		sequence.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>7</td>
 *	<td>Integer</td>
 *	<td>Size of sequence n-grams to use when calculating pairwise similarity.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize output edges to the range [0,1].</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize output edges to z-scores (subtract mean, divide by standard deviation).</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only genes in the list.</td>
 * </tr></table>
 */
