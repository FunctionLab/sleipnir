#include "stdafx.h"

/*!
 * \page Seqs2Ngrams Seqs2Ngrams
 * 
 * Seqs2Ngrams reads an input FASTA file and uses it to compute pairwise sequence similarity scores between
 * genes.  This is done simplistically by computing the similarity of each gene pair's vectors of k-mer
 * counts, but could easily be modified to do something more complex.
 * 
 * \section sec_overview Overview
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
