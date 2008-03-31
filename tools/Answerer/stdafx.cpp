#include "stdafx.h"

/*!
 * \page Answerer Answerer
 * 
 * Answerer generates gold standard DAT or DAB files for machine learning or evaluation.  These are usually
 * DAB files in which each functionally related gene pair is given the value 1, each unrelated gene pair
 * the value 0, and any uncertain gene pairs are left with missing values (NaNs).
 * 
 * \section sec_overview Overview
 * 
 * Given sets of known related genes - pathways, complexes, GO terms, etc. - Answerer generates a gold
 * standard Sleipnir::CDat.  By considering every pair of genes coannotated to one of these sets to be
 * related, the gold standard answers will include a collection of known functionally related pairs.
 * If Answerer is only provided with positive gene sets, it splits the world of gene pairs up into two
 * values, positive (related) and negative (unrelated), modulo any uncertain pairs introduced by the
 * \c overlap option (see below).
 * 
 * In addition to these positive gene sets, Answerer optionally also takes one or more negative sets.
 * These represent "minimally related" genes, such that gene pairs \e not coannotated to a negative set
 * are definitely unrelated.  This is intended to give two threshholds: positive gene sets should be fairly
 * specific (e.g. "mitotic cell cycle" or "aldehyde metabolism" in GO), such that genes coannotated to
 * these processes are known to be doing something biologically similar.  Negative gene sets should be
 * fairly general (e.g. "physiological process" or "translation" in GO), such that any genes not similar
 * enough to be coannotated at this level are definitely unrelated.  Any genes coannotated with specificity
 * "between" these levels (i.e. above the positive level but below the negative level) are uncertain and
 * not included in Answerer's gold standard.
 * 
 * For example, suppose Answerer got the following two positive sets:
 * \code
 * A
 * B
 * C
 * \endcode
 * and
 * \code
 * A
 * D
 * E
 * \endcode
 * and one negative set:
 * \code
 * A
 * B
 * C
 * E
 * \endcode
 * Then it would generate the answer file:
 * \code
 * A	B	1
 * A	C	1
 * A	D	1
 * A	E	1
 * B	C	1
 * B	D	0
 * C	D	0
 * D	E	1
 * \endcode
 * Note that the pairs B E and C E are missing from this answer file: they are neither positive nor negative,
 * since B, C, and E are all coannotated to a negative set.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Answerer -p <positives_dir> [-n <negatives_dir>] -o <answers.dab>
 * \endcode
 * 
 * Generates the answer file \c answers.dab from the positives gene sets (text files, one gene per line) in
 * \c positives_dir and, optionally, the negative sets in \c negatives_dir.  If \c -o is omitted, answers
 * are saved as a DAT on standard output.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Answerer/Answerer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * <td>-p</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Directory containing related (positive) gene lists.  Each gene list is a text file containing one
 * 		systematic gene ID per line.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>File containing known related pairs.  If given, positive gene pairs will be drawn directly from the
 * 		given Sleipnir::CDat rather than calculated from coannotation to gene sets.  Any gene pair with a
 *		non-missing, non-zero value in the given Sleipnir::CDat will become a positive pair.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Directory containing unrelated (negative) gene lists.  Each gene list is a text file containing one
 * 		systematic gene ID per line.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>None</td>
 *	<td>Double</td>
 *	<td>Expected number of positive functional relationships per gene.  If given, negative gene pairs are
 *		chosen at random from the non-positive pairs, with probability equal to the prior of functional
 *		relationship times the size of the genome divided by the requested number of positive
 *		interactions.  For example, if yeast has 6000 genes and you want a gene pair to have a 5% chance
 *		of being functionally related, choose an interaction number 0.05 * 6000 = 300.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Gene file</td>
 *	<td>Text file containing gene IDs, one per line.  Only genes in this list will be used; gene pairs
 *		containing genes not in the list will be ignored.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>0</td>
 *	<td>Double</td>
 *	<td>If nonzero, this fraction of the genome is omitted as genes for future holdout (exclusion) sets.
 *		In practice, this will omit to standard out (if \c -o is given) or standard error (if it is not)
 *		a list of the requested number of genes.  These can be saved to a file and later used as a
 *		holdout/test set.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>0</td>
 *	<td>Double</td>
 *	<td>If nonzero, genes coannotated to positive sets with hypergeometric p-value of overlap less than this
 *		value will be considered uncertain (i.e. missing, NaN) in the output answers instead of unrelated.
 *		In other words, if genes A and B are annotated to two different positive sets, but these two sets
 *		have significant overlap, the gene pair will be neither positive nor negative in the output answer
 *		file.  0.05 is a good value for producing generally sane answer sets.</td>
 * </tr></table>
 */
