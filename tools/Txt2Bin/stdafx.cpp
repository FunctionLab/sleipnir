#include "stdafx.h"

/*!
 * \page Txt2Bin Txt2Bin
 * 
 * Txt2Bin converts between SVM Light compatible text feature files, Sleipnir binary feature files
 * containing the same information (Sleipnir::CSVM::Learn and Sleipnir::CSVM::Evaluate), and DAT/DAB
 * files containing data.  This is primarily useful for condensing multiple DAT/DAB files into a binary
 * SVM example file for learning/evaluation or for displaying the contents of such a file as human-readable
 * text.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Txt2Bin -f dat -t bin -w <answers.dab> -o <examples.bin> <data.dab>*
 * \endcode
 * 
 * Generates the binary SVM examples/features file \c examples.bin using the labels for gene pairs in
 * \c answers.dab and the data (features) for those gene pairs in the one or more \c data.dab files.
 * 
 * \code
 * Txt2Bin -f txt -t bin -i <examples.txt> -o <examples.bin>
 * \endcode
 * 
 * Convert the SVM Light text-based examples file \c examples.txt into a Sleipnir binary SVM examples file
 * \c examples.bin.
 * 
 * \code
 * Txt2Bin -f bin -t txt -i <examples.bin> -o <examples.txt>
 * \endcode
 * 
 * Convert the Sleipnir binary SVM examples file \c examples.bin into an SVM Light text-based examples file
 * \c examples.txt.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Txt2Bin/Txt2Bin.ggo
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
 * 	<td>When converting from DAT/DAB files to SVM features, the data files to be read.  Each gene pair
 *		becomes one example in the output and each DAT/DAB one feature.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>Text or binary feature file</td>
 *	<td>When converting from text or binary feature files, the feature file to be opened.  Can be in
 *		SVM Light text format or Sleipnir binary format.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>stdin</td>
 *	<td>When converting from DAT/DAB data files, the gold standard answer file used to label the
 *		resulting training examples.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>txt</td>
 *	<td>txt, dat, or bin</td>
 *	<td>Format of input file to be converted from.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>bin</td>
 *	<td>txt or bin</td>
 *	<td>Format of output file to be converted to.</td>
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
 * </tr></table>
 */
