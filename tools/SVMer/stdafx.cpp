#include "stdafx.h"

/*!
 * \page SVMer SVMer
 * 
 * SVMer learns and evaluates support vector machine models from DAT/DAB datasets in a variety of ways.  If
 * given PCL inputs, SVMer will construct one example per gene pair by concatenating the two genes'
 * expression vectors to create features.  If given DAT/DAB inputs, SVMer will construct one example per
 * gene pair using each dataset as a feature.  In genewise mode, SVMer will learn one model per gene,
 * with one example constructed for each pair in which that gene participates.
 * 
 * \section sec_overview Overview
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include SVMer/SVMer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * 	<td>None</td>
 * 	<td>None</td>
 * 	<td>PCL or DAT/DAB files</td>
 * 	<td>Input data files from which features are constructed, either PCLs from which expression vectors are
 *		concatenated or DAT/DABs from which pairwise values are read.</td>
 * </tr><tr>
 * 	<td>-i</td>
 * 	<td>stdin</td>
 * 	<td>DAT/DAB file</td>
 * 	<td>If given, functional gold standard for learning.  Should consist of gene pairs with scores of 0
 *		(unrelated), 1 (related), or missing (NaN).  If not given, evaluation is assumed and SVM model(s)
 *		is/are read from \c -m.</td>
 * </tr><tr>
 * 	<td>-o</td>
 * 	<td>stdout</td>
 * 	<td>DAT/DAB file</td>
 * 	<td>Output predictions from the SVM model(s) for each available gene pair.</td>
 * </tr><tr>
 * 	<td>-m</td>
 * 	<td>None</td>
 * 	<td>SVM model file or directory</td>
 * 	<td>In standard mode, output learned SVM model file (if \c -i is given) or input SVM model file to be
 *		evaluated (if it is not).  If genewise mode, directory containing output learned or input evaluated
 *		SVM model files.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, assume input files are PCLs from which features are constructed by concatenation of
 *		expression vectors.  If off, assume input files are DAT/DABs from which one feature is drawn per
 *		dataset for each gene pair example.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>None</td>
 *	<td>Binary feature file</td>
 *	<td>If given, ignore other inputs and assume the given binary file is to be used for model evaluation
 *		(if \c -o is specified) or learning (if it is not).</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, learn/evaluate one SVM model per gene, using only the gene pairs including that gene (and thus
 *		each example represents one other gene).  If off, learn/evaluate one global SVM model in which each
 *		feature represents a gene pair.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>None</td>
 *	<td>Gene text file</td>
 *	<td>If given, in genewise mode, learn/evaluate models only for genes in the given gene set.</td>
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
 *	<td>-k</td>
 *	<td>linear</td>
 *	<td>linear, poly, or rbf</td>
 *	<td>SVM kernel type: linear, polynomial, or radial basis function.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>40</td>
 *	<td>Integer (MB)</td>
 *	<td>SVM cache size in megabytes.</td>
 * </tr><tr>
 *	<td>-C</td>
 *	<td>None</td>
 *	<td>Float</td>
 *	<td>SVM tradeoff between misclassification and margin; an appropriate default is calculated if no value
 *		is given.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>1</td>
 *	<td>Float</td>
 *	<td>Gamma parameter for RBF kernel.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>3</td>
 *	<td>Integer</td>
 *	<td>Degree parameter for polynomial kernel.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>None</td>
 *	<td>Alphas file</td>
 *	<td>If given, SVM Light alphas file used to initialize the SVM model.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>100000</td>
 *	<td>Integer</td>
 *	<td>Maximum number of iterations to run per SVM learning epoch.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip in any PCL data files between the initial ID column and the experimental
 *		data columns.  Must be the same number for all PCL files.</td>
 * </tr></table>
 */
