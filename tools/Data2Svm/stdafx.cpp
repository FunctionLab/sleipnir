#include "stdafx.h"

/*!
 * \page Data2Svm Data2Svm
 * 
 * Data2Svm learns a support vector machine model classifying individual genes in or out (positive or negative)
 * of a given gene set.  Constructs features for each example (gene) based on data in an input PCL.  Similar
 * to \ref SVMer.
 * 
 * \section sec_overview Overview
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Data2Svm -i <data.pcl> -m <learned.svm> -g <context.txt> -G <holdout.txt>
 * \endcode
 * 
 * Learn a support vector machine model (saved as \c learned.svm) using the microarray expression values in
 * \c data.pcl as features, labeling the genes in \c context.txt as positive examples (and all other genes
 * as negatives), and holding out the genes in \c holdout.txt from training.  Outputs (to standard output)
 * the predicted SVM classifications for all genes after learning.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Data2Svm/Data2Svm.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>PCL text file</td>
 *	<td>Input PCL file from which features will be drawn to construct SVM examples.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>stdout</td>
 *	<td>SVM model file</td>
 *	<td>Output learned SVM model.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Gene text file</td>
 *	<td>Set of genes to be labeled as positive examples.</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>None</td>
 *	<td>Gene text file</td>
 *	<td>If given, set of genes to be held out of training and evaluated as test examples.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, evaluate and output SVM predictions only for test genes; if off, evaluate all genes.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, randomize input feature values within each row (gene).</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, randomize output SVM prediction labels across all genes.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>40</td>
 *	<td>Integer (MB)</td>
 *	<td>SVM cache size in megabytes.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>linear</td>
 *	<td>linear, poly, or rbf</td>
 *	<td>SVM kernel type: linear, polynomial, or radial basis function.</td>
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
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to z-scores (subtract mean, divide by standard deviation) before
 *		processing.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr></table>
 */
