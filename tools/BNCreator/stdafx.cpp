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
 * \page BNCreator BNCreator
 * 
 * BNCreator will construct a naive Bayesian classifier from data and learn its parameters or, given an
 * existing classifier and data, evaluate the classifier to predict probabilities of functional relationship.
 * The can be performed extremely efficiently and in a context-specific manner.
 * 
 * \section sec_overview Overview
 * 
 * Given one or more biological datasets (stored as Sleipnir::CDat objects in DAT/DAB/etc. files) and a
 * functional gold standard, BNCreator will construct a naive Bayesian classifier to probabilistically
 * integrate the given data.  This implies the construction of a Bayesian network with a single class
 * node (corresponding to functional relationships drawn from the gold standard) with one child node per
 * dataset.  The values taken by these child nodes are determined from the discretization (QUANT files)
 * of the given datasets.
 * 
 * BNCreator behaves similarly to \ref BNConverter, but in a much simpler and more efficient manner.
 * Moreover, <a href="http://www.ncbi.nlm.nih.gov/pubmed/17369653">Huttenhower et al 2006</a> has shown
 * that more complex models (such as those learned by \ref BNConverter) have essentially no benefits for
 * integration of this type, so BNCreator is generally the way to go.
 * 
 * Suppose you have a directory containing data files of the form:
 * \code
 * MICROARRAY.dab
 * MICROARRAY.quant
 * TF.dab
 * TF.quant
 * CUR_COMPLEX.dab
 * CUR_COMPLEX.quant
 * ...
 * SYNL_TRAD.dab
 * SYNL_TRAD.quant
 * \endcode
 * Each data file is a Sleipnir::CDat, either a DAT or a DAB, containing experimental results.  Each QUANT
 * file describes how to discretize that data for use with the Bayesian network (the number of bins in the
 * QUANT must be the same as the number of values taken by the corresponding node in the Bayesian network).
 * Once we've placed all of these files in a directory (e.g. \c ./data/) and assembled a functional gold
 * standard (e.g. \c ANSWERS.dab, possibly constucted by \ref Answerer), we can learn a naive classifier
 * as follows:
 * \code
 * BNCreator -w ANSWERS.dab -o learned.xdsl ./data/*.dab
 * \endcode
 * This produces a "learned" classifier with probabilities that model (as accurately as possible) the
 * relationship between the given data and functional gold standard.  This model can be used to predict
 * functional relationships between new genes not in the standard (but with experimental data) by
 * running BNCreator again in evaluation mode:
 * \code
 * BNCreator -i learned.xdsl -o predicted_relationships.dab -d ./data/
 * \endcode
 * 
 * The \c predicted_relationships.dab file now containins a Sleipnir::CDat in which each pairwise score
 * represents a probability of functional relationship, and it can be mined with tools such as \ref Dat2Dab
 * or \ref Dat2Graph.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNCreator -w <answers.dab> -o <learned.xdsl> <files.dab>*
 * \endcode
 * 
 * Construct a naive classifier, learn its parameters, and store these in \c learned.xdsl, based on one
 * or more given data files (\c files.dab, which must have associated QUANT files) and the functional gold
 * standard in \c answers.dab.
 * 
 * \code
 * BNCreator -d <data_dir> -i <learned.xdsl> -o <predictions.dab>
 * \endcode
 * 
 * Saves predicted probabilities of functional relationships in \c predictions.dab, based on the parameters
 * in the classifier \c learned.xdsl and the data in \c data_dir.
 * 
 * More realistically, a classifier can be learned and evaluated as:
 * \code
 * BNCreator -w <answers.dab> -o <learned.xdsl> -b <defaults.xdsl> -Z <zeros.txt> -m <files.dab>*
 * BNCreator -d <data_dir> -i <learned.xdsl> -o <predictions.dab> -Z <zeros.txt> -m
 * \endcode
 * 
 * This reads default data values for missing gene pairs from \c zeros.txt and default probability
 * distributions for sparse datasets from \c defaults.xdsl, and it memory maps data files using the \c -m
 * flag for increased efficiency.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNCreator/BNCreator.ggo
 * 
 * <table><tr>
 * 	<th>Flag</th>
 * 	<th>Default</th>
 * 	<th>Type</th>
 * 	<th>Description</th>
 * </tr><tr>
 * 	<td>None</td>
 * 	<td>None</td>
 * 	<td>DAT/DAB files</td>
 * 	<td>Datasets used when learning a Bayesian classifier.  One node will be created in the learned
 *		classifier for each DAT/DAB file on the command line, with the ID of the given filename (which
 *		should be alphanumeric).  Must be accompanied by appropriate QUANT files.</td>
 * </tr><tr>
 * 	<td>-w</td>
 * 	<td>None</td>
 * 	<td>DAT/DAB file</td>
 * 	<td>Functional gold standard for learning.  Should consist of gene pairs with scores of 0 (unrelated),
 * 		1 (related), or missing (NaN).</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>Naive classifier for evaluation.  BNCreator will look for data files in the given directory
 *		with filenames corresponding to the given network's node IDs.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>(X)DSL or DAT/DAB file</td>
 *	<td>During learning, (X)DSL file into which the Bayesian classifier and learned parameters are stored.
 *		During evaluation, DAT or DAB file in which predicted probabilities of functional relationship are
 *		saved.</td>
 * </tr><tr>
 * <td>-d</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Used only during evaluation as the directory containing data files.  Must be DAB, DAT, or DAS files
 *		with associated QUANT files and names corresponding to the network node IDs.  For learning, files
 *		must be given directly on the command line.</td>
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
 *	<td>-C</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs passing an "edge" filter against the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>If present during learning, parameters from the given (X)DSL file are used instead of learned
 *		parameters for probability tables with too few examples.  For details, see
 *		Sleipnir::CBayesNetSmile::SetDefault.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, assume that all missing gene pairs in all datasets have a value of 0 (i.e. the first bin).</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>None</td>
 *	<td>Tab-delimited text file</td>
 *	<td>If given, argument must be a tab-delimited text file containing two columns, the first node
 *		IDs (e.g. \c MICROARRAY or \c TF in the example above) and the second bin numbers (zero indexed).
 *		For each node ID present in this file, missing values will be substituted with the given bin number.
 *		For example, if a zeros file contained the line <tt>TF	2</tt>, each missing gene pair in the
 *		\c TF data file (probably \c TF.dab) would be assumed to have the discretized value 2 during
 *		learning/evaluation.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip in any PCL data files between the initial ID column and the experimental
 *		data columns.  Must be the same number for all PCL files.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, convert any PCL data files to z-scores instead of z-transformed Pearson correlations.
 *		Only used if PCL data files (instead of DAT/DAB/etc.) are present.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>If given, learn multiple context-specific Bayesian classifiers from the given data and answers,
 *		assuming that each file in the given directory is a text gene set specifying a functional
 *		context.  This is much more easily done using \ref BNWeaver.</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, group identical examples into one heavily weighted example.  This greatly improves
 *		efficiency, and there's essentially never a reason to deactivate it.</td>
 * </tr></table>
 */
