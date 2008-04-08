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
 * \page BNTester BNTester
 * 
 * BNTester evaluates arbitrarily structured Bayesian networks to produce predicted probabilities of
 * functional relationship; these networks are often learned with \ref BNConverter.
 * 
 * \section sec_overview Overview
 * 
 * Although \ref BNConverter can perform some limited evaluation on the Bayesian networks it learns,
 * BNTester is the primary tool for performing Bayesian inference on arbitrary network structures.
 * Like \ref BNConverter, BNTester can handle networks with complex structure, unobserved (hidden) nodes,
 * or continuous values.
 * 
 * Bayesian integration generally entails assigning one Bayesian network node to each available biological
 * dataset.  Groups of related datasets (e.g. all physical binding datasets) can be collected under a
 * single unobserved "parent" node, and the network is capped by a single Functional Relationship (FR)
 * node representing whether a particular observation (e.g. gene pair) is functionally related.  This
 * process is detailed in <a href="http://www.ncbi.nlm.nih.gov/pubmed/12826619">Troyanskaya et al 2003</a>.
 * 
 * Given such a network (usually as a <a href="http://genie.sis.pitt.edu/">SMILE</a> DSL or XDSL file) with
 * reasonable parameter values and a collection of discretized biological datasets (usually Sleipnir::CDat s
 * stored as DAT or DAB files with associated QUANT files), BNTester will infer probabilities of functional
 * relationship for (optionally) all gene pairs in the genome.  This can be used to evaluate performance for
 * known related/unrelated gene pairs or to predict function and/or functional relationships for
 * uncharacterized genes/gene pairs.
 * 
 * For example, consider the Bayesian network used in
 * <a href="http://www.ncbi.nlm.nih.gov/pubmed/16420673">Myers et al 2005</a>, which we'll assume we've
 * saved as \c biopixie_learned.xdsl with reasonable parameter values (either learned or assigned manually):
 * \image html bn_pixie.png
 * The name of each node is displayed, and SMILE associated an ID with each name; this might be \c FR for
 * FunctionalRelationship, \c COREGULATION for Coregulation, \c MICROARRAY for Microarray Correlation, and
 * so forth.  Each leaf node corresponds to a single dataset, and each non-leaf node is an unobserved
 * (hidden) and has no associated dataset.  To evaluate this network, we should assemble a directory of data
 * files:
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
 * Once we've placed all of these files in a directory (e.g. \c ./data/), we can perform Bayesian inference
 * to evaluate the network:
 * \code
 * BNTester -d ./data/ -i biopixie_learned.xdsl -o predicted_relationships.dab
 * \endcode
 * The \c predicted_relationships.dab file now containins a Sleipnir::CDat in which each pairwise score
 * represents a probability of functional relationship, and it can be mined with tools such as \ref Dat2Dab
 * or \ref Dat2Graph.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNConverter -d <data_dir> -i <network.xdsl> -o <learned.xdsl> -w <answers.dab> -t <frac>
 *		-e <test_predictions.dab> -E <train_predictions.dab>
 * \endcode
 * 
 * Saves learned parameters for the network \c network.xdsl in the new network \c learned.xdsl, based on
 * the data in \c data_dir (containing files with names corresponding to the network node IDs) and the
 * functional gold standard in \c answers.dab.  Hold \c frac fraction of the gene pairs out of training; store
 * predicted probabilities of functional relationship for these pairs in \c test_predictions.dab and the
 * remaining inferred probabilities in \c train_predictions.dab.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNConverter/BNConverter.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * <td>-d</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Directory containing data files.  Must be DAB, DAT, DAS, or PCL files with associated QUANT files
 *		(unless a continuous network is being learned) and names corresponding to the network node IDs.</td>
 * </tr><tr>
 * <td>-D</td>
 *	<td>None</td>
 *	<td>DAD file</td>
 *	<td>DAD file containing data and/or answers for Bayesian learning or evaluation.  Generally constructed
 *		using \ref Dab2Dad.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>File from which Bayesian network structure and/or parameters are determined.  During learning,
 *		only the structure is used; during evaluation, both structure and parameters are used.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>(X)DSL or DAT/DAB file</td>
 *	<td>During learning, (X)DSL file into which a copy of the Bayesian network with learned parameters is
 *		stored.  During evaluation, DAT or DAB file in which predicted probabilities of functional
 *		relationship are saved.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Functional gold standard for learning.  Should consist of gene pairs with scores of 0 (unrelated),
 *		1 (related), or missing (NaN).</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>If present during learning, parameters from the given (X)DSL file are used instead of learned
 *		parameters for probability tables with too few examples.  For details, see
 *		Sleipnir::CBayesNetSmile::SetDefault.</td>
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
 *	<td>-a</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, randomize all parameters before learning or evaluation.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>None</td>
 *	<td>Integer</td>
 *	<td>If given, randomize the parameters of the network node at the given index.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>0</td>
 *	<td>Double</td>
 *	<td>Fraction of available gene pairs to randomly withhold from training and use for evaluation.</td>
 * </tr><tr>
 *	<td>-E</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>If given, save predicted probabilities of functional relationship for the training gene pairs in the
 *		requested file.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>If given, save predicted probabilities of functional relationship for the test gene pairs in the
 *		requested file.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, assume that all missing gene pairs in all datasets have a value of 0 (i.e. the first bin).</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, use the Extended Logistic Regression (ELR) algorithm for learning (due to
 *		<a href="http://www.cs.ualberta.ca/~greiner/ELR/">Greiner and Zhou 2005</a>) in place of EM.  This
 *		will learn a discriminative model, whereas EM will learn a generative one.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, use Intel's <a href="http://www.intel.com/technology/computing/pnl/">PNL</a> library for
 *		Bayesian network manipulation rather than <a href="http://genie.sis.pitt.edu/">SMILE</a>.  Note
 *		that Sleipnir must be compiled with PNL support for this to function correctly!</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, assume the given (X)DSL file represents a custom function-fitting Bayesian network.  For
 *		details, see Sleipnir::CBayesNetFN.</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, group identical examples into one heavily weighted example.  This greatly improves
 *		efficiency, and there's essentially never a reason to deactivate it.</td>
 * </tr></table>
 */
