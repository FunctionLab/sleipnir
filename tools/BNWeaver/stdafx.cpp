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
 * \page BNWeaver BNWeaver
 * 
 * BNWeaver is a multithreaded tool for learning context-specific naive Bayesian classifiers from datasets
 * (DAT/DAB files).  It is generally paired with \ref BNUnraveler to learn classifiers from data and then
 * infer context-specific functional relationships.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNWeaver -w <answers.dab> -o <contexts_dir> -d <data_dir> [-b <global.xdsl>] [-t <threads>]
 *		<contexts.txt>*
 * \endcode
 * 
 * Using the gold standard answers in \c answers.dab and all DAT/DAB files in \c data_dir, create one
 * context-specific Bayesian classifier in \c contexts_dir for each biological context \c contexts.txt;
 * optionally, use probability tables from \c global.xdsl as fallbacks when insufficient data is available
 * for context-specific learning, and use \c threads parallel threads.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNWeaver/BNWeaver.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * 	<td>None</td>
 * 	<td>None</td>
 * 	<td>Gene text files</td>
 * 	<td>Gene sets representing biological contexts (sets of related genes) for which Bayesian classifiers
 *		will be learned.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>Functional gold standard for learning.  Should consist of gene pairs with scores of 0 (unrelated),
 *		1 (related), or missing (NaN).</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory into which learned naive Bayesian classifiers ((X)DSL files) are placed.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory from which data files are read.  Must be DAT/DAB files with names from which the node IDs
 *		of the Bayesian classifiers can be created.</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which neither gene is in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs including at least one gene from the given set  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, randomly shuffle all data values (by gene pair) before learning.</td>
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
 *		IDs (see \ref BNCreator) and the second bin numbers (zero indexed).  For each node ID present in
 *		this file, missing values will be substituted with the given bin number.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>1</td>
 *	<td>Integer</td>
 *	<td>Number of simultaneous threads to use for individual CPT learning.  Threads are per classifier node
 *		(dataset), so the number of threads actually used is the minimum of \c -t and the number of
 *		datasets.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, assume XDSL files will be used instead of DSL files.</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, group identical examples into one heavily weighted example.  This greatly improves
 *		efficiency, and there's essentially never a reason to deactivate it.</td>
 * </tr></table>
 */
