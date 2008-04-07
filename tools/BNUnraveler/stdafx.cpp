#include "stdafx.h"

/*!
 * \page BNUnraveler BNUnraveler
 * 
 * BNUnraveler is a multithreaded tool for performing inference on many context-specific Bayesian classifiers
 * to produce predicted functional relationship networks.  It is generally paired with \ref BNWeaver to
 * learn classifiers from data and then infer context-specific functional relationships.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNUnraveler -i <contexts_dir> -d <data_dir> -o <networks_dir> [-t <threads>]
 *		<contexts.txt>*
 * \endcode
 * 
 * For each biological context \c contexts.txt, load a context-specific Bayesian classifier (of the
 * same name) from \c contexts_dir and use Bayesian inference with the data (named identically to the node
 * IDs in the classifier) from \c data_dir to create a context-specific functional relationship
 * network in \c networks_dir; optionally, use \c threads parallel threads.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNUnraveler/BNUnraveler.ggo
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
 *		have been learned.  Must have filenames corresponding to the (X)DSL files to be loaded, e.g. if
 *		"mitotic_cell_cycle.txt" is given on the command line, "mitotic_cell_cycle.xdsl" must exist in
 *		\c -i.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory from which (X)DSL Baysian classifier files are read.  Must be naive classifiers with
 *		identical structure (but presumably different parameters).</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory into which inferred functional relationship networks (DAB files) are placed.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory from which data files are read.  Must be DAT/DAB files with names identical to the
 *		nodes of the Bayesian classifiers.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, predict relationship probabilities for all genes with any data, regardless of context.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>None</td>
 *	<td>DAT/DAB file</td>
 *	<td>If given, predict relationship probabilities only for gene pairs in the given answer file.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, predict relationship probabilities for all gene pairs for which both genes are in the list
 *		(regardless of context).</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, predict relationship probabilities only for gene pairs for which both genes are in the list
 *		(in addition to context filtering).</td>
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
 *	<td>Number of simultaneous threads to use for individual CPT inferences.  Threads are per classifier node
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
