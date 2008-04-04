#include "stdafx.h"

/*!
 * \page BNEvaluator BNEvaluator
 * 
 * BNEvaluator will perform Bayesian inference on a given network by reading feature (node) values directly
 * from a PCL file.  Rather than accessing a collection of biological datasets stored as separate files,
 * BNEvaluator reads examples from rows of a PCL file, assumes each experimental column is a feature
 * (node), discretizes the data in each cell with an associated QUANT file, and infers probabilities for
 * each possible value of a single class node (usually binary).
 * 
 * \section sec_overview Overview
 * 
 * BNEvaluator can be used to perform Bayesian inference on arbitrary Bayes networks using data read
 * directly from a PCL File (see Sleipnir::CPCL)  Unlike most of Sleipnir's Bayesian analysis tools,
 * which assume each gene pair is an example and each dataset a feature, BNEvaluator is agnostic to the
 * form and semantics of the PCL's data: each row is an example, each experimental column a feature, and
 * values in the cells are discretized using a single associated QUANT file (e.g. \c data.quant for a
 * data file \c data.pcl ).
 * 
 * For example, suppose we have a data file \c data.pcl containing:
 * \code
 * ID	DATA1	DATA2	DATA3
 * THING1	1.1	2.2	3.3
 * THING2	2.1	3.2	1.3
 * THING3	3.0	2.0	1.0
 * \endcode
 * In the same directory, we have a QUANT file containing:
 * \code
 * 1.5	2.5	3.5	4.5
 * \endcode
 * Finally, we have a Bayesian network \c network.xdsl with reasonable parameters and the structure:
 * \image html bn_evaluator.png
 * Each of the three data nodes can take four values (to agree with our QUANT file), and the class node
 * \c FR is binary.
 * 
 * We now run (noting that our \c data.pcl has zero skip columns):
 * \code
 * BNEvaluator -i network.xdsl -d data.pcl -s 0 -o results.pcl
 * \endcode
 * This will evaluate the probability of the \c FR class node's values using the three examples [0, 1, 2],
 * [1, 2, 0], and [2, 1, 0] for the data nodes and produce a \c results.pcl file that resembles:
 * \code
 * ID	FR0	FR1
 * THING1	0.3	0.7
 * THING2	0.9	0.1
 * THING3	0.4	0.6
 * \endcode
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNEvaluator -i <network.xdsl> -d <data.pcl> -o <results.pcl>
 * \endcode
 * 
 * Saves inferred probabilities for the class (\c FR) node from the network \c network.xdsl in the
 * output \c results.pcl file, based on the data (columns) for each example (row) in \c data.pcl.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNEvaluator/BNEvaluator.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>File from which Bayesian network structure and parameters are determined.  Columns of the input PCL
 *		are mapped by name to node IDs in the given network.</td>
 * </tr><tr>
 * <td>-d</td>
 *	<td>None</td>
 *	<td>PCL text file</td>
 *	<td>PCL file containing examples (rows) from which data features (columns) are read and quantized for
 *		use during evaluation of the given Bayes net.  For a non-continuous Bayes net, must have an
 *		associated QUANT file of the same name (e.g. \c data.quant in the same location as a data file
 *		\c data.pcl).</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>PCL text file</td>
 *	<td>PCL file in which inferred class probabilities are saved.  Each row is an example from the input
 *		file, each column a possible value of the class node (zero indexed).</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only genes in the list.</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene not in the list.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, assume that all missing gene pairs in all datasets have a value of 0 (i.e. the first bin).</td>
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
 *	<td>-a</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>ID of Bayesian inference algorithm to use (passed directly to SMILE).  The default is almost
 *		always the most efficient and accurate option.</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, group identical examples into one heavily weighted example.  This greatly improves
 *		efficiency, and there's essentially never a reason to deactivate it.</td>
 * </tr></table>
 */
