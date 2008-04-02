#include "stdafx.h"

/*!
 * \page BNTruster BNTruster
 * 
 * BNTruster evaluates the "influence" each node in a Bayesian network has on the posterior classification
 * probability.  This can be used to calculate the weight of each dataset during Bayesian integration.
 * 
 * \section sec_overview Overview
 * 
 * BNTruster offers a collection of ways of evaluating the weight of each node in a Bayesian network.  These
 * methods are all variations on determining how much influence each node has on the classification
 * posterior, i.e. how "different" each dataset's probability distributions are for the different class
 * values.  Most commonly, this is used to answer the question, "Given a Bayesian classifier that I've
 * learned to integrate many datasets, how much weight is given to each dataset?"  If the classifier(s)
 * being evaluated are context-specific, this provides a measure of the functional activity of each dataset
 * within each biological context.
 * 
 * The most principled trust calculation provides, for each dataset, the average weighted change in posterior
 * probability of functional relationship over all values in that dataset.  For a dataset D taking one or
 * more discretized values d, this can be written as:
 * \code
 * trust(D) = sum( P(D = d) * |P(FR) - P(FR|D = d)|, d in D )
 * \endcode
 * That is, the trust or "influence" of dataset D is the sum over all its values of the probability of that
 * value times the difference in posterior when that value is detected.
 * 
 * Suppose we have two datasets, a microarray dataset MA quantized into five bins and a physical binding
 * dataset PB quantized into two bins.  We've learned two context-specific naive Bayesian classifiers
 * for these datasets using \ref BNCreator, one for translation (\c translation.xdsl) and one for MAPK
 * signaling (\c mapk.xdsl):
 * \image html bn_trust_translation.png
 * \image html bn_trust_mapk.png
 * These networks have the following conditional probability tables:
 * <center><table>
 *	<tr>
 *		<th colspan="6">Translation</th>
 *	</tr><tr>
 *		<th></th>
 *		<th>FR</th>
 *		<th colspan="2">MA</th>
 *		<th colspan="2">PB</th>
 *	</tr><tr>
 *		<th>Value</th>
 *		<td></td>
 *		<td>FR 0 (No)</td>
 *		<td>FR 1 (Yes)</td>
 *		<td>FR 0 (No)</td>
 *		<td>FR 1 (Yes)</td>
 *	</tr><tr>
 *		<th>0</th>
 *		<td>0.9</td>
 *		<td>0.1</td>
 *		<td>0.05</td>
 *		<td>0.4</td>
 *		<td>0.3</td>
 *	</tr><tr>
 *		<th>1</th>
 *		<td>0.1</td>
 *		<td>0.3</td>
 *		<td>0.2</td>
 *		<td>0.6</td>
 *		<td>0.7</td>
 *	</tr><tr>
 *		<th>2</th>
 *		<td></td>
 *		<td>0.4</td>
 *		<td>0.3</td>
 *		<td></td>
 *		<td></td>
 *	</tr><tr>
 *		<th>3</th>
 *		<td></td>
 *		<td>0.15</td>
 *		<td>0.3</td>
 *		<td></td>
 *		<td></td>
 *	</tr><tr>
 *		<th>4</th>
 *		<td></td>
 *		<td>0.05</td>
 *		<td>0.15</td>
 *		<td></td>
 *		<td></td>
 *	</tr>
 * </table>
 * <table>
 *	<tr>
 *		<th colspan="6">MAPK Signaling</th>
 *	</tr><tr>
 *		<th></th>
 *		<th>FR</th>
 *		<th colspan="2">MA</th>
 *		<th colspan="2">PB</th>
 *	</tr><tr>
 *		<th>Value</th>
 *		<td></td>
 *		<td>FR 0 (No)</td>
 *		<td>FR 1 (Yes)</td>
 *		<td>FR 0 (No)</td>
 *		<td>FR 1 (Yes)</td>
 *	</tr><tr>
 *		<th>0</th>
 *		<td>0.95</td>
 *		<td>0.1</td>
 *		<td>0.1</td>
 *		<td>0.9</td>
 *		<td>0.2</td>
 *	</tr><tr>
 *		<th>1</th>
 *		<td>0.05</td>
 *		<td>0.3</td>
 *		<td>0.3</td>
 *		<td>0.1</td>
 *		<td>0.8</td>
 *	</tr><tr>
 *		<th>2</th>
 *		<td></td>
 *		<td>0.35</td>
 *		<td>0.3</td>
 *		<td></td>
 *		<td></td>
 *	</tr><tr>
 *		<th>3</th>
 *		<td></td>
 *		<td>0.15</td>
 *		<td>0.2</td>
 *		<td></td>
 *		<td></td>
 *	</tr><tr>
 *		<th>4</th>
 *		<td></td>
 *		<td>0.1</td>
 *		<td>0.1</td>
 *		<td></td>
 *		<td></td>
 *	</tr>
 * </table></center>
 * 
 * First, we find the trust of each dataset in the translation-specific network:
 * \code
 * trust(MA) = sum( P(MA = d) * |P(FR) - P(FR|MA = d)|, d = 0 to 4 )
 *		= P(MA = 0) * |P(FR) - P(FR|MA = 0)| + ... + P(MA = 5) * |P(FR) - P(FR|MA = 4)|
 *		= 0.095*|0.1 - 0.053| + 0.29*|0.1 - 0.069| + 0.39*|0.1 - 0.077| +
 *			0.165*|0.1 - 0.182| + 0.06*|0.1 - 0.25|
 *		= 0.045
 * trust(PB) = sum( P(PB = d) * |P(FR) - P(FR|PB = d)|, d = 0 to 1 )
 *		= P(PB = 0) * |P(FR) - P(FR|PB = 0)| + P(PB = 1) * |P(FR) - P(FR|PB = 1)|
 *		= 0.39*|0.1 - 0.077| + 0.61*|0.1 - 0.115|
 *		= 0.018
 * \endcode
 * In the context of translation, the microarray data is somewhat more informative.  However, in the MAPK
 * signaling network:
 * \code
 * trust(MA) = sum( P(MA = d) * |P(FR) - P(FR|MA = d)|, d = 0 to 4 )
 *		= P(MA = 0) * |P(FR) - P(FR|MA = 0)| + ... + P(MA = 5) * |P(FR) - P(FR|MA = 4)|
 *		= 0.1*|0.05 - 0.05| + 0.3*|0.05 - 0.05| + 0.348*|0.05 - 0.043| +
 *			0.153*|0.05 - 0.066| + 0.1*|0.05 - 0.05|
 *		= 0.0049
 * trust(PB) = sum( P(PB = d) * |P(FR) - P(FR|PB = d)|, d = 0 to 1 )
 *		= P(PB = 0) * |P(FR) - P(FR|PB = 0)| + P(PB = 1) * |P(FR) - P(FR|PB = 1)|
 *		= 0.865*|0.05 - 0.012| + 0.135*|0.05 - 0.296|
 *		= 0.066
 * \endcode
 * Unsurprisingly (since the example was completely cooked), the process of MAPK signaling is very active
 * in our physical binding dataset.
 * 
 * To have BNTruster do all of this hard work for us, you can just run:
 * \code
 * BNTruster translation.xdsl mapk.xdsl
 * \endcode
 * 
 * Other trust calculations include sums, calculated as:
 * \code
 * trust(D) = sum( |P(D = d|FR) - P(D = d|~FR)|, d in D )
 * \endcode
 * and ratios, calculated as:
 * \code
 * trust(D) = log( prod( max{P(D = d|FR), P(D = d|~FR)} / min{P(D = d|FR), P(D = d|~FR)}, d in D ) )
 * \endcode
 * 
 * Finally, trust scores can also be calculated for individual bins (i.e. dataset values).  These are only
 * available when using the posteriors trust calculation method, and the value represents the fraction of
 * possible change away from prior.  Positive values represent an increase from prior to posterior, and
 * negative values a decrease.  That is:
 * \code
 * trust(d) = ( P(FR|D = d) - P(FR) ) / ( ( P(FR|D = d) > P(FR) ) ? P(~FR) : P(FR) )
 * \endcode
 * Note that this does not incorporate the probability of observing \c d, i.e. P(D = d).
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNTruster <bayesnets.xdsl>*
 * \endcode
 * 
 * Writes trust scores (by default, weighted average of difference in posterior) to standard output for
 * all of the (X)DSL files given on the command line.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNTruster/BNTruster.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * <td>-t</td>
 *	<td>posteriors</td>
 *	<td>posteriors, sums, or ratios</td>
 *	<td>Type of trust score to calculate, as described above.</td>
 * </tr><tr>
 * <td>-b</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If true, output individual bins' influence scores as described above.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>1</td>
 *	<td>Integer</td>
 *	<td>Number of simultaneous threads to use for posterior trust calculations.  Threads are per-(X)DSL, so
 *		the number of threads actually used is the minimum of \c -t and the number of input files.</td>
 * </tr></table>
 */
