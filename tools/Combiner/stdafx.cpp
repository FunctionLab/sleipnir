#include "stdafx.h"

/*!
 * \page Combiner Combiner
 * 
 * Combiner combines things!  Where things include combining PCLs into one large PCL, DAT/DABs into one
 * averaged DAT/DAB, or DAT/DABs into a DAD dataset file.  Multiple PCLs are combined into a single wide
 * PCL by aligning genes' expression vectors and inserting missing values for genes not present in some
 * dataset(s).  DAT/DABs are combined into a single DAT/DAB by averaging pairwise scores in a configurable
 * manner.  DAT/DABs are combined into a DAD by ordering their individual scores as detailed in
 * Sleipnir::CDatasetCompact::Save.
 * 
 * \section sec_overview Overview
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Combiner -t pcl -o <combined.pcl> <data.pcl>*
 * \endcode
 * 
 * Create a new PCL file \c combined.pcl containing all genes in the microarray PCL files \c data.pcl,
 * with new expression vectors consisting of the concatenation of all data from these input files.  In other
 * words, take the input PCLs, line up each gene's values, smoosh them all together, and plop the result
 * into the output file.
 * 
 * \code
 * Combiner -t dat -o <combined.dab> -n <data.dab>*
 * \endcode
 * 
 * Create a new DAT/DAB file \c combined.dab in which each gene pair's score is the average of the
 * normalized (by z-scoring) scores from all input DAT/DAB files \c data.dab.  The combination method can be
 * modified using \c -m.
 * 
 * \code
 * Combiner -t dab -o <combined.dad> <data.dab>*
 * \endcode
 * 
 * Create a new DAD file \c combined.dad containing the discretized gene pair scores from all input DAT/DAB
 * files \c data.dab, which must be accompanied by appropriate QUANT files.  This is equivalent to
 * \ref Dab2Dad.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Combiner/Combiner.ggo
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
 * 	<td>Input files to be combined; must all be of an appropriate type for the requested output.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>pcl</td>
 *	<td>pcl, dat, or dab</td>
 *	<td>Type of combination to perform: \c pcl combines PCLs into a PCL, \c dat combines DAT/DABs into a
 *		DAT/DAB, and \c dad combines DAT/DABs into a DAD.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>mean</td>
 *	<td>mean, min, max, gmean, or hmean</td>
 *	<td>Type of DAT/DAB combination to perform when combining pairwise scores.  Options are to calculate the
 *		arithmetic, geometric, or harmonic mean, or to retain only the minimum or maximum value for each
 *		gene pair.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>PCL, DAT/DAB, or DAD file</td>
 *	<td>Output file of the type specified by \c -t.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>If nonzero, process input DAT/DABs in subsets of the requested size as described in
 *		Sleipnir::CDataSubset.</td>
 * </tr></table>
 */
