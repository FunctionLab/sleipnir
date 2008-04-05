#include "stdafx.h"

/*!
 * \page Randomizer Randomizer
 * 
 * Randomizer removes all values from the given DAT/DAB and randomly inserts the same number (or a requested
 * number) of ones.  This can be used to quickly create a random gold standard with a known number of positive
 * examples.
 * 
 * \section sec_overview Overview
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Randomizer -i <data.dab>
 * \endcode
 * 
 * Removes all non-missing values in \c data.dab and assigns an equivalent number of 1s randomly across all
 * gene pairs.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Randomizer/Randomizer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>DAT/DAB file</td>
 *	<td>Input DAT/DAB file to be randomized.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>0</td>
 *	<td>Integer</td>
 *	<td>If nonzero, desired number of non-missing values in the output file.</td>
 * </tr></table>
 */
