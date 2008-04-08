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
 * \page Dat2Dab Dat2Dab
 * 
 * Dat2Dab converts tab-delimited text DAT files into binary DAB files and vice versa.  It can also convert
 * PCL and DAS files (see Sleipnir::CDat), perform a variety of normalizations or filters during the
 * conversion process, or lookup individual genes' or gene pairs' values from DAB files.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Dat2Dab -i <data.dab> -o <data.dat>
 * \endcode
 * 
 * Convert the input binary DAB file \c data.dab into the output tab-delimited text DAT file \c data.dat.
 * 
 * \code
 * Dat2Dab -o <data.dab> -n -f -d < <data.dat>
 * \endcode
 * 
 * Read a text DAT file \c data.dat from standard input, allowing duplicates, normalize all scores to the
 * range [0,1], then invert them and save the results to the binary DAB file \c data.dab.
 * 
 * \code
 * Dat2Dab -i <data.dab> -m -l <gene1> -L <gene2>
 * \endcode
 * 
 * Open the binary DAB file \c data.dab using memory mapping and output the score for the gene pair \c gene1
 * and \c gene2.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Dat2Dab/Dat2Dab.ggo
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
 *	<td>Input DAT, DAB, DAS, or PCL file.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>DAT/DAB file</td>
 *	<td>Output DAT, DAB, or DAS file.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output one minus the input's values.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to the range [0,1] before processing.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, normalize input edges to z-scores (subtract mean, divide by standard deviation) before
 *		processing.</td>
 * </tr><tr>
 *	<td>-r</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, transform input values to integer ranks before processing.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Text gene list</td>
 *	<td>If given, use only gene pairs for which both genes are in the list.  For details, see
 *		Sleipnir::CDat::FilterGenes.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>None</td>
 *	<td>Double</td>
 *	<td>If given, remove all input edges below the given cutoff (after optional normalization).</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, replace all missing values with zeros.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, allow (with a warning) duplicate pairs in text-based input.</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, only print list of genes that would be included in the normal output file.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>None</td>
 *	<td>String</td>
 *	<td>If given, lookup all values for pairs involving the requested gene.</td>
 * </tr><tr>
 *	<td>-L</td>
 *	<td>None</td>
 *	<td>String</td>
 *	<td>If given with \c -l, lookup all values for the requested gene pair.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>None</td>
 *	<td>Gene text file</td>
 *	<td>If given with \c -l, lookup all pairs between \c -l and the given gene set.  If given alone,
 *		lookup all pairs between genes in the given set.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>None</td>
 *	<td>Gene pair text file</td>
 *	<td>Tab-delimited text file containing two columns, both gene IDs.  If given, replace each gene ID
 *		from the first column with the corresponding ID in the second column.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, produce output in a tab-delimited half matrix table.  Not recommended for DAT/DABs with
 *		more than a few dozen genes!</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr></table>
 */
