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
 * \page DBCombiner DBCombiner
 * 
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * DBCombiner -i <genes.txt> -x <db_list.txt> -d <input_dir> -D <output_dir> [-s]
 * \endcode
 * Combines multiple \c .db files listed in the \c db_list.txt into one \c .db file.
 * Perhaps for space reason, it is sometime not feasible to generate the CDatabase 
 * covering all datasets on one machine or one partition. Consequently, people 
 * generate separate CDatabases on different machines
 * first, and then join them one .db file at a time with DBCombiner.
 * Thus DBCombiner accepts CDatabase's that are generated separately. In order for
 * DBCombiner to work, only the \c db files covering the same genes may be
 * combined. 
 * The resulted \c .db file will cover datasets in the same order as the order of 
 * DB files in \c db_list.txt.
 * The \c -s option performs a splitting after the combining is done. 
 * It splits the combined CDatabaselet to one gene per \c .db file. 
 * This is required for SeekMiner, SeekServer.
 * 
 * Sample lines from the \c genes.txt file:
 * \code
 * 1    1
 * 2    10
 * 3    100
 * 4    1000
 * 5    10000
 * 6    100008589
 * \endcode
 *
 * Sample lines from the \c db_list.txt file:
 * \code
 * /x/y/database1/00000004.db
 * /x/y/database2/00000004.db
 * /x/y/database3/00000004.db
 * \endcode
 * Note that \c database1, \c database2, \c database3 are three CDatabase's 
 * generated for different datasets.
 * Note how we use the same ID 00000004 to ensure the \c db files cover
 * the same genes.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include DBCombiner/DBCombiner.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>Tab-delimited text file containing two columns, numerical gene IDs (one-based) and unique gene
 *		names (matching those in the input DAT/DAB files).</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Input directory containing DB files</td>
 * </tr><tr>
 *	<td>-D</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Output directory in which database files will be stored.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>Input file containing a list of CDatabaselets to combine</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>None</td>
 *	<td>off</td>
 *	<td>If enabled, split the combined CDatabaselet to one gene per .db file</td>
 * </tr></table>
 */
