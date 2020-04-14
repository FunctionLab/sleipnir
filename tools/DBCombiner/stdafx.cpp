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
 * Combines a set of DB files generated from different Sleipnir::CDatabase's
 * into one DB file. 
 *
 * Perhaps for space reason, it is sometime not feasible to generate a Sleipnir::CDatabase
 * covering all datasets on one machine or one partition. Consequently, people
 * generate separate Sleipnir::CDatabase's on different machines
 * first, and then join them into one CDatabase instance with the help of DBCombiner.
 * DBCombiner performs the joining on a per DB-file basis, so users still need to repeat the joining
 * for all DB files in the database.
 *
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * DBCombiner -i <genes.txt> -x <db_list.txt> -d <input_dir> -D <output_dir> [-s]
 * \endcode
 *
 * Combines the DB files listed in the \c db_list.txt into one DB.
 *
 * DBCombiner accepts DB files that are generated from different Sleipnir::CDatabase instances, as long as the same
 * gene map was used. In order for
 * DBCombiner to work, only DB files covering the same genes may be
 * combined. This can be ensured by using only DB files with the same ID in the file name (see some sample lines in \c db_list.txt below).
 * The final joined DB will have datasets listed in the order defined by \c db_list.txt.
 *
 * The \c -s option further splits the combined Sleipnir::CDatabaselet into one gene per \c DB file.
 * This \c -s must be enabled for Seek coexpression integrations. (\ref SeekMiner, \ref SeekServer).
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
 * Note that \c database1, \c database2, \c database3 are three Sleipnir::CDatabase's
 * generated for different datasets.
 *
 * Note how we use the same ID \c 00000004 to ensure that the DB files cover
 * the same genes.
 *
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
 *	<td>Input directory containing \c *.db files</td>
 * </tr><tr>
 *	<td>-D</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Output directory in which database files will be stored.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>Input file containing a list of Sleipnir::CDatabaselet's to combine</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>None</td>
 *	<td>off</td>
 *	<td>If enabled, split the combined Sleipnir::CDatabaselet to one gene per \c DB file</td>
 * </tr></table>
 */
