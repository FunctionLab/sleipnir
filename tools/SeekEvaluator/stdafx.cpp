/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"

/*!
 * \page SeekEvaluator SeekEvaluator
 * 
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * SeekEvaluator ...
 * \endcode
 * 
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include SeekEvaluator/SeekEvaluator.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Text file</td>
 *	<td>Tab-delimited text file containing two columns, numerical gene IDs (one-based) and unique gene
 *		names (matching those in the input DAT/DAB files).</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Input directory containing DB files</td>
 * </tr><tr>
 *	<td>-D</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Output directory in which database files will be stored.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>.</td>
 *	<td>Text file</td>
 *	<td>Input file containing list of CDatabaselets to combine</td>
 * </tr></table>
 */
