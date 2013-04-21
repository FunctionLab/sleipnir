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
 * \page SVMperfer SVMperfer
 * 
 * SVMperfer performs SVM learning using the SVMperf library.  It supports cross validation and
 * reading from binary PCL files created by PCL2Bin.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * SVMperfer -l <labels_file> -p <params_file> -i <data.bin> -o <output_directory> -a
 * \endcode
 * 
 * The labels file is of the format (NOTE WELL: IN ALL THE FOLLOWING FORMATS DELIMITERS ARE TABS --
 * doxygen converts them to spaces automatically).
 * \code
 * ACTA2	-1
 * ACTN4	1
 * ADAM10	-1
 * AGRN	1
 * AGTR1	-1
 * ALDOB	-1
 * ALOX12	1
 * ANGPT2	1
 * APOA4	1
 * AQP1	1
 * \endcode
 * where -1 indicates negative and 1 indicates positive.  The examples must be separated with tabs.
 * 
 * Output is of the format
 * \code
 * IGHV1-69	0	1.94073
 * DAG1	1	1.9401
 * FNDC3B	0	1.93543
 * HPGD	-1	1.93181
 * TPSAB1	0	1.92928
 * CLIC5	1	1.92759
 * \endcode
 * where the first column is the example name, the second column is the gold standard status (matching labels)
 * and the third column is the prediction from the SVM.
 * 
 * The params_file is of the format
 * \code
 * 10	0.1	0.5
 * 10	0.01	0.5
 * 10	0.001	0.5
 * 10	0.0001	0.5
 * 10	0.00001	0.5
 * 10	0.000001	0.5
 * \endcode
 * where the first column represents the error function, the second column represents the tradeoff constant
 * and the third column represents k_value (for precision at k recall, but unused for the AUC error function
 * in the example above.
 * 
 * SVMperfer can also be used to output a model or learn a network, although currently those features are undocumented.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include SVMperfer/SVMperfer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>PCL/BIN file</td>
 *	<td>Input PCL file</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>Output directory.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>None</td>
 *	<td>Labels file</td>
 *	<td>The file with examples formatted as noted above.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>None</td>
 *	<td>Model file</td>
 *	<td>If present, output the learned model to this file.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on output predictions for all genes in the PCL.</td>
 * </tr><tr>
 *	<td>-S</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, use slack rescaling.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>int</td>
 *	<td>Number of columns to skip from PCL file.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>Normalize PCL to 0 mean, 1 variance.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>5</td>
 *	<td>int</td>
 *	<td>Number of cross validation intervals.</td>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>10</td>
 *	<td>int</td>
 *	<td>Which loss function should be used? (options: 0, 1, 2, 3, 4, 5, 10).</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>0.5</td>
 *	<td>float</td>
 *	<td>value of k for precision@k or recall@k.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>1</td>
 *	<td>float</td>
 *	<td>SVM tradeoff constant C (note that this differs from the version in SVM light by a constant factor, check SVMPerf docs for details).</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>None</td>
 *	<td>Filename</td>
 *	<td>Parameters file (to test with multiple parameters).</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>Memory map binary input PCLs (BIN files).</td>
 * </tr></table>
 */
