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
 * reading from binary PCL files created by PCL2Bin.  SVMperfer has been used for Network inference and
 * gene prediction studies.
 *
 * NOTE: Delimiters are tabs in all the following formats -- doxygen converts them to spaces automatically.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * SVMperfer -l <labels_file> -p <params_file> -i <data.bin> -o <output_directory> -a
 * SVMperfer -l <labels_file> -c 5 -t 50 -i <PCL/Dat file> -o <output_file> 
 * \endcode
 * 
 * The label file and Test label file is assumed to have a example name (i.e. row name of input file) and
 * its known label (-1 for negative examples and 1 for positive examples) separated with tabs.  Genes 
 * are examples in the following example.
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
 * 
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
 * where the first column is the example name, the second column is the known label (given in the label 
 * file) and the third column is the SVM prediction (soft value).  Unlabelled examples are given a label
 * of 0.  Examples are sorted by their predicted SVM output soft value.
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
 * 
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
