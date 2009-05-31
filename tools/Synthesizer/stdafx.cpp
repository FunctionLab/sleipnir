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
 * \page Synthesizer Synthesizer
 * 
 * Synthesizer creates synthetic microarray data (PCLs) and gene sequences (FASTAs) based on a simple
 * coexpression model.  It can spike in zero or more transcription factors, which will influence both
 * the expression and sequence composition (by insertion of discrete binding sites) of their synthetic
 * targets.  Synthesizer was created primarily to interact with \ref COALESCE, but it can be used
 * independently, e.g. for testing expression clustering algorithms.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Synthesizer -o <data.pcl> -O <data.fasta> > <description.txt>
 * \endcode
 * 
 * Create synthetic expression data \c data.pcl and gene promoter sequences \c data.fasta, storing a
 * description of the spiked TFs in \c description.txt.
 * 
 * \code
 * Synthesizer -o <data.pcl> -O <data.fasta> -g <genes> -c <conditions>
 *		-n <tfs> -f <sequence.fasta> > <description.txt>
 * \endcode
 * 
 * Create synthetic expression data \c data.pcl with \c conditions conditions and gene promoter sequences
 * \c data.fasta for \c genes genes based on an HMM model of the sequence in \c sequence.fasta, storing a
 * description of the \c tfs spiked TFs in \c description.txt.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Synthesizer/Synthesizer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 * 	<td>-o</td>
 * 	<td>None</td>
 * 	<td>PCL file</td>
 * 	<td>If given, synthetic gene expression data.  Contains the number of genes specified by \c -g, conditions
 *		\c -c, and spiked TF modules \c -n, each influencing gene expression as per \c -M.</td>
 * </tr><tr>
 * 	<td>-O</td>
 * 	<td>None</td>
 * 	<td>FASTA file</td>
 * 	<td>If given, synthetic gene sequence data.  Contains the number of genes specified by \c -g, conditions
 *		\c -c, and spiked TF modules \c -n, each influencing binding sites as per \c -p and \c -P.</td>
 * </tr><tr>
 * 	<td>-g</td>
 * 	<td>5000</td>
 * 	<td>Integer</td>
 * 	<td>Number of synthetic genes to create in PCL and FASTA outputs.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>100</td>
 *	<td>Integer</td>
 *	<td>Number of synthetic expression conditions to create in PCL output.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>10</td>
 *	<td>Integer</td>
 *	<td>Number of synthetic transcription factors to create, influencing expression targets in PCL output
 *		and binding sites in FASTA sequence output.</td>
 * </tr><tr>
 *	<td>-q</td>
 *	<td>0.01</td>
 *	<td>Double (probability)</td>
 *	<td>Probability of a synthetic TF targeting a gene.</td>
 * </tr><tr>
 *	<td>-Q</td>
 *	<td>0.1</td>
 *	<td>Double (probability)</td>
 *	<td>Probability of a synthetic TF being active in an expression condition.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>5</td>
 *	<td>Integer (base pairs)</td>
 *	<td>Minimum length of the randomly generated synthetic motif associated with a TF.</td>
 * </tr><tr>
 *	<td>-T</td>
 *	<td>12</td>
 *	<td>Integer (base pairs)</td>
 *	<td>Maximum length of the randomly generated synthetic motif associated with a TF.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>0</td>
 *	<td>Double</td>
 *	<td>Mean of baseline (i.e. not influenced by any TF) randomly generated expression values.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Standard deviation of baseline (i.e. not influenced by any TF) randomly generated expression
 *		values.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>2</td>
 *	<td>Double</td>
 *	<td>Mean of synthetic TF effect on expression.  Actual expression effects are chosen randomly
 *		as either positive or negative from a normal distribution with this average.</td>
 * </tr><tr>
 *	<td>-S</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Standard deviation of synthetic TF effect on expression.  Actual expression effects are chosen randomly
 *		as either positive or negative from a normal distribution with this standard deviation.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>None</td>
 *	<td>FASTA file</td>
 *	<td>Input sequence file from which an HMM will be built to generate the output synthetic sequences.
 *		Degree of the HMM is specified by \c -d.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>3</td>
 *	<td>Integer</td>
 *	<td>Degree of the HMM used to generate synthetic output sequences; built from the input FASTA file \c -f.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>1000</td>
 *	<td>Integer</td>
 *	<td>Minimum length per gene of randomly generated output sequences.</td>
 * </tr><tr>
 *	<td>-L</td>
 *	<td>3000</td>
 *	<td>Integer</td>
 *	<td>Maximum length per gene of randomly generated output sequences.</td>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>1</td>
 *	<td>Integer</td>
 *	<td>Minimum number of synthetic TFBSs present in a gene's sequence once it has been determined to be a
 *		target of a synthetic TF.</td>
 * </tr><tr>
 *	<td>-P</td>
 *	<td>5</td>
 *	<td>Integer</td>
 *	<td>Maximum number of synthetic TFBSs present in a gene's sequence once it has been determined to be a
 *		target of a synthetic TF.</td>
 * </tr><tr>
 *	<td>-y</td>
 *	<td>None</td>
 *	<td>String</td>
 *	<td>If given, comma-separated list of output sequence types that should contain spiked TFBSs.  \c 5 is a
 *		common value to include spiked TFBSs only in upstream flank sequences.</td>
 * </tr><tr>
 *	<td>-w</td>
 *	<td>60</td>
 *	<td>Integer</td>
 *	<td>Wrap width of generated FASTA files.</td>
 * </tr></table>
 */
