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
 * \page COALESCE COALESCE
 * 
 * COALESCE performs regulatory modules prediction (expression biclustering and de novo sequence motif
 * analysis) as described in Huttenhower et al. 2009.  The algorithm consumes gene expression data and,
 * if motif discovery is performed, DNA sequence; it outputs zero or more predicted regulatory modules,
 * each consisting of a set of genes, a subset of conditions under which they are coexpressed, and zero
 * or more (under)enriched sequence motifs.  It predicts modules serially by seeding a set of correlated
 * genes, selecting conditions under which they are significantly coregulated, selecting motifs significantly
 * under/over-enriched, and using Bayesian integration to combine these features add highly probable genes
 * (and remove improbable genes).
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * COALESCE -i <input.pcl>
 * \endcode
 * 
 * Perform expression biclustering on the expression data in \c input.pcl (which may contain missing values).
 * 
 * \code
 * COALESCE -i <input.pcl> -f <input.fasta>
 * \endcode
 * 
 * Perform regulatory module prediction (biclustering plus motif discovery) on the expression data in
 * \c input.pcl and using the sequence data in \c input.fasta.  The latter may consist of several subtypes
 * of sequences (e.g. upstream/downstream flanks) denoted by tab-delimited markers in the FASTA IDs, and
 * introns and exons (or flanks and UTRs) may be differentiated by capitalization, e.g.:
 * \code
 * > GENE1	5
 * aaccGGTT
 * > GENE1
 * TGCAtgcaACGT
 * > GENE1	3
 * TTGGccaa
 * > GENE2	5
 * ...
 * \endcode
 * where \c aacc represents \c GENE1 's upstream flank, \c GGTT its 5' UTR, \c TGCA its first exon, \c TTGG its
 * 3' UTR, and so forth.
 * 
 * \code
 * COALESCE -i <input.pcl> -f <input.fasta> -o <output_dir>
 * \endcode
 * 
 * Perform regulatory module prediction on the expression data in \c input.pcl and using the sequence data in
 * \c input.fasta, depositing PCL and motif files progressively in \c output_dir (in addition to the
 * module summaries printed to standard output during normal execution).
 * 
 * 
 * \code
 * COALESCE -i <input.pcl> -j <modules.txt/modules_dir>
 * \endcode
 * 
 * Postprocess a set of preliminary predicted modules stored in standard output format in \c modules.txt (or
 * in separate files in \c modules_dir) that were generated from data in \c input.pcl.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include COALESCE/COALESCE.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>PCL file</td>
 *	<td>Input expression data to be biclustered; may contain missing values.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>None</td>
 *	<td>FASTA file</td>
 *	<td>If given, input sequence data to be mined for regulatory motifs.  Can contain sub-types of
 *		sequence as described above; only gene IDs also present in \c -i will be analyzed.</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>If given, input description of dataset blocks with expected covariance in \c -i.  Each line of the
 *		file is considered to be one dataset, consisting of a tab-delimited list of one or more condition
 *		identifiers from \c -i.  Conditions not listed in the file will be treated as independent
 *		single-condition datasets.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>If given, output directory into which regulatory module PCL and predicted motif files are placed.
 *		Otherwise, module information is printed to standard output.</td>

 * </tr><tr>
 *	<th colspan="4">Algorithm Parameters</th>
 * </tr><tr>
 *	<td>-p</td>
 *	<td>0.95</td>
 *	<td>Double (probability)</td>
 *	<td>Probability threshhold for including genes in a regulatory module.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>0.05</td>
 *	<td>Double (p-value)</td>
 *	<td>P-value threshhold for including conditions in a regulatory module.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>0.05</td>
 *	<td>Double (p-value)</td>
 *	<td>P-value threshhold for including motifs in a regulatory module.</td>
 * </tr><tr>
 *	<td>-C</td>
 *	<td>0.5</td>
 *	<td>Double (z-score)</td>
 *	<td>Z-score threshhold for including conditions in a regulatory module.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>0.5</td>
 *	<td>Double (z-score)</td>
 *	<td>Z-score threshhold for including motifs in a regulatory module.</td>

 * </tr><tr>
 *	<th colspan="4">Sequence Parameters</th>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>7</td>
 *	<td>Integer</td>
 *	<td>Number of base pairs in minimal k-mer motif seeds; longer motifs are built out of units of
 *		this length.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>0.05</td>
 *	<td>Double (p-value)</td>
 *	<td>P-value threshhold for merging motifs with similar distributions among genes; must also meet \c -G.</td>
 * </tr><tr>
 *	<td>-G</td>
 *	<td>2.5</td>
 *	<td>Double</td>
 *	<td>Edit distance threshhold for merging motifs with similar sequences; must also meet -c \g.</td>
 * </tr><tr>
 *	<td>-y</td>
 *	<td>1</td>
 *	<td>Double</td>
 *	<td>Edit distance penalty for gaps when comparing motifs.</td>
 * </tr><tr>
 *	<td>-Y</td>
 *	<td>2.1</td>
 *	<td>Double</td>
 *	<td>Edit distance penalty for mismatches when comparing motifs.</td>

 * </tr><tr>
 *	<th colspan="4">Performance Parameters</th>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>0.05</td>
 *	<td>Double (p-value)</td>
 *	<td>P-value threshhold for determining significantly correlated genes when seeding a new module.</td>
 * </tr><tr>
 *	<td>-N</td>
 *	<td>100000</td>
 *	<td>Integer</td>
 *	<td>Maximum number of gene pairs to sample when selecting genes to seed a new module.</td>
 * </tr><tr>
 *	<td>-q</td>
 *	<td>None</td>
 *	<td>String</td>
 *	<td>If given, sequence subtypes to use when predicting motifs; otherwise, all sequence types are
 *		analyzed.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>5000</td>
 *	<td>Integer (base pairs)</td>
 *	<td>Resolution in base pairs with which motif frequencies are tracked, i.e. in units of 1 / \c -b.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>5</td>
 *	<td>Integer</td>
 *	<td>Minimum number of genes required for a module to be preserved.</td>
 * </tr><tr>
 *	<td>-E</td>
 *	<td>100</td>
 *	<td>Integer</td>
 *	<td>Maximum number of motifs to be merged during module convergence.</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>1000</td>
 *	<td>Integer</td>
 *	<td>Maximum number of motifs to be associated with a module during convergence.</td>

 * </tr><tr>
 *	<th colspan="4">Postprocessing Parameters</th>
 * </tr><tr>
 *	<td>-j</td>
 *	<td>None</td>
 *	<td>Text file/directory</td>
 *	<td>Input text file of modules (as produced on standard output) or directory of module files (as
 *		produced by \c -o) to be postprocessed.</td>
 * </tr><tr>
 *	<td>-K</td>
 *	<td>None</td>
 *	<td>Text file</td>
 *	<td>If given, text file containing known TFs and motifs.  Each line should be tab-delimited text in
 *		which the first column contains a TF identifier (not necessarily unique), and the subsequent 4n
 *		columns contain the PWM values (as floating point fractions) of the n base pairs of the motif,
 *		e.g. <tt>GATA 0 0 1 0 1 0 0 0 0 0 0 1 1 0 0 0</tt>.</td>
 * </tr><tr>
 *	<td>-F</td>
 *	<td>0.05</td>
 *	<td>Double (p-value)</td>
 *	<td></td>
 * </tr>Score threshhold for labeling a predicted motif as a known TF based on -c -K.<tr>
 *	<td>-S</td>
 *	<td>pvalue</td>
 *	<td>pvalue, rmse, or js</td>
 *	<td>Scoring method for labeling a predicted motif as a known TF; options are p-value of PWM
 *		correlation, root-mean-square error between PWMs, or Jenson-Shannon divergence between PWMs.</td>
 * </tr><tr>
 *	<td>-J</td>
 *	<td>1</td>
 *	<td>Double (fraction)</td>
 *	<td>Minimum overlap fraction for two preliminary modules to be merged.</td>
 * </tr><tr>
 *	<td>-L</td>
 *	<td>0.5</td>
 *	<td>Double (fraction)</td>
 *	<td>Minimum fraction of merged preliminary modules in which a gene must be present in order to be
 *		maintained in the resulting postprocessed module.</td>
 * </tr><tr>
 *	<td>-T</td>
 *	<td>1</td>
 *	<td>Double (z-score)</td>
 *	<td>Z-score threshhold of cocluster frequencies of genes in a preliminary module above which they
 *		must occur to be maintained in the resulting postprocessed module.</td>
 * </tr><tr>
 *	<td>-R</td>
 *	<td>On</td>
 *	<td>Flag</td>
 *	<td>If given, convert reverse complement motifs to a single strand before outputting their PWMs.</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>0.3</td>
 *	<td>Double (bits)</td>
 *	<td>Information threshhold (in bits) for a preliminary motif to be preserved in a postprocessed module.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>2500</td>
 *	<td>Integer</td>
 *	<td>Maximum number of motifs to be merged exactly into a postprocessed module; motif counts above this
 *		threshhold are merged heuristically.</td>

 * </tr><tr>
 *	<th colspan="4">Miscellaneous Parameters</th>
 * </tr><tr>
 *	<td>-e</td>
 *	<td>Off</td>
 *	<td>Flag</td>
 *	<td>If given, automatically detect and normalize single-channel conditions by log-transforming them
 *		against the per-gene median.</td>
 * </tr><tr>
 *	<td>-O</td>
 *	<td>On</td>
 *	<td>Flag</td>
 *	<td>If given, generate standard output progressively as modules are finalized.</td>

 * </tr><tr>
 *	<th colspan="4">Standard Parameters</th>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>1</td>
 *	<td>Integer</td>
 *	<td>Number of simultaneous threads to use for applicable stages of module formation.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr></table>
 */
