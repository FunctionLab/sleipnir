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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"

/*!
 * \page Counter Counter
 * 
 * Counter performs three related tasks:
 * <ol>
 *	<li>Given data in the form of DAB/QUANT file pairs and a gold standard answer file, it outputs the raw
 *		counts that would be used to generate conditional probability tables in a Bayesian classifier.  This
 *		is analogous to the learning mode of \ref BNCreator.</li>
 *	<li>Given raw counts and (optionally) a prior weight for each dataset, it outputs the resulting Bayesian
 *		classifiers with properly normalized probability tables (also analogous to \ref BNCreator
 *		learning, although much faster since it starts with pre-counted data).</li>
 *	<li>Given a collection of Bayesian classifiers, it performs inference to output global or context-specific
 *		probabilities of functional relationship for all gene pairs (analogous to the evaluation mode of
 *		\ref BNCreator).</li>
 * </ol>
 * 
 * \section sec_overview Overview
 * 
 * Counter provides fine-grained control over the process of counting values in data, constructing probability
 * tables from these counts, and using these probability tables for Bayesian inference.  Like \ref BNCreator,
 * given one or more biological datasets (stored as Sleipnir::CDat objects in DAT/DAB/etc. files) and a
 * functional gold standard, Counter can construct naive Bayesian classifiers to probabilistically integrate
 * the given data.  However, it does this by providing the intermediate count information describing the number
 * of times each data value is present, and these counts can be weighted (i.e. regularized, see
 * Heckerman/Geiger/Chickering 1995 or Steck/Jaakkola 2002) before they are normalized into probability
 * distributions.  This allows datasets known to provide more diverse or novel information to be upweighted,
 * and it provides a way of doing very rapid Bayesian inference using naive classifiers independent of the
 * SMILE library.
 * 
 * To use Counter, suppose you have a directory containing data files of the form:
 * \code
 * MICROARRAY.dab
 * MICROARRAY.quant
 * CUR_COMPLEX.dab
 * CUR_COMPLEX.quant
 * TF.dab
 * TF.quant
 * ...
 * SYNL_TRAD.dab
 * SYNL_TRAD.quant
 * \endcode
 * Each data file is a Sleipnir::CDat, either a DAT or a DAB, containing experimental results.  Each QUANT
 * file describes how to discretize that data for use with the Bayesian network (the number of bins in the
 * QUANT must be the same as the number of values taken by the corresponding node in the Bayesian network).
 * These files should all be placed in the same directory (e.g. \c ./data/); in a different location, you
 * should assemble a functional gold standard (e.g. \c ANSWERS.dab, possibly constucted by \ref Answerer).
 * 
 * Counter is most useful for context-specific learning and evaluation, since it allows parallelization and
 * storage of counts from many datasets over many biological contexts.  Each context generally represents a
 * pathway, process, or other biological area in which different datasets are expected to behave differently
 * (e.g. a microarray dataset will be informative regarding protein translation, since ribosomes are highly
 * transcriptionally regulated, but it will not measure post-transcriptional regulation such as phosphorylation
 * in MAPK cascades).  Each context is provided to Sleipnir as a single text file containing one gene per line,
 * all collected in the same directory (e.g. \c ./contexts/):
 * \code
 * DNA_catabolism.txt
 * DNA_integration.txt
 * ...
 * mitochondrion_organization_and_biogenesis.txt
 * mitotic_cell_cycle.txt
 * translation.txt
 * \endcode
 * 
 * To generate context-specific data counts from this information, you might create an empty \c ./output/
 * directory and run:
 * \code
 * Counter -w ANSWERS.dab -d ./data/ -m -t 4 -o ./output/ ./contexts/*.txt
 * \endcode
 * where \c -m indicates that the input data files should be memory mapped (generally improving performance)
 * and <tt>-t 4</tt> uses four threads in parallel.  This will generate one file per context in the
 * \c ./output/ directory, each of the form:
 * \code
 * DNA_catabolism	5
 * 79000	190
 * MICROARRAY
 * 1000	4900	9100	5000	930	130	7
 * 3	15	28	18	2	0	0
 * CUR_COMPLEX
 * 79000	2
 * 190	0
 * ...
 * SYNL_TRAD
 * 79000	10
 * 190	15
 * ...
 * \endcode
 * 
 * Here, each file is named for a context (e.g. \c DNA_catabolism.txt) and begins with a header giving its name
 * (e.g. \c DNA_catabolism) and the number of datasets it contains (e.g. 5).  The values below this header give
 * the total counts for unrelated and related pairs in the relevant subset of the answer file (e.g. 79000 and
 * 190 in the context of DNA catabolism).  The appropriate number of dataset blocks follow this, each giving
 * the count of values found in each of that dataset's discretized bins (based on its QUANT file) for the
 * unrelated and related pairs, respectively.  To generate a \c global.txt file for the global context (i.e.
 * for the entire answer file), run Counter with no context arguments on the command line:
 * \code
 * Counter -w ANSWERS.dab -d ./data/ -m -o .
 * \endcode
 * 
 * Now, given a directory with count files for each context, you can create regularized Bayesian classifiers
 * from them, either in human-readable (X)DSL format (for use with SMILE/GeNIe) or in a compact binary format
 * for rapid inference.  To generate (X)DSL files, create an empty \c ./networks/ directory and run:
 * \code
 * Counter -c ./output/ -o ./networks/ -s datasets.txt -b ./global.txt -l
 * \endcode
 * This will generate one (X)DSL file per context in the \c ./networks/ directory (including \c global.xdsl).
 * To instead store these classifiers in a binary format for later Bayesian inferernce, run:
 * \code
 * Counter -c ./output/ -o ./networks.bin -s datasets.txt -b ./global.txt
 * \endcode
 * In these commands, \c datasets.txt is a tab-delimited text file containing two columns, the first a
 * one-based integer index and the second a list of each dataset's name:
 * \code
 * 1	MICROARRAY
 * 2	CUR_COMPLEX
 * ...
 * 5	SYNL_TRAD
 * \endcode
 * This is the same format as is used with other tools such as \ref BNServer.
 * 
 * One of Counter's unique features is the ability to regularize the parameters of the Bayesian classifiers
 * constructed from data counts.  This means that each dataset's probability distributions can be weighted
 * according to a prior trust in that dataset: probability distributions from trusted datasets will be used
 * unchanged, but less trusted datasets will be made closer to a uniform distribution (and thus have less
 * impact on the eventual predictions).  Weights are provided as a combination of a pseudocount parameter,
 * which determines the effective number of counts in each CPT, and a file of alphas, which give the weight
 * (relative to the pseudocounts) to give to a uniform prior for each dataset.  For example, suppose we
 * normalize to a pseudocount total of 100.  Then we might have a tab-delimited text file \c alphas.txt
 * containing:
 * \code
 * MICROARRAY	100
 * CUR_COMPLEX	0
 * ...
 * SYNL_TRAD	10
 * \endcode
 * This means that the probabilities for the \c MICROARRAY node will be equally weighted between the actual
 * counts and a uniform prior, those for the \c CUR_COMPLEX node will be drawn entirely from the data, and
 * those for the \c SYNL_TRAD node will use ~91% information from the data (100/110) and ~9% a uniform prior.
 * To generate Bayesian classifiers reflecting these weights, run:
 * \code
 * Counter -c ./output/ -o ./networks.bin -s datasets.txt -b ./global.txt -p 100 -a ./alphas.txt
 * \endcode
 * For more information, see Sleipnir::CBayesNetMinimal::OpenCounts.
 * 
 * Finally, to use your learned Bayesian classifiers and genomic data to infer functional relationships, you
 * can create an empty \c ./predictions/ directory and run:
 * \code
 * Counter -n ./networks.bin -o ./predictions/ -d ./data/ -s datasets.txt -g genes.txt -m -t 4 ./contexts/*.txt
 * \endcode
 * This will generate one DAB file per context (e.g. \c ./predictions/DNA_catabolism.dab), each containing the
 * probability of functional relationship for each gene pair predicted from the datasets in \c ./data/ and the
 * classifiers in \c ./networks.bin.  The \c genes.txt file is of the same format as \c datasets.txt and lists
 * each gene in the genome, e.g.
 * \code
 * 1	YAL068C
 * 2	YAL066W
 * 3	YAL065C
 * ...
 * \endcode
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Counter -w <answers.dab> -d <data_dir> -o <output_dir> <contexts.txt>*
 * \endcode
 * 
 * For each context \c contexts.txt, generate a counts file in \c output_dir summarizing the number of each
 * data value for the DAT/DAB files in \c data_dir (which must have associated QUANT files) relative to the
 * functional gold standard in \c answers.dab.  A global count file is generated if no contexts are provided
 * on the command line.
 * 
 * \code
 * Counter -c <counts_dir> -o <networks.bin> -s <datasets.txt> -b <global.txt> -p <pseudocounts> -a <alphas.txt>
 * \endcode
 * 
 * Using the counts previously output to \c counts_dir and the global counts file \c global.txt, save a set of
 * Bayesian classifiers in \c networks.bin (one per context plus a global default classifier) each containing
 * one node per dataset as specified in \c datasets.txt.  Probability distributions can be optionally
 * regularized using the effective pseudocount number \c pseudocounts and the relative weight of a uniform
 * prior for each node given in \c alphas.txt.
 * 
 * \code
 * Counter -n <networks.bin> -o <output_dir> -d <data_dir> -s <datasets.txt> -g <genes.txt> <contexts.txt>*
 * \endcode
 * 
 * Performs Bayesian inference for each classifier previously saved in \c networks.bin, producing one predicted
 * functional relationship network DAT/DAB file per \c contexts.txt in the directory \c output_dir, using
 * data from \c data_dir, the dataset list in \c datasets.txt, and the genome list in \c genes.txt.  A global
 * relationship network is produced if no contexts are provided on the command line.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Counter/Counter.ggo
 * 
 * <table><tr>
 * 	<th>Flag</th>
 * 	<th>Default</th>
 * 	<th>Type</th>
 * 	<th>Description</th>
 * </tr><tr>
 * 	<td>None</td>
 * 	<td>None</td>
 * 	<td>Text files</td>
 * 	<td>Contexts used for calculating context-specific counts or producing context-specific functional
 *		relationship predictions.  Each is a text file containing one gene per line.  When no contexts are
 *		provided, global (context-independent) calculations are performed.</td>
 * </tr><tr>
 * 	<td>-w</td>
 * 	<td>None</td>
 * 	<td>DAT/DAB file</td>
 * 	<td>Activates count generation mode.  Functional gold standard for counting.  Should consist of gene pairs
 *		with scores of 0 (unrelated), 1 (related), or missing (NaN).</td>
 * </tr><tr>
 * 	<td>-c</td>
 * 	<td>None</td>
 * 	<td>Directory</td>
 * 	<td>Activates Bayesian classifier generation mode.  Directory containing previously calculated count files,
 *		one text file per context.  Should not contain the global count file.</td>
 * </tr><tr>
 * 	<td>-n</td>
 * 	<td>None</td>
 * 	<td>Binary file</td>
 * 	<td>Activates Bayesian inference (functional relationship prediction) mode.  Binary file containing
 *		previously calculated Bayesian classifiers.</td>
 * </tr><tr>
 * 	<td>-o</td>
 * 	<td>None</td>
 * 	<td>Directory or binary file</td>
 * 	<td>In count generation mode, directory in which data value count files are placed.  In classifier
 *		generation mode, file in which binary classifiers are saved or directory in which (X)DSL files are
 *		saved.  In inference mode, directory in which context-specific functional relationship DAT/DAB files
 *		are created.</td>
 * </tr><tr>
 * 	<td>-d</td>
 * 	<td>.</td>
 * 	<td>Directory</td>
 * 	<td>Directory from which data DAT/DAB files (and accompanying QUANT files) are read.</td>
 * </tr><tr>
 * 	<td>-s</td>
 * 	<td>None</td>
 * 	<td>Text file</td>
 * 	<td>Tab-delimited text file containing two columns, the first a one-based integer index and the second
 *		the name of each dataset to be used (excluding the DAT/DAB suffix, e.g. \c MICROARRAY,
 *		\c CUR_COMPLEX, etc.)</td>
 * </tr><tr>
 * 	<td>-g</td>
 * 	<td>None</td>
 * 	<td>Text file</td>
 * 	<td>Tab-delimited text file containing two columns, the first a one-based integer index and the second
 *		the unique identifier of each gene in the genome.</td>
 * </tr><tr>
 * 	<td>-b</td>
 * 	<td>None</td>
 * 	<td>Text file</td>
 * 	<td>Count file containing default (global) values for global inference or for fallback in contexts with
 *		too little data.</td>
 * </tr><tr>
 *	<td>-Z</td>
 *	<td>None</td>
 *	<td>Tab-delimited text file</td>
 *	<td>If given, argument must be a tab-delimited text file containing two columns, the first node
 *		IDs (see \ref BNCreator) and the second bin numbers (zero indexed).  For each node ID present in
 *		this file, missing values will be substituted with the given bin number.</td>
 * </tr><tr>
 * 	<td>-p</td>
 * 	<td>-1</td>
 * 	<td>Integer</td>
 * 	<td>If not -1, the effective number of pseudocounts to use relative to the weights in the given alphas file
 *		(if any).</td>
 * </tr><tr>
 * 	<td>-a</td>
 * 	<td>None</td>
 * 	<td>Text file</td>
 * 	<td>If given, tab-delimited text file containing dataset IDs and the relative weight given to a uniform
 *		prior for each dataset.</td>
 * </tr><tr>
 * 	<td>-y</td>
 * 	<td>.</td>
 * 	<td>Directory</td>
 * 	<td>Directory in which temporary files are generated during inference mode.</td>
 * </tr><tr>
 * 	<td>-l</td>
 * 	<td>off</td>
 * 	<td>Flag</td>
 * 	<td>If on, SMILE (X)DSL files are created in classifier generation mode rather than a single binary
 *		file.</td>
 * </tr><tr>
 * 	<td>-x</td>
 * 	<td>on</td>
 * 	<td>Flag</td>
 * 	<td>If on, XDSL files are generated instead of DSL files.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>1</td>
 *	<td>Integer</td>
 *	<td>Number of simultaneous threads to use for individual CPT learning.  Threads are per classifier node
 *		(dataset), so the number of threads actually used is the minimum of \c -t and the number of
 *		datasets.</td>
 * </tr></table>
 */
