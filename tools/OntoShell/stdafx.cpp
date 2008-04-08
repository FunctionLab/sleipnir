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
 * \page OntoShell OntoShell
 * 
 * OntoShell provides a uniform command line interface (mimicking a traditional UNIX shell) to multiple
 * functional catalogs such as GO, KEGG, and MIPS.  Ontology terms are laid out like directories in a file
 * system, with gene annotations serving as files.  The structure of the catalogs can be listed and
 * explored, individual genes' annotations can be catted, and TermFinder queries can be processed either
 * interactively or from the command line.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * OntoShell -o gene_ontology.obo -g <gene_association.sgd> -k ko -K <org>
 *		-m funcat-2.0_scheme -a <funcat-2.0_data_18052006>
 * \endcode
 * 
 * Initiate an interactive command prompt exploring functional catalogs and annotations from the
 * Gene Ontology (\c gene_ontology.obo and \c gene_association.sgd), KEGG (\c ko and the three-letter
 * organism code \c org, e.g. SCE, HSA, etc.), and MIPS (\c funcat_scheme and \c funcat_data).
 * 
 * \code
 * OntoShell -o gene_ontology.obo -g <gene_association.sgd> -k ko -K <org>
 *		-m funcat-2.0_scheme -a <funcat-2.0_data_18052006> -x "<command>"
 * \endcode
 * 
 * Initiate a non-interactive session using the requested functional catalogs; executes the command string
 * \c command (which may contain multiple commands separated by semicolons) and outputs the results to
 * standard output.
 * 
 * \code
 * OntoShell -C <configuration.ini>
 * \endcode
 * 
 * Executes OntoShell using the <a href="http://www.gnu.org/software/gengetopt/gengetopt.html">gengetopt</a>
 * command line arguments stored in the configuration file \c configuration.ini.  This can be extremely
 * useful for storing complex command lines, e.g. the paths to all of the functional catalog files.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include OntoShell/OntoShell.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>None</td>
 *	<td>String</td>
 *	<td>If given, rather than providing the user with an interactive command prompt, execute the given
 *		command (with results going to standard output) and exit.  This is useful for tasks such as saving
 *		the results of a TermFinder query to a file.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, open GO and KEGG using alternate gene IDs (the first synonym rather than the systematic ID)
 *		suitable for processing human (rather than yeast etc.) annotations.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, include annotation database IDs in the gene synonym lists (suitable for WormBase and
 *		FlyBase IDs).</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>OBO text file</td>
 *	<td>OBO file containing the structure of the Gene Ontology.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Annotation text file</td>
 *	<td>Gene Ontology annotation file for the desired organism.</td>
 * </tr><tr>
 *	<td>-k</td>
 *	<td>None</td>
 *	<td>KEGG orthology text file</td>
 *	<td>\c ko file containing the structure and annotations of the KEGG orthology.</td>
 * </tr><tr>
 *	<td>-K</td>
 *	<td>SCE</td>
 *	<td>KEGG organism code</td>
 *	<td>Three letter organism code of annotations to be read from the \c ko file.  Options include SCE for
 *		yeast, HSA for human, DME for fly, CEL for worm, and MMU for mouse.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>None</td>
 *	<td>MIPS schema text file</td>
 *	<td>File containing the schema (structure) of the MIPS functional catalog.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>None</td>
 *	<td>MIPS annotation text file</td>
 *	<td>File containing the annotations to be used with the MIPS functional catalog.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>None</td>
 *	<td>MIPS phenotypes schema text file</td>
 *	<td>File containing the schema (structure) of the MIPS phenotype catalog.</td>
 * </tr><tr>
 *	<td>-A</td>
 *	<td>None</td>
 *	<td>MIPS phenotype annotation text file</td>
 *	<td>File containing the annotations to be used with the MIPS phenotype catalog.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>None</td>
 *	<td>SGD features text file</td>
 *	<td>If given, use gene names and information from the given SGD_features.tab file to annotate genes.</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, perform tab-completion on terms containing zero annotations; otherwise, exclude terms with
 *		no genes from the tab-completion list.</td>
 * </tr></table>
 */
