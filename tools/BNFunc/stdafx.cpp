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
 * \page BNFunc BNFunc
 * 
 * BNFunc produces gene sets (and, optionally, answer files) from functional catalog slims (i.e. lists of
 * individual ontology terms).  These gene sets are essentially lists of all genes annotated under the
 * given terms, and can be used with tools such as \ref Answerer to build more complex functional gold
 * standards.
 * 
 * \section sec_overview Overview
 * 
 * BNFunc provides a quick and easy way to retrieve gene sets from a functional catalog.  It consumes a
 * slim file, which is a list of functional catalog terms; it produces one gene set for each term in the
 * slim, each containing the list of genes annotated at or below the appropriate term.  BNFunc has
 * limited abilities to produce a functional gold standard directly (by marking coannotated gene pairs
 * as related, 1, an non-coannotated gene pairs as unrelated, 0), or these gene sets can be used with
 * \ref Answerer to construct a gold standard.
 * 
 * Suppose we've obtained the functional gold slim from <a href="http://function.princeton.edu/grifn">GRIFn</a>
 * and transformed it into a file \c positive_slim.txt in the appropriate format:
 * \code
 * translation	GO:0043037
 * cytoskeleton organization and biogenesis	GO:0007010
 * transcription from RNA polymerase II promoter	GO:0006366
 * ...
 * boron transport	GO:0046713
 * \endcode
 * 
 * You should first create a directory to hold the results (e.g. \c ./positive_sets/) and download the
 * <a href="http://www.geneontology.org/GO.downloads.ontology.shtml">Gene Ontology structure</a> and
 * <a href="http://www.geneontology.org/GO.current.annotations.shtml">yeast annotation</a> files.  Then run:
 * \code
 * BNFunc -i positive_slim.txt -d ./positive_sets/ -y gene_ontology.obo -g gene_association.sgd
 * \endcode
 * This will produce one gene file per term in the slim, e.g. \c translation containing:
 * \code
 * YOR335C
 * YJR047C
 * YGL105W
 * ...
 * \endcode
 * Down through <tt>boron transport</tt> containing:
 * \code
 * YNL275W
 * \endcode
 * 
 * If you want to go on to create a functional gold standard as in
 * <a href="http://www.ncbi.nlm.nih.gov/pubmed/16420673">Myers et al 2005</a> or
 * <a href="http://www.ncbi.nlm.nih.gov/pubmed/17005538">Huttenhower et al 2006</a>, you'll need negative
 * gene sets (i.e. the "minimally related" gene sets) in addition to positive ones.  You can obtain a
 * negative functional slim from the <a href="http://avis.princeton.edu/mefit/mefit/download">MEFIT download
 * site</a>, and transform it into a tab-delimited file \c negative_slim.txt of the proper form:
 * \code
 * development	GO:0007275
 * nitrogen compound metabolism	GO:0006807
 * catabolism	GO:0009056
 * ...
 * regulation of biological process	GO:0050789
 * \endcode
 * Now create another new directory -c ./negative_sets/ and run:
 * \code
 * BNFunc -i negative_slim.txt -d ./negative_sets/ -y gene_ontology.obo -g gene_association.sgd
 * \endcode
 * 
 * Now you have positive and negative gene sets that you can easily use with \ref Answerer.
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNFunc -i <slim.txt> -d <output_dir> -y gene_ontology.obo -g <gene_association.sgd>
 *		-k ko -K <ORG> -m funcat-2.0_scheme -a <funcat-2.0_data_18052006>
 * \endcode
 * 
 * Saves gene lists for the terms specified in \c slim.txt into the directory \c output_dir.  The slim
 * file must list IDs from one or more of the provided functional catalogs.  Only a subset of these need
 * be used: the Gene Ontology (arguments \c -y and \c -g), the KEGG orthology (arguments \c -k and \c -K,
 * with organism codes SCE, HSA, etc.), or the MIPS funcat (\c -m and \c -a).
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNFunc/BNFunc.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>None</td>
 *	<td>Slim text file</td>
 *	<td>Tab-delimited text file containing two columns with one functional catalog term per line.  The
 *		first column specifies a text description of the term, and the second column gives the ID
 *		of the term (e.g. GO:0006796, ko00361, 01.04.01, etc.)</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory in which gene set files are created.  Large slims can create lots of files; use the
 *		default (current directory .) with caution!</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>None</td>
 *	<td>DAB file</td>
 *	<td>If given, produce a functional gold standard from the given positive (and, optionally, negative)
 *		slim in addition to outputting gene lists.  For details, see \ref Answerer.</td>
 * </tr><tr>
 *	<td>-I</td>
 *	<td>None</td>
 *	<td>Slim text file</td>
 *	<td>If given, use the given slim as negative (minimally related) gene sets when producing a functional
 *		gold standard.  For details, see \ref Answerer.</td>
 * </tr><tr>
 *	<td>-y</td>
 *	<td>None</td>
 *	<td>OBO text file</td>
 *	<td>OBO file containing the structure of the Gene Ontology.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>Annotation text file</td>
 *	<td>Gene Ontology annotation file for the desired organism.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>bp</td>
 *	<td>bp, mf, or cc</td>
 *	<td>Gene Ontology namespace to be used for term ID lookups.</td>
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
 *	<td>-s</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output the first common gene name (if present) rather than the systematic name from the
 *		annotation file.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output each gene's database ID rather than the systematic name from the annotation file.</td>
 * </tr><tr>
 *	<td>-l</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, output every available ID and synonym for each gene (tab delimited, one gene per line).</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>0</td>
 *	<td>Double</td>
 *	<td>Fraction of genes to be reserved for testing.  If nonzero, genes will randomly be selected, not
 *		placed in any output set, and printed to standard output.  These can be saved and used for later
 *		holdout evaluation.</td>
 * </tr><tr>
 *	<td>-q</td>
 *	<td>None</td>
 *	<td>SQL text file</td>
 *	<td>If given, gene sets will also be saved as SQL tables in addition to the text file lists.</td>
 * </tr><tr>
 *	<td>-N</td>
 *	<td>None</td>
 *	<td>Directory</td>
 *	<td>If given, gene sets indicating which genes are functionally unrelated to each slim term will be
 *		produced in the requested directory.  This is rarely directly useful, since it's easier to
 *		produce negative gene sets from a second slim file.</td>
 * </tr><tr>
 *	<td>-L</td>
 *	<td>0.05</td>
 *	<td>Double</td>
 *	<td>If two input terms have a hypergeometric p-value of overlap below this threshhold, genes annotated
 *		to the two terms cannot be considered unrelated.  They will either be related (if coannotated to
 *		some other term) or missing (neither related nor unrelated).  Only applies if \c -o or \c -N is
 *		given.</td>
 * </tr></table>
 */
