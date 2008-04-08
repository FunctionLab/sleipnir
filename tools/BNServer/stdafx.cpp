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
 * \page BNServer BNServer
 * 
 * BNServer is a complex tool provided a multithreaded TCP/IP interface to real-time Bayesian data
 * integration and inference.  A running BNServer can service client requests over a network for
 * values from specific biological datasets, predicted functional relationships, queries into a functional
 * relationship network, and graph visualization using <a href="http://www.graphviz.org/">Graphviz</a>.
 * 
 * \section sec_overview Overview
 * 
 * As the name implies, BNServer is a network server which can service client requests for information
 * using a simple, binary TCP/IP protocol based on Sleipnir::CServer.  The server loads a variety of
 * information at startup, including one or more biological datasets (stored in a Sleipnir::CDatabase
 * rather than standard Sleipnir::CDat files) and one or more naive classifiers, and can provide various
 * related pieces of information to a client:
 * - Given two genes, retrieve all data values for the pair from all datasets in the database.
 * - Given a single gene, perform inference and return probabilities of functional relationship for all
 *	gene pairs involving that gene.
 * - Given a single gene, indicate which contexts it's most functionally related to.
 * - Given a set of genes, perform inference for all of them and return the portion of the resulting
 *	functional relationship network most related to the query set.
 * - Given a set of genes, perform a <a href="http://search.cpan.org/dist/GO-TermFinder/">TermFinder</a>
 *	functional enrichment query.
 * 
 * Since BNServer is a very complex program, let's first go over the pieces of data necessary to start
 * it running.  First, convert a collection of biological datasets (usually DAB/QUANT file pairs) into a
 * Sleipnir::CDatabase directory.  Let's call that directory \c ./db/, which will contain a bunch of
 * files:
 * \code
 * 00000000.db
 * 00000001.db
 * 00000002.db
 * ...
 * \endcode
 * You'll also need a tab-delimited text file listing each biological context of interest (usually a
 * functional catalog term) in three columns, an integer index (one-based), a textual description, and an ID.
 * Let's call this file \c contexts.txt:
 * \code
 * 1	activation of NF-kappaB transcription factor	GO:0051092
 * 2	adult locomotory behavior	GO:0008344
 * 3	aging	GO:0007568
 * ...
 * \endcode
 * Next, create a tab-delimited text file with two columns, a one-based integer index and a gene ID.  This
 * file (\c genes.txt) will contain every gene of interest for your organism, e.g.:
 * \code
 * 1	YPL149W
 * 2	YHR171W
 * 3	YBR217W
 * \endcode
 * Now tie these two files together by creating a relational mapping between the two lists.  You'll need a
 * tab-delimited \c contexts_genes.txt file with two columns.  The first column is a context ID, the second
 * a gene ID, and there should be one row for each gene annotation to a context term.  For example, if
 * context #1 contains genes #3, 9, and 16, and context #2 contains genes #3 and 4, your \c contexts_genes.txt
 * file might start:
 * \code
 * 1	3
 * 1	9
 * 1	16
 * 2	3
 * 2	4
 * ...
 * \endcode
 * 
 * If you have fully predicted, context-specific functional relationship networks produced by
 * \ref BNUnraveler, you can use \ref Hubber to generate background connectivities for each gene from
 * them.  Let's call this binary file of gene "hubbiness" \c backgrounds.bin.  If you don't have such
 * a file, don't worry; it's optional.
 * 
 * If you have knowledge of disease-associated genes, you can also create a tab-delimited \c diseases_genes.txt
 * file providing relational mappings between disease and gene IDs.  The diseases don't need names, so
 * there's no separate diseases file.  Like the contexts/genes mapping, \c diseases_genes.txt should be a
 * tab-delimited text file; unlike the other file, this one's optional.  If you have it, it should look like:
 * \code
 * 1	18719
 * 2	7473
 * 2	19634
 * 3	4117
 * ...
 * \endcode
 * 
 * Now, for each context you listed in your contexts file, you should have a context-specific Bayes net
 * stored in an (X)DSL file.  This must be a naive classifier with parameters that have already been
 * learned; \ref BNWeaver is ideal for this.  Put all of these files in a single directory,
 * e.g. \c ./contexts/, named the same as their textual description in \c contexts.txt with non-alphanumeric
 * characters substituted with underscores.  For example, \c ./contexts/ might contain:
 * \code
 * activation_of_NF_kappaB_transcription_factor.xdsl
 * adult_locomotory_behavior.xdsl
 * aging.xdsl
 * ...
 * \endcode
 * Note that \e every non-alphanumeric character must be replaced with an underscore.  You should also have
 * one more non-context-specific, global Bayesian network containing the "default" probabilities to be
 * used outside of a specific context.  This can be generated with \ref BNWeaver or, better yet,
 * \ref BNCreator.  Call that file \c default.xdsl.
 * 
 * We're on the home stretch here.  Collect Gene Ontology
 * <a href="http://www.geneontology.org/GO.downloads.ontology.shtml">structure</a> and
 * <a href="http://www.geneontology.org/GO.current.annotations.shtml">annotation</a> files and the KEGG
 * orthology <a href="ftp://ftp.genome.jp/pub/kegg/genes/">ko</a> file.  Make sure
 * <a href="http://www.graphviz.org/">Graphviz</a> is installed and in your current path, and you can
 * finally run:
 * \code
 * BNServer -d ./db/ -i contexts.txt -c contexts_genes.txt -a backgrounds.bin -s diseases_genes.txt
 *		-n ./contexts/ -b default.xdsl -g gene_ontology.obo -G gene_association.sgd -k ko -K SCE
 * \endcode
 * 
 * This will start a server running on the default port, using the default path to Graphviz and generating
 * output files in the current directory (all of which can be modified using other options; see below).
 * What can you do with such a thing?  BNServer will listen on the specified port and service incoming
 * network requests.  As detailed in Sleipnir::CServer, all Sleipnir network communication is prefixed by
 * a byte count for the incoming message.  BNServer accepts several different message types; each should
 * consist of the byte count for the whole message followed by the opcode (a single byte) and appropriate
 * one- or four-byte arguments:
 * <table>
 *	<tr>
 *		<th>Opcode</th>
 *		<th>Name</th>
 *		<th>Arguments</th>
 *		<th>Description</th>
 *	</tr><tr>
 *		<td>0</td>
 *		<td>Inference</td>
 *		<td>Four-byte context ID, zero or more four-byte gene IDs</td>
 *		<td>For each input gene, returns one four-byte floating point value per gene in the genome, each
 *			representing the probability of functional relationship with the input gene in the given
 *			context.  Uninferrable pairs are marked with NaNs.  For example, in a three-gene genome, an input
 *			request of <tt>5 0 2</tt> would return probabilities for the second gene (ID #2) in the default
 *			context (which has no ID, hence #0): <tt>24 0.1 NaN 0.9</tt></td>
 *	</tr><tr>
 *		<td>1</td>
 *		<td>Data</td>
 *		<td>Two four-byte gene IDs</td>
 *		<td>Returns discretized data values for the given gene pair across each dataset in the database.
 *			For example, suppose the database contained six datasets.  Requesting data for genes #1 and #2,
 *			<tt>9 1 1 2</tt>, would return something of the form <tt>6 0 1 0 2 3 1</tt>, with each dataset's
 *			value encoded in a one-indexed byte (and zero representing a missing value).</td>
 *	</tr><tr>
 *		<td>2</td>
 *		<td>Graph</td>
 *		<td>One-byte boolean, four-byte context ID, four-byte neighbor count, zero or more four-byte
 *			gene IDs</td>
 *		<td>For the given gene set, perform Bayesian inference in the given context and retrieve the
 *			requested number of neighbors most related to the query set in the resulting functional
 *			relationship graph.  If the given boolean is true, the resulting graph will be saved in DOT
 *			format in the server's file directory and the filename returned over the network.  If it's
 *			false, the contents of the DOT themselves will be sent back to the caller.</td>
 *	</tr><tr>
 *		<td>3</td>
 *		<td>Contexts</td>
 *		<td>Zero or more four-byte gene IDs</td>
 *		<td>For each gene in the request, return two four-byte floating point values per context indicating
 *			the gene's in-connectivity and background-connectivity in each context.  A gene's in-connectivity
 *			is its average probability of functional relationship with a gene in the context; its
 *			background-connectivity is its average probability of functional relationship with any gene.
 *			For example, for a server with three contexts, a query of <tt>5 3 2</tt> would retrieve gene #2's
 *			in- and background-connectivities with each context: <tt>24 0.9 0.1 0.15 0.12 0.3 0.35</tt>.
 *	</tr><tr>
 *		<td>4</td>
 *		<td>TermFinder</td>
 *		<td>Four-byte ontology ID, four-byte floating point p-value, zero or more four-byte gene IDs</td>
 *		<td>Return the given gene set's functional enrichments in the given ontology below the given p-value
 *			cutoff.  Ontology IDs are 0 for GO BP, 1 for GO MF, 2 for GO CC, and 3 for KEGG.  Results are
 *			prepended with the number of terms followed by, for each term, the null-terminated ID string,
 *			null-terminated description string, four-byte floating point p-value, four-byte integers in-term
 *			hits, term size, query size, background size, number of genes annotated to term, IDs of genes
 *			annotated to term.  Thus, a query of the form <tt>21 4 0 0.05 1 2 3</tt> might result in two
 *			enriched terms: <tt>149 2 GO:0006412 translation 0.01 3 773 3 7455 5 1 2 3 4 5 GO:0006081
 *			aldehyde metabolic process 0.02 2 26 3 7455 4 1 2 8 9</tt></td>
 *	</tr><tr>
 *		<td>5</td>
 *		<td>Diseases</td>
 *		<td>Four-byte context ID, zero or more four-byte gene IDs</td>
 *		<td>For each gene in the request, return two four-byte floating point values per disease indicating
 *			the gene's in-connectivity and background-connectivity to each disease in the requested context.
 *			A gene's in-connectivity is its average probability of functional relationship with a gene in the
 *			disease; its background-connectivity is its average probability of functional relationship with
 *			any gene.  For example, for a server with three diseases, a query of <tt>9 0 2</tt> would
 *			retrieve gene #2's in- and background-connectivities with each gene in the default context, e.g.
 *			<tt>24 0.9 0.1 0.15 0.12 0.3 0.35</tt>.
 *	</tr>
 * </table>
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * BNServer -d <database_dir> -i <contexts.txt> -c <contexts_genes.txt> -n <contexts_dir> -b <default.xdsl>
 *		-g gene_ontology.obo -G <gene_association.sgd> -k ko -K <ORG>
 * \endcode
 * 
 * Starts a server on the default port using data from a Sleipnir::CDatabase in \c database_dir, a context
 * list from \c contexts.txt, context/gene mapping defitions from \c contexts_genes.txt, context-specific
 * Bayesian classifiers from \c contexts_dir, a default global classifier \c default.xdsl, the given
 * Gene Ontology structure and annotation files, and the KEGG orthology \c ko with the given organism code.
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include BNServer/BNServer.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory from which Sleipnir::CDatabase data files are read.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Context text file</td>
 *	<td>Tab-delimited text file from which context indices, names, and string IDs are read.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>None</td>
 *	<td>Mapping text file</td>
 *	<td>Tab-delimited text file from which context indices and associated gene indices are read.</td>
 * </tr><tr>
 *	<td>-a</td>
 *	<td>None</td>
 *	<td>Binary file</td>
 *	<td>Binary file from which each gene's background connectivity for each context is read.  See
 *		\ref Hubber for details.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>None</td>
 *	<td>Mapping text file</td>
 *	<td>Tab-delimited text file from which disease indices and associated gene indices are read.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory from which context-specific Bayesian classifiers ((X)DSL files) are read.</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>Bayesian classifier for the default (global) context.</td>
 * </tr><tr>
 *	<td>-x</td>
 *	<td>on</td>
 *	<td>Flag</td>
 *	<td>If on, assume XDSL files will be used instead of DSL files.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, read binary stored Bayesian classifiers from a file specified by \c -n; if off, assume
 *		\c -n specifies a directory of (X)DSL files.</td>
 * </tr><tr>
 *	<td>-M</td>
 *	<td>None</td>
 *	<td>Binary file</td>
 *	<td>If given, store Bayesian classifiers in a custom binary format in the given filename.  Cannot be
 *		used with \c -m on.  Loading the classifiers from a binary file can be faster than loading several
 *		hundred separate (X)DSLs.  If your classifiers aren't changing, you can load them from (X)DSL
 *		files once, leaving \c -m off and saving them with \c -M, then on subsequent runs turn \c -M off and
 *		load the binary file with \c -m on.</td>
 * </tr><tr>
 *	<td>-g</td>
 *	<td>None</td>
 *	<td>OBO text file</td>
 *	<td>OBO file containing the structure of the Gene Ontology.</td>
 * </tr><tr>
 *	<td>-G</td>
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
 *	<td>-p</td>
 *	<td>1234</td>
 *	<td>Integer</td>
 *	<td>TCP/IP port on which the server should listen.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>100</td>
 *	<td>Integer</td>
 *	<td>Millisecond timeout between server listen polls.  The default should always be fine.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Directory into which the server will place all generated files (e.g. DOTs from graph queries).</td>
 * </tr><tr>
 *	<td>-z</td>
 *	<td>fdp</td>
 *	<td>Executable file</td>
 *	<td>Path to the Graphviz executable used to render DOT files for graph queries.  If \c fdp is in your
 *		path, the default should work fine; otherwise, you should provide an absolute path, e.g.
 *		\c /usr/bin/fdp.  Either \c neato or \c fdp should work interchangeably.</td>
 * </tr></table>
 */
