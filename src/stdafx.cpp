#include "stdafx.h"

/*! \mainpage The Sleipnir Library for Computational Functional Genomics
 * 
 * Greetings, and thanks for your interest in the Sleipnir library!  Sleipnir is a C++ library enabling
 * efficient analysis, integration, mining, and machine learning over genomic data.  This includes a
 * particular focus on microarrays, since they make up the bulk of available data for many organisms,
 * but Sleipnir can also integrate a wide variety of other data types, from pairwise physical interactions
 * to sequence similarity or shared transcription factor binding sites.  All analysis is done with attention
 * to speed and memory usage, enabling the integration of hundreds of datasets covering tens of thousands
 * of genes.  In addition to the core library, Sleipnir comes with a variety of pre-made tools, providing
 * solutions to common data processing tasks and examples to help you use Sleipnir in your own programs.
 * Sleipnir is free, open source, fully documented, and ready to be used by itself or as a component in your
 * computational biology analyses.
 * 
 * \section sec_download Download
 * 
 * - <a href="sleipnir-current.tar.gz">sleipnir-current.tar.gz</a>, the current source code.
 * - <a href="sleipnir-doc-current.tar.gz">sleipnir-doc-current.tar.gz</a>, the current documentation.
 * 
 * Sleipnir and its associated tools are provided as source code that can be compiled under Linux (using
 * gcc), Windows (using Visual Studio or cygwin), or MacOS (using gcc).  For more information, see
 * \ref sec_building "Building Sleipnir".
 * 
 * \section sec_citation Citation
 * 
 * If you use Sleipnir, please cite our publication:
 * 
 * <b>Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya "The Sleipnir library for computational
 * functional genomics" ***, 2008</b>
 * <a href="http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=***">PMID ***</a>
 * 
 * \section sec_building Building Sleipnir
 * 
 * ***
 * 
 * \section sec_uses Example Uses
 * 
 * Sleipnir can be used to satisfy a variety of needs in bioinformatic data processing, from simple
 * data normalization to complex integration and machine learning.  The tools provided with Sleipnir can
 * be used by themselves, or you can integrate the Sleipnir library into your own tools.
 * 
 * \subsection ssec_uses_tools Tools
 * 
 * \subsubsection sssec_uses_tools_microarray Microarray Processing
 * 
 * You're investigating four different knockout strains of yeast.  To assay their transcriptional response
 * to nutrient limitation, you've grown the four cultures on media containing nothing but cheetos for two
 * days, resulting in four two-color microarray time courses.  Rather than using a pooled reference, you've
 * used the zero time point of each time course as its reference.  This leaves you with four PCL datasets,
 * each containing twelve conditions, and each using a different reference.  Your microarray technique is
 * good but not great, so there are some missing values, and the different reference channels make it
 * difficult to compare the different datasets.  What can you do?
 * 
 * - Let's assume you have the datasets in four PCL files named for the gene knockouts: \c pza1.pcl,
 *	\c ber1.pcl, \c rmn1.pcl, and \c cke1.pcl.  First, use \ref answerer KNNImputer to impute and remove
 *	missing values for each file:
 * \code
 * KNNImputer -i pza1.pcl -o pza1_imputed.pcl
 * \endcode
 * - Now, you could concatenate all four PCLs into one using \ref combiner Combiner :
 * \code
 * Combiner -o combined_imputed.pcl *_imputed.pcl
 * \endcode
 * - You could hierarchically cluster this combined file (or the individual files) to view with
 *	<a href="http://jtreeview.sourceforge.net/">Java TreeView</a> using \ref mcluster MCluster :
 * \code
 * MCluster -o combined_imputed.gtr -i combined_imputed.pcl > combined_imputed.cdt
 * \endcode
 * - If you were still concerned about the different references, you could cluster the microarrays in
 *	normalized correlation space instead.  First, generate normalized correlations for each file
 *	individually using \ref distancer Distancer :
 * \code
 * Distancer -i pza1_imputed.pcl -o pza1_imputed.dab
 * \endcode
 * - Next, combine the four resulting DAB files using \ref combiner Combiner, this time averaging pairwise
 *	correlation scores rather than concatenating PCL files:
 * \code
 * Combiner -t dat -o combined_imputed_normalized.dab *_imputed.dab
 * \endcode
 * - Finally, cluster the expression data using these normalized scores rather than the non-normalized
 *	correlations \ref mcluster MCluster used in the last example:
 * \code
 * MCluster -o combined_imputed_normalized.gtr -i combined_imputed_normalized.dab
 *		< combined_imputed.pcl > combined_imputed_normalized.cdt
 * \endcode
 * 
 * \subsubsection sssec_uses_tools_aneuploidy Clustering With Aneuploidies
 * 
 * In your previous microarray experiment, you discover that your ber1 knockout strain developed an
 * aneuploidy halfway through your time course.  The end of the right arm of chromosome one was duplicated
 * in the last six conditions, artificially doubling the expression level of all of its genes.  How can you
 * keep this huge upregulation from driving your clustering?
 * 
 * First, create a PCL file of weights for every gene in every condition.  Let's assume your original
 * \c ber1_imputed.pcl file looks like this:
 * \code
 * ORF	NAME	GWEIGHT	TIME1	TIME2	...	TIME12
 * EWEIGHT			1	1	...	1
 * YAL001C	TFC3	1	0.1	0.2	...	0.12
 * YAL002W	VPS8	1	-0.1	-0.2	...	-0.12
 * ...
 * YAR070C	YAR070C	1	1.1	1.2	...	1.12
 * YAR071W	PHO11	1	2.1	2.2	...	2.12
 * YAR073W	IMD1	1	-1.1	-1.2	...	-1.12
 * YAR075W	YAR075W	1	-2.12	-2.11	...	-2.1
 * ...
 * YPR203W	YPR203W	1	0.12	0.11	...	0.1
 * YPR204W	YPR204W	1	-0.12	-0.11	...	-0.1
 * \endcode
 * The four YAR genes listed here have been duplicated, and their expression levels are correspondingly
 * high.  Create a \b weights PCL file with exactly the same structure, save that the expression values
 * are all replaced by the desired weights of each gene in each condition.  A weight of 1.0 means that the
 * gene should be counted normally, a weight of 0.5 means that it should contribute half as much weight,
 * 2.0 twice as much, and so forth:
 * \code
 * ORF	NAME	GWEIGHT	TIME1	TIME2	...	TIME12
 * EWEIGHT			1	1	...	1
 * YAL001C	TFC3	1	1.0	1.0	...	1.0
 * YAL002W	VPS8	1	1.0	1.0	...	1.0
 * ...
 * YAR070C	YAR070C	1	1.0	1.0	...	0.5
 * YAR071W	PHO11	1	1.0	1.0	...	0.5
 * YAR073W	IMD1	1	1.0	1.0	...	0.5
 * YAR075W	YAR075W	1	1.0	1.0	...	0.5
 * ...
 * YPR203W	YPR203W	1	1.0	1.0	...	1.0
 * YPR204W	YPR204W	1	1.0	1.0	...	1.0
 * \endcode
 * Each of the four duplicated YAR genes should be assigned a weight of 0.5 in the conditions where it was
 * duplicated; thus, the whole row for PHO11 should be:
 * \code
 * YAR071W	PHO11	1	1.0	1.0	1.0	1.0	1.0	1.0	0.5	0.5	0.5	0.5	0.5	0.5
 * \endcode
 * Let's name this file \c ber1_weights.pcl.  Now, run \ref mcluster MCluster with the expression file and
 * the weights file:
 * \code
 * MCluster -o ber1_weighted.gtr -w ber1_weights.pcl -i ber1_imputed.pcl > ber1_weighted.cdt
 * \endcode
 * The resulting cluster output will still contain the doubled expression values, so you can see what the
 * genes' actual expression levels were, but they won't contribute abnormally much to the clustering.
 * 
 * \subsubsection sssec_uses_tools_catalogs Exploring Functional Catalogs
 * 
 * Suppose you've just downloaded the latest and greatest versions of the
 * <a href="http://www.geneontology.org/">Gene Ontology</a>,
 * <a href="http://mips.gsf.de/projects/funcat">MIPS Funcat</a>, and
 * <a href="http://www.genome.jp/kegg/">KEGG Orthology</a>.  You're still chasing down information on your
 * four knockout yeast strains, so you also get the
 * <a href="http://www.geneontology.org/GO.current.annotations.shtml">GO yeast annotations</a> and
 * <a href="ftp://ftpmips.gsf.de/catalogue/annotation_data">Funcat yeast annotations</a>.  This should
 * give you five files:
 * - \c gene_ontology.obo and \c gene_association.sgd for GO.
 * - \c funcat-2.0_scheme and \c funcat-2.0_data_18052006 (or something similar) for MIPS.
 * - \c ko for KEGG (the orthology file can be hard to find; it's on their FTP site).
 * 
 * Let's load them into \ref ontoshell OntoShell and look around:
 * \code
 * OntoShell -o gene_ontology.obo -g gene_assication.sgd -m funcat-2.0_scheme -a funcat-2.0_data_18052006
 *		-k ko -K SCE
 * \endcode
 * This should produce a command line from which you can explore the three ontologies simultaneously:
 * \code
/> ls
- ROOT
O KEGG  1517
O GOBP  6462
O GOMF  6310
O GOCC  6434
O MIPS  6773
O MIPSP 0
/> cat PHO11
YAR071W (PHO11)
One of three repressible acid phosphatases, a glycoprotein that is transported t
o the cell surface by the secretory pathway
KEGG: ko00361            Metabolism; Xenobiotics Biodegradation and Metab...
      ko00740            Metabolism; Metabolism of Cofactors and Vitamins...
GOBP: GO:0006796         phosphate metabolic process
GOCC: GO:0005576         extracellular region
GOMF: GO:0003993         acid phosphatase activity
MIPS: 01.04.01           phosphate utilization
      01.05.01           C-compound and carbohydrate utilization
      01.07              metabolism of vitamins, cofactors, and prostheti...
/> ls -g GOBP/GO:0007624
- GO:0007624         1     0     ultradian rhythm
P GO:0048511         0     1     rhythmic process
 YGL181W(GTS1,FHT1,LSR1)
 * \endcode
 * For more information on specific \ref ontoshell OntoShell commands and capabilities, please see its
 * documentation.
 * 
 * Suppose you've discovered four genes showing unusual activity during your cheeto time courses.  Create a
 * gene list text file for those four genes:
 * \code
 * YAR014C
 * YNL161W
 * YKL189W
 * YOR353C
 * \endcode
 * Suppose this is named \c cheeto_genes.txt.  We can test for functional enrichment among this gene set
 * across all three catalogs in \ref ontoshell OntoShell:
 * \code
/> find -g -l cheeto_genes.txt 0.01
KEGG:
ko04150            0.00791035    1    1    12   1517 Environmental Information Processing; Signal Transduction; ...
GOBP:
GO:0000903         6.55773e-011  4    4    8    6462 cellular morphogenesis during vegetative growth
GO:0016049         4.14193e-006  4    4    103  6462 cell growth
GO:0008361         1.13194e-005  4    4    132  6462 regulation of cell size
...
GOCC:
GO:0030427         6.15338e-006  4    4    153  6434 site of polarized growth
GO:0005933         7.37199e-006  4    4    160  6434 cellular bud
GO:0043332         1.07503e-005  3    4    34   6434 mating projection tip
...
MIPS:
40.01              4.53959e-005  4    4    239  6773 cell growth / morphogenesis
40                 7.86743e-005  4    4    274  6773 CELL FATE
40.01.03           0.00519136    2    4    37   6773 directional cell growth (morphogenesis)
 * \endcode
 * So it looks like eating nothing but cheetos has something to do with vegetative growth!  You could run
 * this same command directly from the command line to save the output in a file for later reference:
 * \code
 * OntoShell -o gene_ontology.obo -g gene_assication.sgd -m funcat-2.0_scheme -a funcat-2.0_data_18052006
 *		-k ko -K SCE -x 'find -g -l cheeto_genes.txt 0.01' > cheeto_genes_enriched_terms.txt
 * \endcode
 * 
 * \subsubsection sssec_uses_tools_bayesian Bayesian Data Integration
 * 
 * You've done about as much by-hand analysis of your cheeto time courses as you can, so you're ready to
 * throw some machine learning algorithms at them.  Suppose you want to construct a predicted functional
 * relationship network specific to your four datasets and the process of "cellular morphogenesis during
 * vegetative growth".
 * - First, you need to assemble a gold standard listing all (or at least some approximation of) known
 *	functional relationships in yeast.  The
 *	<a href="http://avis.princeton.edu/GRIFn/data/GO_functional_slim.txt">GRIFn functional slim</a> lists
 *	several hundred GO terms known to be informative; that is, genes coannotated to these terms are likely
 *	to be functionally related.  Let's turn this functional slim into a gold standard answer file.  Using a
 *	text editor or Excel, flip the file's columns around and trim out extraneous data so that it has two
 *	tab-delimited columns, first the GO term name and then its ID, one per lime:
 * \code
 * translation	GO:0043037
 * cytoskeleton organization and biogenesis	GO:0007010
 * transcription from RNA polymerase II promoter	GO:0006366
 * ...
 * boron transport	GO:0046713
 * \endcode
 * - Now, we'll turn these GO term IDs into gene sets by dumping the genes annotated to each term.  Create
 *	a directory named \c positives and run \ref bnfunc BNFunc:
 * \code
 * BNFunc -o gene_ontology.obo -a gene_assocation.sgd -i GO_functional_slim.txt -d positives
 * \endcode
 * - Using \ref answerer Answerer, we can turn these gene sets into a gold standard DAB file:
 * \code
 * Answerer -p positives -o answers.dab
 * \endcode
 * - You can automatically create and learn a Bayesian network that will integrate your four datasets in a
 *	context-specific manner.  Assuming you've generated the four files \c pza1_imputed.dab and so forth in the
 *	example above, and let's add four quantization QUANT files for them and one for the answer file.  The
 *	answer file's easy; it just contains 0s and 1s, so create a text file \c answers.quant containing one line:
 * \code
 * 0.5	1.5
 * \endcode
 *	Since your four data files are all normalized similarity scores generated from microarrays, you can
 *	discretize them using evenly spaced bins around zero.  Make four files \c pza1_imputed.quant and so forth,
 *	each containing the single tab-delimited line:
 * \code
 * -1.5    -0.5    0.5     1.5     2.5     3.5     4.5
 * \endcode
 * - Once the five QUANT files are created, you can automatically learn a context-specific Bayesian network
 *	using \ref bncreator BNCreator.  Still got the \c cheeto_genes.txt file you created earlier?  Assuming all
 *	of your files are together in the current directory, run:
 * \code
 * BNCreator -w answers.dab -o cheeto_network.xdsl -d . -c cheeto_genes.txt
 * \endcode
 * - Congratulations!  You've learned a Bayesian network by mining your four datasets, looking specifically
 *	for functional information pertaining to the genes of interest in your study.  Let's turn it into a
 *	predicted functional interaction network by running \ref bncreator BNCreator one more time:
 * \code
 * BNCreator -i cheeto_network.xdsl -o cheeto_network.dab -d .
 * \endcode
 * - The pairwise scores now stored in \c cheeto_network.dab each represent a probability of functional
 *	interaction between each gene pair.  You can do all sorts of interesting analyses on this network,
 *	including visualizing portions of it using the <a href="http://function.princeton.edu/pixie">bioPIXIE</a>
 *	algorithm.  If you want to see what the portion of the network around your original genes of interest
 *	looks like, use \ref dat2graph Dat2Graph:
 * \code
 * Dat2Graph -i cheeto_network.dab -q cheeto_genes.txt -k 5 > cheeto_genes_subnetwork.dot
 * \endcode
 * - This DOT file can be converted into a picture using the <a href="http://www.graphviz.org/">Graphviz</a>
 *	tools from AT&T.  It should look something like this:
 *	\image html cheeto_genes_subnetwork.png
 * - If you want to do follow-up analysis on the predicted interaction network, you can either use the
 *	Sleipnir library to build your own tools (see the \ref ssec_uses_library "Core Library" section below),
 *	or you can dump the network as text using \ref dat2dab Dat2Dab:
 * \code
 * Dat2Dab -i cheeto_network.dab -o cheeto_network.dat
 * \endcode
 *	The resulting DAT file is plain text which looks something like this:
 * \code
 * YAL001C YAL040C 0.0155878
 * YAL001C YAL041W 0.242001
 * YAL001C YAL056W 0.345961
 * ...
 * \endcode
 * - Or you might like to cluster your original microarray data using the predicted functional relationship
 *	probabilities as a similarity measure, in place of Pearson correlation or Euclidean distance.  This will
 *	give you a visual representation of what transcriptional patterns look like for genes predicted to be
 *	similar based on their behaviors in your data and their relationships to your four genes of interest.
 *	If you've generated an appropriate combined PCL (\ref sssec_uses_tools_microarray see above), you can
 *	run:
 * \code
 * MCluster -o combined_imputed_fr_predictions.gtr -i cheeto_network.dab
 *		< combined_imputed.pcl > combined_imputed_fr_predictions.cdt
 * \endcode
 * - Or you might like to find dense subgraphs (i.e. clusters) in the predicted functional network using
 *	\ref cliquer Cliquer:
 * \code
 * Cliquer -i cheeto_network.dab -r 3 -w 0.33 > cheeto_network_clusters.txt
 * \endcode
 * which will produce clusters and confidences resembling:
 * \code
 * 16.203  YAR071W YBR093C YHR215W YBR092C YDL106C
 * 14.1403 YBR066C YDR043C YHL027W YDR477W YBR112C YJL089W
 * 15.3548 YBR093C YML077W YML123C YHR136C
 * 14.1379 YIL108W YOR032C YFR028C YGL003C YGR225W
 * 13.2328 YMR238W YOL131W YLR295C YKL046C YDR309C YHR061C YMR055C
 * ...
 * \endcode
 * - Finally, you can assign predicted gene funtions by examining which genes in the network are strongly
 *	connected to known processes.  Say you want to see the top 100 genes predicted to be most associated
 *	with your process of interest, the four active genes listed in \c cheeto_genes.txt.  Then use
 *	\ref hubber Hubber to run:
 * \code
 * Hubber -i cheeto_network.dab -g 100 cheeto_genes.txt > cheeto_gene_predictions.txt
 * \endcode
 * This will result in a list of the 100 genes most confidently predicted to function with your four active
 * genes, of the form:
 * \code
 * YER033C|22.51|0	YDR389W|17.9219|0	YPL204W|16.5705|0 ...
 * \endcode
 * There's plenty more you can do with the Sleipnir tools to analyze interaction networks and datasets of
 * many types - look for genetic hubs, predict associations between pathways, try out different ways of
 * measuring protein domain similarity, or find out what biological processes are activated in your
 * experimental conditions.  To really get into the nitty gritty, keep reading to find out how you can
 * integrate the Sleipnir library into your own computational biology tools and programs.
 * 
 * \subsection ssec_uses_library Core Library
 * 
 * ***
 * 
 * \section sec_philosophy Philosophy
 * 
 * ***
 * 
 * \section sec_dependencies External Dependencies
 * 
 * Although Sleipnir can be built without any additional libraries, it requires three for complete
 * functionality:
 * 
 * - <a href="http://genie.sis.pitt.edu/">SMILE</a>, a graphical models library from the Decision Systems
 *   Laboratory at the University of Pittsburgh.
 * - <a href="http://svmlight.joachims.org/">SVM Light</a>, a support vector machine library by
 *   Thorsten Joachims.
 * - <a href="http://log4cpp.sourceforge.net/">log4cpp</a>, a library of logging and output functions.
 * 
 * In addition to these core libraries, various Sleipnir tools also use the following:
 * 
 * - <a href="http://www.gnu.org/software/gengetopt/gengetopt.html">gengetopt</a>, a GNU tool for
 *   automatically processing command line arguments.
 * - <a href="http://sourceware.org/pthreads-win32/">pthreads-w32</a>, a Windows versions of the standard
 *   Linux pthreads library.
 * - <a href="http://tiswww.case.edu/php/chet/readline/rltop.html">readline</a>, a GNU library for
 *   handling command line history and editing; pre-installed on Linux systems, but necessary for Windows.
 * - <a href="http://www.boost.org/">Boost</a>, a massive library of C++ utilities of which Sleipnir uses
 *   only one tiny (but critical) part in one specific tool.
 * 
 * For information on integrating these tools into your Sleipnir environment, please see the
 * \ref sec_building "Building Sleipnir" section.
 * 
 * \section sec_history Version History
 * 
 * - <a href="sleipnir-0.9.tar.gz">0.9</a>, ***-08 <br>
 * First version made available to reviewers.
 * 
 * \section sec_license License
 * 
 * Sleipnir is provided under the Creative Commons Attribution 3.0 license.
 *
 * You are free to share, copy, distribute, transmit, or adapt this work PROVIDED THAT you attribute the work
 * to the authors listed above.  For more information, please see the following web page:
 * <a href="http://creativecommons.org/licenses/by/3.0/">http://creativecommons.org/licenses/by/3.0/</a>
 */

namespace Sleipnir {

const char	c_szSleipnir[]	= "Sleipnir";
#ifdef USE_LOG4CPP_STUB
Category	g_CatSleipnir;
#else // USE_LOG4CPP_STUB
Category&	g_CatSleipnir	= Category::getInstance( c_szSleipnir );
#endif // USE_LOG4CPP_STUB

}
