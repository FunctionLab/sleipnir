#ifndef ANNOTATION_H
#define ANNOTATION_H

namespace Sleipnir {

/*!
 * \brief
 * Encapsulates the hypergeometric functional enrichment of a query against one ontology term.
 * 
 * Generated by an IOntology::TermFinder call, each STermFound struct represents the hypergeometric
 * enrichment of one ontology term for a set of query genes.  Documentation assumes a set of query
 * genes Q and a set of term genes T.
 * 
 * \see
 * IOntology
 */
struct STermFound {
	/*!
	 * \brief
	 * ID of the ontology term for this enrichment.
	 */
	size_t	m_iID;
	/*!
	 * \brief
	 * Hypergeometric p-value of this enrichment.
	 */
	double	m_dP;
	/*!
	 * \brief
	 * The number of genes in Q intersect T.
	 */
	size_t	m_iHitsTerm;
	/*!
	 * \brief
	 * The number of genes in T.
	 */
	size_t	m_iSizeTerm;
	/*!
	 * \brief
	 * The number of genes in Q.
	 */
	size_t	m_iHitsTotal;
	/*!
	 * \brief
	 * The total number of genes in the background set (genome).
	 */
	size_t	m_iSizeTotal;

	/*!
	 * \brief
	 * Construct a new structure with the given parameter values.
	 * 
	 * \param iID
	 * ID of the ontology term for this enrichment.
	 * 
	 * \param dP
	 * Hypergeometric p-value of this enrichment.
	 * 
	 * \param iHitsTerm
	 * The number of genes in Q intersect T.
	 * 
	 * \param iSizeTerm
	 * The number of genes in T.
	 * 
	 * \param iHitsTotal
	 * The number of genes in Q.
	 * 
	 * \param iSizeTotal
	 * The total number of genes in the background set (genome).
	 */
	STermFound( size_t iID, double dP, size_t iHitsTerm, size_t iSizeTerm, size_t iHitsTotal, size_t iSizeTotal ) :
		m_iID(iID), m_dP(dP), m_iHitsTerm(iHitsTerm), m_iSizeTerm(iSizeTerm), m_iHitsTotal(iHitsTotal),
		m_iSizeTotal(iSizeTotal) { }
};

}

#include "annotationi.h"

namespace Sleipnir {

/*!
 * \brief
 * Encapsulates a functional catalog/hierarchy/ontology such as GO, KEGG, or MIPS.
 * 
 * IOntology provides a uniform interface able to capture ontological functional catalogs such as the
 * Gene Ontology, hierarchical catalogs such as MIPS, and (essentially) flag functional groupings such
 * as KEGG.  IOntology's structure is modeled on GO's, being the most general: each term in the
 * functional catalog has zero or more parents (more general terms), zero or more children (more specific
 * terms), and zero or more directly annotated genes.  A gene annotated to some term is also implicitly
 * annotated to that term's ancestors.
 * 
 * \see
 * COntologyKEGG | COntologyGO | COntologyMIPS | COntologyMIPSPhenotypes | CSlim
 */
class IOntology {
public:
	/*!
	 * \brief
	 * Returns string identifier of the encapsulated ontology.
	 * 
	 * \returns
	 * String identifier of the ontology.
	 */
	virtual const std::string& GetID( ) const = 0;
	/*!
	 * \brief
	 * Returns the number of nodes (terms) in the ontology.
	 * 
	 * \returns
	 * Number of nodes (terms) in the ontology.
	 */
	virtual size_t GetNodes( ) const = 0;
	/*!
	 * \brief
	 * Returns the ontology-specific ID string of the requested term.
	 * 
	 * \param iTerm
	 * Index of ontology term.
	 * 
	 * \returns
	 * Ontology-specific ID of the requested term (e.g. "GO:0007093").
	 */
	virtual const std::string& GetID( size_t iTerm ) const = 0;
	/*!
	 * \brief
	 * Returns the ontology-specific description of the requested term.
	 * 
	 * \param iTerm
	 * Index of ontology term.
	 * 
	 * \returns
	 * Ontology-specific description of the requested term (e.g. "mitotic cell cycle").
	 */
	virtual const std::string& GetGloss( size_t iTerm ) const = 0;
	/*!
	 * \brief
	 * Returns the number of parents of the requested term.
	 * 
	 * \param iTerm
	 * Index of ontology term.
	 * 
	 * \returns
	 * Number of parents of the requested term.
	 */
	virtual size_t GetParents( size_t iTerm ) const = 0;
	/*!
	 * \brief
	 * Returns the ontology term index of the requested parent of the requested term.
	 * 
	 * \param iTerm
	 * Index of ontology term.
	 * 
	 * \param iParent
	 * Parent to retrieve.
	 * 
	 * \returns
	 * Ontology term index of the requested parent.
	 * 
	 * \remarks
	 * Requested parent must be less than IOntology::GetParents.
	 */
	virtual size_t GetParent( size_t iTerm, size_t iParent ) const = 0;
	/*!
	 * \brief
	 * Returns the number of children of the requested term.
	 * 
	 * \param iTerm
	 * Index of ontology term.
	 * 
	 * \returns
	 * Number of children of the requested term.
	 */
	virtual size_t GetChildren( size_t iTerm ) const = 0;
	/*!
	 * \brief
	 * Returns the ontology term index of the requested child of the requested term.
	 * 
	 * \param iTerm
	 * Index of ontology term.
	 * 
	 * \param iChild
	 * Child to retrieve.
	 * 
	 * \returns
	 * Ontology term index of the requested child.
	 * 
	 * \remarks
	 * Requested child must be less than IOntology::GetChildren.
	 */
	virtual size_t GetChild( size_t iTerm, size_t iChild ) const = 0;
	/*!
	 * \brief
	 * Returns the number of genes annotated to or (optionally) below this term.
	 * 
	 * \param iTerm
	 * Ontology term index.
	 * 
	 * \param fRecursive
	 * If true, include gene annotations to descendant terms.
	 * 
	 * \returns
	 * Number of genes annotated to or below this term.
	 */
	virtual size_t GetGenes( size_t iTerm, bool fRecursive = false ) const = 0;
	/*!
	 * \brief
	 * Returns the requested gene annotated to or below the given term.
	 * 
	 * \param iTerm
	 * Ontology term index from which to retrieve gene.
	 * 
	 * \param iGene
	 * Index of gene within the requested term.
	 * 
	 * \returns
	 * Gene annotation at the requested index.
	 * 
	 * \remarks
	 * If iGene is less than the number of genes annotated directly to the term, a direct annotation is
	 * returned.  Otherwise, descendant terms are searched and the requested gene index is retrieved from
	 * this set.  iGene must be less than the value of IOntology::GetGenes.
	 */
	virtual const CGene& GetGene( size_t iTerm, size_t iGene ) const = 0;
	/*!
	 * \brief
	 * Indicates whether the given gene is annotated to or (optionally) below the given term.
	 * 
	 * \param iTerm
	 * Ontology term index.
	 * 
	 * \param Gene
	 * Gene for which annotation is tested.
	 * 
	 * \param fRecursive
	 * If true, include annotations to descendants of the given term.
	 * 
	 * \returns
	 * True if the given gene is annotated to or (optionally) below the given term.
	 */
	virtual bool IsAnnotated( size_t iTerm, const CGene& Gene, bool fRecursive = true ) const = 0;
	/*!
	 * \brief
	 * Returns the ontology term index corresponding to the given ontology-specific ID string.
	 * 
	 * \param strID
	 * Ontology term ID to retrieve (e.g. "GO:0007093").
	 * 
	 * \returns
	 * The index of the requested term, -1 if not found.
	 */
	virtual size_t GetNode( const std::string& strID ) const = 0;
	/*!
	 * \brief
	 * Obtain the primary gene names for all genes in the ontology.
	 * 
	 * \param vecstrGenes
	 * Outputs primary gene names for all genes in the ontology.
	 */
	virtual void GetGeneNames( std::vector<std::string>& vecstrGenes ) const = 0;
	/*!
	 * \brief
	 * Uses the hypergeometric distribution to find functional enrichments of the given gene set.
	 * 
	 * \param Genes
	 * Gene set to test for functional enrichments.
	 * 
	 * \param vecsTerms
	 * Output statistics for the given genes over all terms in the ontology.
	 * 
	 * \param fBonferroni
	 * If true, p-values are Bonferroni corrected.
	 * 
	 * \param fRecursive
	 * If true, annotations to descendants are used when calculating term overlaps.
	 * 
	 * \param fGenome
	 * If true, use all genes in the genome as background; otherwise, use only genes with at least one
	 * annotation in the ontology.
	 * 
	 * \param pBackground
	 * If non-null, use the given gene set as the background.
	 * 
	 * TermFinder uses the GO::TermFinder technique of Boyle, Sherlock, et al to calculate ontology
	 * terms (i.e. functions) enriched in the given gene set.  This involves significance testing the
	 * overlap of the given gene set with annotations to every ontology term using the hypergeometric
	 * test.
	 */
	virtual void TermFinder( const CGenes& Genes, std::vector<STermFound>& vecsTerms, bool fBonferroni = true,
		bool fRecursive = true, bool fGenome = false, const CGenes* pBackground = NULL ) const = 0;
};

// BUGBUG: These should really be templated instead of duplicated like this...

/*!
 * \brief
 * Implements IOntology for the KEGG orthology.
 * 
 * COntologyKEGG parses the "ko" file from the Kyoto Encyclopedia of Genes and Genomes to extract
 * organism-independent function annotations from KEGG.
 */
class COntologyKEGG : COntologyKEGGImpl, public IOntology {
public:
	COntologyKEGG( );
	bool Open( std::istream& istm, CGenome& Genome, const std::string& strOrganism, bool fSynonyms = false );

	void GetGeneNames( std::vector<std::string>& vecstrGenes ) const {

		return COntologyImpl::GetGeneNames( vecstrGenes ); }

	void TermFinder( const CGenes& Genes, std::vector<STermFound>& vecsTerms, bool fBonferroni = true,
		bool fRecursive = true, bool fGenome = false, const CGenes* pBackground = NULL ) const {

		return COntologyImpl::TermFinder( Genes, vecsTerms, fBonferroni, fRecursive, fGenome, pBackground ); }

	size_t GetNode( const std::string& strID ) const {

		return COntologyImpl::GetNode( strID ); }

	bool IsAnnotated( size_t iTerm, const CGene& Gene, bool fRecursive ) const {

		return COntologyImpl::IsAnnotated( iTerm, Gene, fRecursive ); }

	size_t GetNodes( ) const {

		return COntologyImpl::GetNodes( ); }

	const std::string& GetID( ) const {

		return COntologyImpl::GetID( ); }

	const std::string& GetID( size_t iTerm ) const {

		return COntologyImpl::GetID( iTerm ); }

	const std::string& GetGloss( size_t iTerm ) const {

		return COntologyImpl::GetGloss( iTerm ); }

	size_t GetParents( size_t iTerm ) const {

		return COntologyImpl::GetParents( iTerm ); }

	size_t GetParent( size_t iTerm, size_t iParent ) const {

		return COntologyImpl::GetParent( iTerm, iParent ); }

	size_t GetChildren( size_t iTerm ) const {

		return COntologyImpl::GetChildren( iTerm ); }

	size_t GetChild( size_t iTerm, size_t iChild ) const {

		return COntologyImpl::GetChild( iTerm, iChild ); }

	size_t GetGenes( size_t iTerm, bool fRecursive ) const {

		return COntologyImpl::GetGenes( iTerm, fRecursive ); }

	const CGene& GetGene( size_t iTerm, size_t iGene ) const {

		return COntologyImpl::GetGene( iTerm, iGene ); }
};

/*!
 * \brief
 * Implements IOntology for the Gene Ontology.
 * 
 * COntologyGO parses OBO and GO annotation files to obtain the structure and annotations to the Gene
 * Ontology.  Aspects of the ontology (e.g. biological process, molecular function, cellular compartment)
 * can be opened individually to obtain independent IOntology objects.
 */
class COntologyGO : COntologyGOImpl, public IOntology {
public:
	/*!
	 * \brief
	 * Gene Ontology aspect/namespace.
	 */
	enum ENamespace {
		/*!
		 * \brief
		 * Cellular compartment
		 */
		ENamespaceCC,
		/*!
		 * \brief
		 * Biological process
		 */
		ENamespaceBP,
		/*!
		 * \brief
		 * Molecular function
		 */
		ENamespaceMF
	};

	static bool Open( std::istream& istmOntology, std::istream& istmAnnotations, CGenome& Genome,
		COntologyGO& OntoBP, COntologyGO& OntoMF, COntologyGO& OntoCC, bool fDatabaseIDs = false,
		bool fSynonyms = false );

	COntologyGO( );
	bool Open( std::istream& istmOntology, std::istream& istmAnnotations, CGenome& Genome,
		ENamespace eNamespace, bool fDatabaseIDs = false, bool fSynonyms = false );

	void GetGeneNames( std::vector<std::string>& vecstrGenes ) const {

		return COntologyImpl::GetGeneNames( vecstrGenes ); }

	void TermFinder( const CGenes& Genes, std::vector<STermFound>& vecsTerms, bool fBonferroni = true,
		bool fRecursive = true, bool fGenome = false, const CGenes* pBackground = NULL ) const {

		return COntologyImpl::TermFinder( Genes, vecsTerms, fBonferroni, fRecursive, fGenome, pBackground ); }

	size_t GetNode( const std::string& strID ) const {

		return COntologyImpl::GetNode( strID ); }

	bool IsAnnotated( size_t iTerm, const CGene& Gene, bool fRecursive ) const {

		return COntologyImpl::IsAnnotated( iTerm, Gene, fRecursive ); }

	size_t GetNodes( ) const {

		return COntologyImpl::GetNodes( ); }

	const std::string& GetID( ) const {

		return COntologyImpl::GetID( ); }

	const std::string& GetID( size_t iTerm ) const {

		return COntologyImpl::GetID( iTerm ); }

	const std::string& GetGloss( size_t iTerm ) const {

		return COntologyImpl::GetGloss( iTerm ); }

	size_t GetParents( size_t iTerm ) const {

		return COntologyImpl::GetParents( iTerm ); }

	size_t GetParent( size_t iTerm, size_t iParent ) const {

		return COntologyImpl::GetParent( iTerm, iParent ); }

	size_t GetChildren( size_t iTerm ) const {

		return COntologyImpl::GetChildren( iTerm ); }

	size_t GetChild( size_t iTerm, size_t iChild ) const {

		return COntologyImpl::GetChild( iTerm, iChild ); }

	size_t GetGenes( size_t iTerm, bool fRecursive ) const {

		return COntologyImpl::GetGenes( iTerm, fRecursive ); }

	const CGene& GetGene( size_t iTerm, size_t iGene ) const {

		return COntologyImpl::GetGene( iTerm, iGene ); }
};

/*!
 * \brief
 * Implements IOntology for the MIPS functional catalog.
 * 
 * COntologyMIPS parses the "funcat" structure and annotation files from the Munich Information center for
 * Protein Sequences.
 */
class COntologyMIPS : protected COntologyMIPSImpl, public IOntology {
public:
	COntologyMIPS( );
	bool Open( std::istream& istmOntology, std::istream& istmAnnotations, CGenome& Genome );

	void GetGeneNames( std::vector<std::string>& vecstrGenes ) const {

		return COntologyImpl::GetGeneNames( vecstrGenes ); }

	void TermFinder( const CGenes& Genes, std::vector<STermFound>& vecsTerms, bool fBonferroni = true,
		bool fRecursive = true, bool fGenome = false, const CGenes* pBackground = NULL ) const {

		return COntologyImpl::TermFinder( Genes, vecsTerms, fBonferroni, fRecursive, fGenome, pBackground ); }

	size_t GetNode( const std::string& strID ) const {

		return COntologyImpl::GetNode( strID ); }

	bool IsAnnotated( size_t iTerm, const CGene& Gene, bool fRecursive ) const {

		return COntologyImpl::IsAnnotated( iTerm, Gene, fRecursive ); }

	size_t GetNodes( ) const {

		return COntologyImpl::GetNodes( ); }

	const std::string& GetID( ) const {

		return COntologyImpl::GetID( ); }

	const std::string& GetID( size_t iTerm ) const {

		return COntologyImpl::GetID( iTerm ); }

	const std::string& GetGloss( size_t iTerm ) const {

		return COntologyImpl::GetGloss( iTerm ); }

	size_t GetParents( size_t iTerm ) const {

		return COntologyImpl::GetParents( iTerm ); }

	size_t GetParent( size_t iTerm, size_t iParent ) const {

		return COntologyImpl::GetParent( iTerm, iParent ); }

	size_t GetChildren( size_t iTerm ) const {

		return COntologyImpl::GetChildren( iTerm ); }

	size_t GetChild( size_t iTerm, size_t iChild ) const {

		return COntologyImpl::GetChild( iTerm, iChild ); }

	size_t GetGenes( size_t iTerm, bool fRecursive ) const {

		return COntologyImpl::GetGenes( iTerm, fRecursive ); }

	const CGene& GetGene( size_t iTerm, size_t iGene ) const {

		return COntologyImpl::GetGene( iTerm, iGene ); }
};

/*!
 * \brief
 * Extends COntologyMIPS to include the (apparently defunct) "phencat" phenotype hierarchy.
 * 
 * COntologyMIPS parses the "phencat" structure and annotation files from the Munich Information center for
 * Protein Sequences.
 */
class COntologyMIPSPhenotypes : public COntologyMIPS {
public:
	COntologyMIPSPhenotypes( );

protected:
	/*!
	 * \brief
	 * String identifier for the MIPS Phenotype ontology.
	 */
	static const char	c_szMIPSPhen[];
};

/*!
 * \brief
 * Represents a set of ontology terms.
 * 
 * Modeled after the "GO Slim" sets of GO terms, CSlim represents any set of IOntology terms.  A slim is
 * generally opened from a text file listing ontology IDs and can be used in a variety of settings, e.g.
 * functional enrichment or gold standard generation.
 * 
 * \see
 * CDat::Open | CDataPair::Open
 */
class CSlim : CSlimImpl {
public:
	bool Open( std::istream& istmSlim, const IOntology* pOntology );
	void GetGeneNames( std::vector<std::string>& vecstrGenes ) const;

	/*!
	 * \brief
	 * Returns the gene at the requested index below the requested slim term.
	 * 
	 * \param iSlim
	 * Index of the slim term for which a gene are retrieved.
	 * 
	 * \param iGene
	 * Index of the gene to retrieve.
	 * 
	 * \returns
	 * Gene annotated at the requested index below the requested slim term.
	 * 
	 * \remarks
	 * iSlim must be less than CSlim::GetSlims; iGene must be less than CSlim::GetGenes.  All genes annotated
	 * recursively to descendants of the slim terms are considered.
	 * 
	 * \see
	 * IOntology::GetGenes
	 */
	const CGene& GetGene( size_t iSlim, size_t iGene ) const {

		return *m_vecvecpGenes[ iSlim ][ iGene ]; }

	/*!
	 * \brief
	 * Returns the number of ontology terms in this slim.
	 * 
	 * \returns
	 * Number of ontology terms in this slim.
	 */
	size_t GetSlims( ) const {

		return m_vecstrSlims.size( ); }

	/*!
	 * \brief
	 * Returns the number of genes annotated below the given slim term.
	 * 
	 * \param iSlim
	 * Index of the slim term for which genes are retrieved.
	 * 
	 * \returns
	 * Number of genes annotated below the given slim term.
	 * 
	 * \remarks
	 * iSlim must be less than CSlim::GetSlims.
	 */
	size_t GetGenes( size_t iSlim ) const {

		return m_vecvecpGenes[ iSlim ].size( ); }

	/*!
	 * \brief
	 * Returns the string ID of the ontology term at the given slim index.
	 * 
	 * \param iSlim
	 * Index of the slim term to identify.
	 * 
	 * \returns
	 * String ID of the requested slim index.
	 * 
	 * \remarks
	 * iSlim must be less than CSlim::GetSlims.
	 */
	const std::string& GetSlim( size_t iSlim ) const {

		return m_vecstrSlims[ iSlim ]; }
};

}

#endif // ANNOTATION_H