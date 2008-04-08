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
#ifndef GENOME_H
#define GENOME_H

#include <fstream>
#include <string>

#include "genomei.h"

namespace Sleipnir {

/*!
 * \brief
 * Represents a gene with a single primary unique identifier, zero or more synonyms, and zero or more
 * functional annotations.
 * 
 * Uniquely identifying genes is a daunting task in many organisms; in yeast, almost every gene has a single
 * transcript, a single protein, a single systematic name, and at most a few common names.  In higher
 * organisms, multiple transcripts, multiple naming schemes, and a lack of standards make consistent
 * gene identification nearly impossible.
 * 
 * Fortunately or unfortunately, Sleipnir ignores most of this complexity.  A Sleipnir gene consists of at
 * least one unique systematic name - but that's it.  A gene's primary identifier should be unique
 * within its genome, and it can have zero or more (optionally unique) synonyms; CGenome does its best to
 * disambiguate overlapping synonyms, but the primary identifier must be unique within a gene set of
 * interest.
 * 
 * Each gene can also have an optional textual description (gloss), zero or more functional annotations
 * tying it to IOntology objects, and reflective of Sleipnir's roots in yeast biology, genes can also be
 * decorated with a few features indicating whether they encode dubious ORFs or RNA genes.
 * 
 * \remarks
 * For clarity, this documentation will refer to a gene's unique primary identifier as its ID and to its
 * optional synonyms as either synonyms or simply names.
 * 
 * \see
 * CGenes
 */
class CGene : CGeneImpl {
public:
	CGene( const std::string& strID );

	bool AddSynonym( const std::string& strName );
	bool AddAnnotation( const IOntology* pOntology, size_t iTerm );
	bool IsAnnotated( const IOntology* pOntology ) const;
	bool IsAnnotated( const IOntology* pOntology, size_t iTerm ) const;

	/*!
	 * \brief
	 * Return the gene's primary identifier.
	 * 
	 * \returns
	 * Gene's primary identifier.
	 */
	const std::string& GetName( ) const {

		return m_strName; }

	/*!
	 * \brief
	 * Return the number of synonyms of the current gene.
	 * 
	 * \returns
	 * Current gene's number of synonyms.
	 * 
	 * \see
	 * GetSynonym
	 */
	size_t GetSynonyms( ) const {

		return m_iSynonyms; }

	/*!
	 * \brief
	 * Return the requested synonym of the current gene.
	 * 
	 * \param iSynonym
	 * Synonym to retrieve.
	 * 
	 * \returns
	 * Requested synonym.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetSynonyms.
	 */
	const std::string& GetSynonym( size_t iSynonym ) const {

		return m_astrSynonyms[ iSynonym ]; }

	/*!
	 * \brief
	 * Sets the human-readable text description of the gene.
	 * 
	 * \param strGloss
	 * Text description of the gene.
	 * 
	 * \see
	 * GetGloss
	 */
	void SetGloss( const std::string& strGloss ) {

		m_strGloss = strGloss; }

	/*!
	 * \brief
	 * Return the human-readable text description of the gene.
	 * 
	 * \returns
	 * Human-readable text description of the gene.
	 * 
	 * \see
	 * SetGloss
	 */
	const std::string& GetGloss( ) const {

		return m_strGloss; }

	/*!
	 * \brief
	 * Sets whether the gene encodes a dubious reading frame.
	 * 
	 * \param fDubious
	 * Value to store for dubious reading frame status.
	 * 
	 * \see
	 * GetDubious
	 */
	void SetDubious( bool fDubious ) {

		m_fDubious = fDubious; }

	/*!
	 * \brief
	 * Return whether the gene encodes a dubious reading frame.
	 * 
	 * \returns
	 * Dubious reading frame status.
	 * 
	 * \see
	 * SetDubious
	 */
	bool GetDubious( ) const {

		return m_fDubious; }

	/*!
	 * \brief
	 * Sets whether the gene encodes an RNA gene.
	 * 
	 * \param fRNA
	 * Value to store for RNA gene status.
	 * 
	 * \see
	 * GetRNA
	 */
	void SetRNA( bool fRNA ) {

		m_fRNA = fRNA; }

	/*!
	 * \brief
	 * Return whether the gene encodes an RNA gene.
	 * 
	 * \returns
	 * RNA gene status.
	 * 
	 * \see
	 * SetRNA
	 */
	bool GetRNA( ) const {

		return m_fRNA; }

	/*!
	 * \brief
	 * Return the number of different ontologies in which this gene is annotated.
	 * 
	 * \returns
	 * Number of ontologies in which this gene is annotated.
	 * 
	 * \see
	 * GetOntology
	 */
	size_t GetOntologies( ) const {

		return m_iOntologies; }

	/*!
	 * \brief
	 * Return the ontology at the requested index.
	 * 
	 * \param iOntology
	 * Index of ontology to return.
	 * 
	 * \returns
	 * Ontology at requested index.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetOntologies.
	 */
	const IOntology* GetOntology( size_t iOntology ) const {

		return m_apOntologies[ iOntology ]; }

	/*!
	 * \brief
	 * Return the number of annotations of this gene within the requested ontology.
	 * 
	 * \param iOntology
	 * Ontology whose annotations should be counted.
	 * 
	 * \returns
	 * Number of annotations of this gene within the requested ontology.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetOntologies.
	 */
	size_t GetAnnotations( size_t iOntology ) const {

		return m_apveciAnnotations[ iOntology ]->size( ); }

	/*!
	 * \brief
	 * Return the ontology term index of the requested annotation.
	 * 
	 * \param iOntology
	 * Ontology index of requested annotation.
	 * 
	 * \param iAnnotation
	 * Requested annotation index.
	 * 
	 * \returns
	 * Ontology term of the requested annotation.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given indices must be smaller than GetOntologies
	 * and GetAnnotations.
	 */
	size_t GetAnnotation( size_t iOntology, size_t iAnnotation ) const {

		return (*m_apveciAnnotations[ iOntology ])[ iAnnotation ]; }
};

/*!
 * \brief
 * Organizes a collection of unique genes representing a background or maximum gene set for some situation.
 * 
 * Ideally, a genome represents a collection of all known genes for some organism, each with a single unique
 * identifier and some number of non-overlapping synonyms.  In practice, this doesn't happen: a genome often
 * represents a background or comprehensive gene set for some situation (e.g. functional enrichment), or the
 * total set of genes in some data file or analysis (e.g. a functional catalog).  CGenome will do its
 * best to disambiguate overlapping gene names, but it boils down to a simple one-to-one map, which will
 * not deal with ambiguous synonyms accurately.  For best results, guarantee that each gene has a unique
 * primary identifier that does not overlap with any synonyms, and look up genes using only those identifiers.
 * 
 * \remarks
 * CGenome uses an internal map to look up genes given a name.  Both primary identifiers and synonyms are
 * placed in this map, which allows synonyms to be looked up rapidly, but at the price of potentially
 * screwing up in any case with overlapping names.  It's a tough problem; Sleipnir's solution is, as usual,
 * to favor efficiency at the expense of restricting input format (i.e. find some way to uniquely identify
 * your genes).
 * 
 * \see
 * CGene | CGenes
 */
class CGenome : CGenomeImpl {
public:
	bool Open( std::istream& istmFeatures );
	bool Open( const std::vector<std::string>& vecstrGenes );
	CGene& AddGene( const std::string& strID );
	size_t FindGene( const std::string& strGene ) const;
	std::vector<std::string> GetGeneNames( ) const;
	size_t CountGenes( const IOntology* pOntology ) const;
	bool AddSynonym( CGene& Gene, const std::string& strName );

	/*!
	 * \brief
	 * Return the gene at the requested index within the genome.
	 * 
	 * \param iGene
	 * Index of gene to retrieve.
	 * 
	 * \returns
	 * Gene at the requested index.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given value must be smaller than GetGenes.
	 */
	CGene& GetGene( size_t iGene ) const {

		return *m_vecpGenes[ iGene ]; }

	/*!
	 * \brief
	 * Return the index of the gene with the given name, or -1 if one cannot be found.
	 * 
	 * \param strGene
	 * Name of gene whose index should be retrieved.
	 * 
	 * \returns
	 * Index of gene with the given name; -1 if one cannot be found.
	 * 
	 * \remarks
	 * This is sensitive to the same name mapping issues discussed in AddSynonym.  For maximum safety,
	 * refer to genes only by their unique primary identifiers.
	 */
	size_t GetGene( const std::string& strGene ) const {
		TMapStrI::const_iterator	iterGene;

		return ( ( ( iterGene = m_mapGenes.find( strGene ) ) == m_mapGenes.end( ) ) ? -1 :
			iterGene->second ); }

	/*!
	 * \brief
	 * Return the number of genes in the genome.
	 * 
	 * \returns
	 * Number of genes in the genome.
	 */
	size_t GetGenes( ) const {

		return m_vecpGenes.size( ); }
};

/*!
 * \brief
 * Represents a simple set of unique genes.
 * 
 * \remarks
 * Genes are represented by index only and not explicitly checked for uniqueness, so most of the naming
 * issues of CGenome are avoided.  Gene comparisons generally assume a constant gene pool drawn from the
 * base CGenome and are thus performed by pointer comparisons for efficiency; in other words, don't
 * expect two different CGene objects with the same primary ID to behave correctly.
 */
class CGenes : CGenesImpl {
public:
	CGenes( CGenome& Genome );

	bool Open( std::istream& istm, bool fCreate = true );
	bool Open( const std::vector<std::string>& vecstrGenes, bool fCreate = true );
	void Filter( const CGenes& GenesExclude );
	size_t CountAnnotations( const IOntology* pOntology, size_t iTerm, bool fRecursive = true,
		const CGenes* pBackground = NULL ) const;
	std::vector<std::string> GetGeneNames( ) const;

	/*!
	 * \brief
	 * Construct a new gene set by loading genes from the given text file, one per line.
	 * 
	 * \param szFile
	 * File containing gene IDs to load, one per line.
	 * 
	 * \param fCreate
	 * If true, add unknown genes to the underlying genome; otherwise, unknown gene IDs are ignored.
	 * 
	 * \returns
	 * True if gene set was constructed successfully.
	 * 
	 * Loads a text file of the form:
	 * \code
	 * GENE1
	 * GENE2
	 * GENE3
	 * \endcode
	 * containing one primary gene identifier per line.  If these gene identifiers are found in the gene set's
	 * underlying genome, CGene objects are loaded from there.  Otherwise, if fCreate is true, new genes are
	 * created from the loaded IDs.  If fCreate is false, unrecognized genes are skipped with a warning.
	 * 
	 * \see
	 * CGenome::AddGene
	 */
	bool Open( const char* szFile, bool fCreate = true ) {
		std::ifstream	ifsm;

		ifsm.open( szFile );
		return ( ifsm.is_open( ) && Open( ifsm, fCreate ) ); }

	/*!
	 * \brief
	 * Return the number of genes in the set.
	 * 
	 * \returns
	 * Number of genes in the set.
	 */
	size_t GetGenes( ) const {

		return m_vecpGenes.size( ); }

	/*!
	 * \brief
	 * Return true if the given name is a primary identifier of a gene in the set.
	 * 
	 * \param strGene
	 * Primary gene identifier for which the set is searched.
	 * 
	 * \returns
	 * True if the set contains a gene with the given primary identifier.
	 * 
	 * \see
	 * GetGene
	 */
	bool IsGene( const std::string& strGene ) const {

		return ( m_mapGenes.find( strGene ) != m_mapGenes.end( ) ); }

	/*!
	 * \brief
	 * Return the gene set's underlying genome.
	 * 
	 * \returns
	 * Gene set's underlying genome.
	 */
	CGenome& GetGenome( ) const {

		return m_Genome; }

	/*!
	 * \brief
	 * Return the gene at the requested index.
	 * 
	 * \param iGene
	 * Gene index to retrieve.
	 * 
	 * \returns
	 * Gene at the requested index.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetGenes.
	 */
	const CGene& GetGene( size_t iGene ) const {

		return *m_vecpGenes[ iGene ]; }

	/*!
	 * \brief
	 * Return the index of the gene with the given primary identifier, or -1 if none exists.
	 * 
	 * \param strGene
	 * Primary gene identifier for which the set is searched.
	 * 
	 * \returns
	 * Index of the gene with the given primary identifier; -1 if none exists.
	 * 
	 * \see
	 * IsGene
	 */
	size_t GetGene( const std::string& strGene ) const {
		TMapStrI::const_iterator	iterGene;

		return ( ( ( iterGene = m_mapGenes.find( strGene ) ) == m_mapGenes.end( ) ) ? -1 :
			iterGene->second ); }
};

}

#endif // GENOME_H
