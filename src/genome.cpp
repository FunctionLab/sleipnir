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
#include "genome.h"
#include "meta.h"
#include "annotation.h"

namespace Sleipnir {

/*!
 * \brief
 * Create a new gene with the given primary identifier.
 * 
 * \param strID
 * Unique primary identifier of the new gene.
 * 
 * \see
 * GetName
 */
CGene::CGene( const std::string& strID ) : CGeneImpl(strID) { }

CGeneImpl::CGeneImpl( const std::string& strName ) : m_strName(strName),
	m_iOntologies(0), m_apOntologies(NULL), m_apveciAnnotations(NULL),
	m_iSynonyms(0), m_astrSynonyms(NULL), m_fRNA(false), m_fDubious(false) { }

CGeneImpl::~CGeneImpl( ) {
	size_t	i;

	if( m_iOntologies ) {
		delete[] m_apOntologies;
		for( i = 0; i < m_iOntologies; ++i )
			delete m_apveciAnnotations[ i ];
		delete[] m_apveciAnnotations; }
	if( m_iSynonyms )
		delete[] m_astrSynonyms; }

CGeneImpl& CGeneImpl::operator=( const CGeneImpl& Gene ) {
	size_t	i, j;

	m_strName = Gene.m_strName;
	if( m_iSynonyms = Gene.m_iSynonyms ) {
		m_astrSynonyms = new string[ m_iSynonyms ];
		for( i = 0; i < m_iSynonyms; ++i )
			m_astrSynonyms[ i ] = Gene.m_astrSynonyms[ i ]; }
	if( m_iOntologies = Gene.m_iOntologies ) {
		m_apOntologies = new IOntology const*[ m_iOntologies ];
		m_apveciAnnotations = new vector<size_t>*[ m_iOntologies ];
		for( i = 0; i < m_iOntologies; ++i ) {
			m_apOntologies[ i ] = Gene.m_apOntologies[ i ];
			m_apveciAnnotations[ i ] = new vector<size_t>( );
			m_apveciAnnotations[ i ]->resize( Gene.m_apveciAnnotations[ i ]->size( ) );
			for( j = 0; j < m_apveciAnnotations[ i ]->size( ); ++j )
				(*m_apveciAnnotations[ i ])[ j ] = (*Gene.m_apveciAnnotations[ i ])[ j ]; } }

	return *this; }

/*!
 * \brief
 * Adds an annotation in the gene object to a specific functional ontology term.
 * 
 * \param pOntology
 * Ontology into which the gene should be annotated.
 * 
 * \param iTerm
 * Ontology term to which the gene should be annotated.
 * 
 * \returns
 * True if the annotation was added successfully.
 * 
 * Adds an annotation directly to the current gene based on the given ontology and functional term.
 * This is usually called from from the appropriate IOntology implementation during parsing and is used
 * to provide a local cache of direct annotations to the gene object.
 * 
 * \remarks
 * This does not inform pOntology of the annotation, so if it is not already reflected in that ontology,
 * the annotation will only be visible from the CGene object, not from the IOntology.
 * 
 * \see
 * IsAnnotated
 */
bool CGene::AddAnnotation( const IOntology* pOntology, size_t iTerm ) {
	size_t	i, iOnto;

	for( iOnto = 0; iOnto < m_iOntologies; ++iOnto )
		if( m_apOntologies[ iOnto ] == pOntology )
			break;
	if( iOnto >= m_iOntologies )
		IncrementOntologies( pOntology );

	for( i = 0; i < m_apveciAnnotations[ iOnto ]->size( ); ++i )
		if( (*m_apveciAnnotations[ iOnto ])[ i ] == iTerm )
			return false;
	m_apveciAnnotations[ iOnto ]->push_back( iTerm );
	return true; }

void CGeneImpl::IncrementOntologies( const IOntology* pOntology ) {
	const IOntology**	apOntologies;
	vector<size_t>**	apveciAnnotations;

	apOntologies = new IOntology const*[ m_iOntologies + 1 ];
	if( m_apOntologies ) {
		memcpy( apOntologies, m_apOntologies, m_iOntologies * sizeof(*m_apOntologies) );
		delete[] m_apOntologies; }
	apOntologies[ m_iOntologies ] = pOntology;
	m_apOntologies = apOntologies;

	apveciAnnotations = new vector<size_t>*[ m_iOntologies + 1 ];
	if( m_apveciAnnotations ) {
		memcpy( apveciAnnotations, m_apveciAnnotations, m_iOntologies *
			sizeof(*m_apveciAnnotations) );
		delete[] m_apveciAnnotations; }
	apveciAnnotations[ m_iOntologies++ ] = new vector<size_t>;
	m_apveciAnnotations = apveciAnnotations; }

/*!
 * \brief
 * Return true if the gene is annotated within the given ontology.
 * 
 * \param pOntology
 * Ontology to test for gene annotation.
 * 
 * \returns
 * True if gene is annotated within the given ontology.
 */
bool CGene::IsAnnotated( const IOntology* pOntology ) const {
	size_t	i;

	for( i = 0; i < m_iOntologies; ++i )
		if( pOntology == m_apOntologies[ i ] )
			return true;

	return false; }

/*!
 * \brief
 * Return true if the gene is directly annotated to the given ontology term.
 * 
 * \param pOntology
 * Ontology to test for gene annotation.
 * 
 * \param iTerm
 * Term to test for direct gene annotation.
 * 
 * \returns
 * True if gene is directly annotated to the given ontology term.
 * 
 * \remarks
 * Does not examine the ontology or term descendants directly, so will only return true if the gene
 * has been directly annotated to the given term using AddAnnotation.
 */
bool CGene::IsAnnotated( const IOntology* pOntology, size_t iTerm ) const {
	size_t	i, j;

	for( i = 0; i < m_iOntologies; ++i )
		if( pOntology == m_apOntologies[ i ] ) {
			for( j = 0; j < m_apveciAnnotations[ i ]->size( ); ++j )
				if( iTerm == (*m_apveciAnnotations[ i ])[ j ] )
					return true;
			break; }

	return false; }

/*!
 * \brief
 * Appends the given synonym to the gene's synonym list.
 * 
 * \param strName
 * Synonym to be added.
 * 
 * \returns
 * True if the synonym was successfully added.
 * 
 * \remarks
 * Addition will fail if the synonym is the same as the gene's name or an existing synonym.
 * 
 * \see
 * GetSynonym
 */
bool CGene::AddSynonym( const std::string& strName ) {
	size_t	i;
	string*	astrSynonyms;

	if( !strName.length( ) )
		g_CatSleipnir( ).warn( "CGene::AddSynonym( %s ) adding null synonym to %s",
			strName.c_str( ), m_strName.c_str( ) );
	if( strName == m_strName )
		return false;
	for( i = 0; i < m_iSynonyms; ++i )
		if( strName == m_astrSynonyms[ i ] )
			return false;

	astrSynonyms = new string[ ++m_iSynonyms ];
	for( i = 0; ( i + 1 ) < m_iSynonyms; ++i )
		astrSynonyms[ i ] = m_astrSynonyms[ i ];
	astrSynonyms[ i ] = strName;

	if( m_astrSynonyms )
		delete[] m_astrSynonyms;
	m_astrSynonyms = astrSynonyms;
	return true; }

const char	CGenomeImpl::c_szDubious[]	= "Dubious";
const char	CGenomeImpl::c_szORF[]		= "ORF";
const char*	CGenomeImpl::c_aszRNA[]		= { "ncRNA", "rRNA", "snRNA", "snoRNA", "tRNA",
	NULL };

CGenomeImpl::~CGenomeImpl( ) {
	size_t	i;

	for( i = 0; i < m_vecpGenes.size( ); ++i )
		delete m_vecpGenes[ i ]; }

/*!
 * \brief
 * Construct a new genome by loading the SGD features file.
 * 
 * \param istmFeatures
 * Stream containing the SGD features information.
 * 
 * \returns
 * True if the genome was loaded successfully.
 * 
 * Loads a (presumably yeast) genome from a file formatted as per the SGD features file (SGD_features.tab).
 * This includes gene IDs, synonyms, glosses, and RNA and dubious tags.
 */
bool CGenome::Open( std::istream& istmFeatures ) {
	static const size_t	c_iBuf	= 4096;
	char			szBuf[ c_iBuf ];
	vector<string>	vecstrLine, vecstrNames;
	size_t			i;

	while( istmFeatures.peek( ) != EOF ) {
		istmFeatures.getline( szBuf, c_iBuf - 1 );
		vecstrLine.clear( );
		CMeta::Tokenize( szBuf, vecstrLine );
		if( vecstrLine.size( ) < 16 )
			return false;

		if( vecstrLine[ 1 ] != c_szORF ) {
			for( i = 0; c_aszRNA[ i ]; ++i )
				if( vecstrLine[ 1 ] == c_aszRNA[ i ] )
					break;
			if( !c_aszRNA[ i ] )
				continue; }

		{
			CGene&	Gene	= AddGene( vecstrLine[ 3 ] );

//	1	Type		ORF, CDS, RNA, etc.
//	2	Qualifier	Dubious, Verified, etc.
//	3	ORF
//	4	Name
//	5	Aliases
//	15	DESC
			Gene.SetRNA( vecstrLine[ 1 ] != c_szORF );
			Gene.SetDubious( vecstrLine[ 2 ] == c_szDubious );
			if( vecstrLine[ 4 ].length( ) )
				AddSynonym( Gene, vecstrLine[ 4 ] );
			if( vecstrLine[ 5 ].length( ) ) {
				vecstrNames.clear( );
				CMeta::Tokenize( vecstrLine[ 5 ].c_str( ), vecstrNames, "|" );
				for( i = 0; i < vecstrNames.size( ); ++i )
					AddSynonym( Gene, vecstrNames[ i ] ); }
			Gene.SetGloss( vecstrLine[ 15 ] );
		} }

	return true; }

/*!
 * \brief
 * Constructs a new genome containing the given gene IDs.
 * 
 * \param vecstrGenes
 * Vector of gene IDs to add to the new genome.
 * 
 * \returns
 * True if the genome was created successfully.
 * 
 * \remarks
 * Genes in the new genome will have no information beyond the provided primary IDs, which should (as
 * usual) be unique.
 */
bool CGenome::Open( const std::vector<std::string>& vecstrGenes ) {
	size_t	i;

	for( i = 0; i < vecstrGenes.size( ); ++i )
		AddGene( vecstrGenes[ i ] );

	return true; }

/*!
 * \brief
 * Adds a new gene with the given primary ID to the genome.
 * 
 * \param strID
 * Gene ID to be added to the genome.
 * 
 * \returns
 * A reference to the newly added gene, or to an existing gene with the given name.
 * 
 * Given a gene name, AddGene will first test to see if any gene in the genome has that ID or synonym; if so,
 * a reference to the existing gene is returned.  Otherwise, an empty gene with the given primary ID is
 * created, and a reference to this new gene is returned.
 * 
 * \remarks
 * A newly created gene will have no information beyond the provided primary ID.
 * 
 * \see
 * FindGene
 */
CGene& CGenome::AddGene( const std::string& strID ) {
	TMapStrI::const_iterator	iterGene;
	CGene*						pGene;

	if( ( iterGene = m_mapGenes.find( strID ) ) != m_mapGenes.end( ) )
		return GetGene( iterGene->second );

	pGene = new CGene( strID );
	m_vecpGenes.push_back( pGene );
	m_mapGenes[ strID ] = m_vecpGenes.size( ) - 1;
	return *pGene; }

/*!
 * \brief
 * Return the index of a gene within the genome, or -1 if it does not exist.
 * 
 * \param strGene
 * Name of gene to be retrieved from the genome.
 * 
 * \returns
 * Index of the requested gene, or -1 if it does not exist.
 * 
 * Search the genome's gene list for a gene with the given name, primary or synonymous, and return its
 * index if found.
 * 
 * \remarks
 * Both the genome's internal name map and the synonyms of every gene are explicitly searched; the latter
 * can be very slow, and the internal map will not always contain synonyms (depending on how the genome
 * was constructed).
 * 
 * \see
 * AddGene
 */
size_t CGenome::FindGene( const std::string& strGene ) const {
	size_t	i, j, iRet;

	if( ( iRet = GetGene( strGene ) ) != -1 )
		return iRet;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = 0; j < GetGene( i ).GetSynonyms( ); ++j )
			if( strGene == GetGene( i ).GetSynonym( j ) )
				return i;

	return -1; }

/*!
 * \brief
 * Return a vector of all primary gene IDs in the genome.
 * 
 * \returns
 * Vector of all primary gene IDs in the genome.
 */
vector<string> CGenome::GetGeneNames( ) const {
	vector<string>	vecstrRet;
	size_t			i;

	vecstrRet.resize( m_vecpGenes.size( ) );
	for( i = 0; i < vecstrRet.size( ); ++i )
		vecstrRet[ i ] = m_vecpGenes[ i ]->GetName( );

	return vecstrRet; }

/*!
 * \brief
 * Returns the number of genes in the genome with annotations in the given ontology.
 * 
 * \param pOntology
 * Ontology to be scanned for annotated genes.
 * 
 * \returns
 * Number of genes in the genome with annotations in the given ontology.
 * 
 * \see
 * CGene::GetOntology
 */
size_t CGenome::CountGenes( const IOntology* pOntology ) const {
	size_t	i, j, iRet;

	for( iRet = i = 0; i < m_vecpGenes.size( ); ++i )
		for( j = 0; j < m_vecpGenes[ i ]->GetOntologies( ); ++j )
			if( pOntology == m_vecpGenes[ i ]->GetOntology( j ) ) {
				iRet++;
				break; }

	return iRet; }

/*!
 * \brief
 * Explicitly add a gene synonym to the gene and to the genome's name map.
 * 
 * \param Gene
 * Gene to which synonym is to be added.
 * 
 * \param strName
 * Synonym to be added to the given gene.
 * 
 * \returns
 * True if the synonym was added successfully.
 * 
 * \remarks
 * Addition will fail if the synonym is the given gene's primary ID or an existing synonym.
 * 
 * \see
 * CGene::AddSynonym
 */
bool CGenome::AddSynonym( CGene& Gene, const std::string& strName ) {

	if( ( strName != Gene.GetName( ) ) && ( Gene.AddSynonym( strName ) ) ) {
		m_mapGenes[ strName ] = m_mapGenes[ Gene.GetName( ) ];
		return true; }

	return false; }

/*!
 * \brief
 * Simultaneously construct multiple new gene sets loaded from the given file, one per line, with tab-delimited genes.
 * 
 * \param szFile
 * File from which gene sets are loaded.
 * 
 * \param Genome
 * Genome containing all genes which might become members of these gene sets.
 * 
 * \param vecstrNames
 * Human-readable identifiers for the loaded gene sets.
 * 
 * \param vecpGenes
 * Vector to which loaded gene sets are appended.
 * 
 * \returns
 * True on success, false otherwise.
 * 
 * Opens multiple gene sets from the given tab-delimited text file.  Each line should contain a single tab-delimited gene
 * set, and the first token on each line should be a human-readable identifier for that line's gene set.
 * 
 * \see
 * Open
 */
bool CGenes::Open( const char* szFile, CGenome& Genome, std::vector<std::string>& vecstrNames, std::vector<CGenes*>& vecpGenes ) {
	ifstream		ifsm;
	vector<char>	veccBuffer;
	istream*		pistm;
	bool			fRet;

	if( szFile ) {
		ifsm.open( szFile );
		pistm = &ifsm; }
	else
		pistm = &cin;
	if( !ifsm.is_open( ) ) {
		g_CatSleipnir( ).error( "CGenes::Open( %s ) could not open file", szFile ? szFile : "stdin" );
		return false; }

	veccBuffer.resize( CFile::GetBufferSize( ) );
	fRet = false;
	while( !pistm->eof( ) ) {
		vector<string>	vecstrLine;

		pistm->getline( &veccBuffer[0], veccBuffer.size( ) - 1 );
		CMeta::Tokenize( &veccBuffer[0], vecstrLine );
		if( vecstrLine.empty( ) )
			continue;
		vecstrNames.push_back( vecstrLine[0] );
		vecstrLine.erase( vecstrLine.begin( ) );
		vecpGenes.push_back( new CGenes( Genome ) );
		if( !( fRet = vecpGenes.back( )->Open( vecstrLine ) ) )
			break; }
	if( szFile )
		ifsm.close( );

	return fRet; }

/*!
 * \brief
 * Construct a new gene set containing genomes drawn from the given underlying genome.
 * 
 * \param Genome
 * Genome containing all genes which might become members of this gene set.
 */
CGenes::CGenes( CGenome& Genome ) : CGenesImpl( Genome ) { }

CGenesImpl::CGenesImpl( CGenome& Genome ) : m_Genome(Genome) { }

/*!
 * \brief
 * Construct a new gene set by loading genes from the given text stream, one per line.
 * 
 * \param istm
 * Stream containing gene IDs to load, one per line.
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
bool CGenes::Open( std::istream& istm, bool fCreate ) {
	static const size_t	c_iBuffer	= 1024;
	char	szBuf[ c_iBuffer ];
	CGene*	pGene;
	size_t	i, iGene;
	char*	pc;

	if( istm.rdstate( ) != ios_base::goodbit )
		return false;

	m_mapGenes.clear( );
	m_vecpGenes.clear( );
	while( istm.peek( ) != EOF ) {
		istm.getline( szBuf, c_iBuffer - 1 );
		if( pc = strchr( szBuf, '\t' ) )
			*pc = 0;
		if( !szBuf[ 0 ] || ( szBuf[ 0 ] == c_cComment ) )
			continue;
		if( fCreate )
			pGene = &m_Genome.AddGene( szBuf );
		else {
			if( ( iGene = m_Genome.FindGene( szBuf ) ) == -1 ) {
				g_CatSleipnir( ).warn( "CGenes::Open( %d ) unknown gene: %s", fCreate, szBuf );
				continue; }
			pGene = &m_Genome.GetGene( iGene ); }
		for( i = 0; i < m_vecpGenes.size( ); ++i )
			if( m_vecpGenes[ i ] == pGene )
				break;
		if( i != m_vecpGenes.size( ) )
			continue;
		m_mapGenes[ pGene->GetName( ) ] = m_vecpGenes.size( );
		m_vecpGenes.push_back( pGene ); }

	return true; }

/*!
 * \brief
 * Return the number of genes in the set annotated at or, optionally, below the given ontology term.
 * 
 * \param pOntology
 * Ontology in which annotations are counted.
 * 
 * \param iTerm
 * Ontology term at or below which annotations are counted.
 * 
 * \param fRecursive
 * If true, count annotations at or below the given term; otherwise, count only direct annotations to the
 * term.
 * 
 * \param pBackground
 * If non-null, count only annotations for genes also contained in the given background set.
 * 
 * \returns
 * Number of genes in the gene set annotated at or below the given ontology term.
 * 
 * \see
 * IOntology::IsAnnotated
 */
size_t CGenes::CountAnnotations( const IOntology* pOntology, size_t iTerm, bool fRecursive,
	const CGenes* pBackground ) const {
	size_t	i, iRet;

	for( iRet = i = 0; i < m_vecpGenes.size( ); ++i )
		if( ( !pBackground || pBackground->IsGene( m_vecpGenes[ i ]->GetName( ) ) ) &&
			pOntology->IsAnnotated( iTerm, *m_vecpGenes[ i ], fRecursive ) )
			iRet++;

	return iRet; }

/*!
 * \brief
 * Construct a new gene set containing the given gene IDs.
 * 
 * \param vecstrGenes
 * Primary identifiers of genes in the new gene set.
 * 
 * \param fCreate
 * If true, add unknown genes to the underlying genome; otherwise, unknown gene IDs are ignored.
 * 
 * \returns
 * True if gene set was constructed successfully.
 * 
 * If the given gene identifiers are found in the gene set's underlying genome, CGene objects are loaded
 * from there.  Otherwise, if fCreate is true, new genes are created from the loaded IDs.  If fCreate is
 * false, unrecognized genes are skipped with a warning.
 * 
 * \see
 * CGenome::AddGene
 */
bool CGenes::Open( const std::vector<std::string>& vecstrGenes, bool fCreate ) {
	size_t	i, iGene;
	CGene*	pGene;

	m_mapGenes.clear( );
	m_vecpGenes.clear( );
	for( i = 0; i < vecstrGenes.size( ); ++i ) {
		if( !fCreate && ( ( iGene = m_Genome.FindGene( vecstrGenes[ i ] ) ) == -1 ) )
			continue;

		pGene = fCreate ? &m_Genome.AddGene( vecstrGenes[ i ] ) : &m_Genome.GetGene( iGene );
		m_mapGenes[ vecstrGenes[ i ] ] = m_vecpGenes.size( );
		m_vecpGenes.push_back( pGene ); }

	return true; }

/*!
 * \brief
 * Remove the given genes from the gene set.
 * 
 * \param GenesExclude
 * Genes to be removed from the current gene set.
 * 
 * \remarks
 * Comparisons are performed using pointers to CGene objects, so both gene sets should use the same
 * underlying CGenome for proper behavior.
 */
void CGenes::Filter( const CGenes& GenesExclude ) {
	size_t	i, j, iSize;

	iSize = m_vecpGenes.size( );
	for( i = 0; i < GenesExclude.GetGenes( ); ++i )
		for( j = 0; j < iSize; ++j )
			if( m_vecpGenes[ j ] == &GenesExclude.GetGene( i ) ) {
				m_mapGenes.erase( m_vecpGenes[ j ]->GetName( ) );
				m_vecpGenes[ j ] = m_vecpGenes[ m_vecpGenes.size( ) - 1 ];
				iSize--;
				break; }
	m_vecpGenes.resize( iSize ); }

/*!
 * \brief
 * Return the primary identifiers of all genes in the set.
 * 
 * \returns
 * Vector of primary identifiers of all genes in the set.
 */
vector<string> CGenes::GetGeneNames( ) const {
	vector<string>	vecstrRet;
	size_t			i;

	vecstrRet.resize( m_vecpGenes.size( ) );
	for( i = 0; i < vecstrRet.size( ); ++i )
		vecstrRet[ i ] = m_vecpGenes[ i ]->GetName( );

	return vecstrRet; }

}
