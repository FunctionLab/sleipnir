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
#include "annotation.h"
#include "genome.h"
#include "meta.h"

namespace Sleipnir {

const char	COntologyOBO::c_szBiologicalProcess[]	= "biological_process";
const char	COntologyOBO::c_szCellularComponent[]	= "cellular_component";
const char	COntologyOBO::c_szMolecularFunction[]	= "molecular_function";

const char	COntologyOBOImpl::c_szAltID[]		= "alt_id: ";
const char	COntologyOBOImpl::c_szID[]			= "id: ";
const char	COntologyOBOImpl::c_szIsA[]			= "is_a: ";
const char	COntologyOBOImpl::c_szIsObsolete[]	= "is_obsolete: ";
const char	COntologyOBOImpl::c_szOBO[]			= "OBO";
const char	COntologyOBOImpl::c_szHUMAN[]		= "_HUMAN";
const char	COntologyOBOImpl::c_szName[]			= "name: ";
const char	COntologyOBOImpl::c_szNamespace[]	= "namespace: ";
const char	COntologyOBOImpl::c_szNOT[]			= "NOT";
const char	COntologyOBOImpl::c_szPartOf[]		= "part_of ";
const char	COntologyOBOImpl::c_szRelationship[]	= "relationship: ";
const char	COntologyOBOImpl::c_szSGD[]			= "SGD";
const char	COntologyOBOImpl::c_szTerm[]			= "[Term]";

COntologyOBO::COntologyOBO( ) {

	m_pOntology = this; }

COntologyOBOImpl::SParserOBO::SParserOBO( std::istream& istm, CGenome& Genome, bool fDBIDs, bool fSynonyms ) :
	m_fObsolete(false), m_szTarget(NULL), m_fDBIDs(fDBIDs), m_fSynonyms(fSynonyms), SParser( istm, Genome ) { }

void COntologyOBOImpl::SParserOBO::Reset( ) {

	m_fObsolete = false;
	m_vecstrIDs.clear( ); }

COntologyOBOImpl::COntologyOBOImpl( ) : COntologyImpl( c_szOBO ) { }

/*!
 * \brief
 * Initializes three aspects of the Gene Ontology using OBO structure and GO annotation files.
 * 
 * \param istmOntology
 * Stream from which OBO structure file is read.
 * 
 * \param istmAnnotations
 * Stream from which gene annotation file is read.
 * 
 * \param Genome
 * Genome into which genes are inserted or read during annotation parsing.
 * 
 * \param OntoBP
 * Ontology initialized with GO's biological process aspect.
 * 
 * \param OntoMF
 * Ontology initialized with GO's molecular function aspect.
 * 
 * \param OntoCC
 * Ontology initialized with GO's cellular compartment aspect.
 * 
 * \param fDatabaseIDs
 * If true, use annotation database IDs as primary gene names.
 * 
 * \param fSynonyms
 * If true, use first gene synonym (if present) as primary gene names.
 * 
 * \returns
 * True if the ontologies were successfully initialized.
 * 
 * Equivalent to calling COntologyOBO::Open three times for the three GO aspects.
 */
bool COntologyOBO::Open( std::istream& istmOntology, std::istream& istmAnnotations, CGenome& Genome,
	COntologyOBO& OntoBP, COntologyOBO& OntoMF, COntologyOBO& OntoCC, bool fDatabaseIDs, bool fSynonyms ) {

	if( !OntoBP.Open( istmOntology, istmAnnotations, Genome, c_szBiologicalProcess, fDatabaseIDs, fSynonyms ) )
		return false;

	istmOntology.clear( );
	istmOntology.seekg( 0, ios_base::beg );
	istmAnnotations.clear( );
	istmAnnotations.seekg( 0, ios_base::beg );
	if( !OntoCC.Open( istmOntology, istmAnnotations, Genome, c_szCellularComponent, fDatabaseIDs, fSynonyms ) )
		return false;

	istmOntology.clear( );
	istmOntology.seekg( 0, ios_base::beg );
	istmAnnotations.clear( );
	istmAnnotations.seekg( 0, ios_base::beg );
	return OntoMF.Open( istmOntology, istmAnnotations, Genome, c_szMolecularFunction, fDatabaseIDs, fSynonyms ); }

/*!
 * \brief
 * Initializes the ontology from one aspect of the Gene Ontology using OBO structure and GO annotation files.
 * 
 * \param istmOntology
 * Stream from which OBO structure file is read.
 * 
 * \param istmAnnotations
 * Stream from which gene annotation file is read.
 * 
 * \param Genome
 * Genome into which genes are inserted or read during annotation parsing.
 * 
 * \param szNamespace
 * Aspect (namespace) of GO to use for this ontology.
 * 
 * \param fDatabaseIDs
 * If true, use annotation database IDs as primary gene names.
 * 
 * \param fSynonyms
 * If true, use first gene synonym (if present) as primary gene names.
 * 
 * \returns
 * True if the ontology was successfully initialized.
 * 
 * COntologyOBO::Open parses the structure of one aspect of the Gene Ontology from the given
 * (organism-independent) OBO file and obtains gene annotations from the given (organism-specific) annotation
 * file.  Terms are identified by GO IDs (e.g. "GO:0007093"), and genes are identified by the annotation name,
 * first synonym, or database ID as specified by the fSynonyms and fDatabaseIDs parameters.  Genes are
 * retrieved from Genome if already present or inserted if not; it is thus important to ensure that the proper
 * primary gene names are used so as to agree with any identifiers already present in Genome.
 * 
 * \remarks
 * For yeast, fSynonyms and fDatabaseIDs should both be false to use ORF IDs as primary gene names (false will
 * use common names).  For human, it is recommended to pre-map annotation IDs to HGNC symbols, as the default
 * IPI and Uniprot peptide IDs can be extremely ambiguous.  For mouse, fDatabaseIDs should be true to use MGI
 * IDs as primary gene names.  For worm, fDatabaseIDs should be true to use Wormbase gene IDs as primary gene
 * names, or fSynonyms should be true to use systematic transcript IDs (both false will use common names).
 * For fly, fDatabaseIDs should be true to use Flybase gene IDs as primary gene names.
 * 
 * \see
 * COntologyOBO::Open
 */
bool COntologyOBO::Open( std::istream& istmOntology, std::istream& istmAnnotations, CGenome& Genome,
	const char* szNamespace, bool fDatabaseIDs, bool fSynonyms ) {
	SParserOBO	sParserOnto( istmOntology, Genome );
	SParserOBO	sParserGene( istmAnnotations, Genome, fDatabaseIDs, fSynonyms );

	m_strID += szNamespace;
	sParserOnto.m_szTarget = szNamespace;
	return ( OpenOntology( sParserOnto ) && OpenGenes( sParserGene ) ); }

bool COntologyOBOImpl::OpenOntology( SParserOBO& sParser ) {
	size_t					i, j, iParent;
	vector<vector<size_t> >	vecveciChildren;

	g_CatSleipnir( ).info( "COntologyOBOImpl::OpenOntology( ) %s", sParser.m_szTarget );
	Reset( );
	sParser.m_vecvecstrParents.clear( );
	sParser.m_vecNodes.clear( );
	if( !( sParser.GetLine( ) && OpenHeader( sParser ) ) )
		return false;

	while( sParser.m_istm.peek( ) != EOF )
		if( !OpenBlock( sParser ) )
			return false;

	m_aNodes = new SNode[ m_iNodes = sParser.m_vecNodes.size( ) ];
	vecveciChildren.resize( m_iNodes );
	for( i = 0; i < m_iNodes; ++i ) {
		m_aNodes[ i ] = sParser.m_vecNodes[ i ];
		if( m_aNodes[ i ].m_iParents = sParser.m_vecvecstrParents[ i ].size( ) ) {
			m_aNodes[ i ].m_aiParents = new size_t[ m_aNodes[ i ].m_iParents ];
			for( j = 0; j < m_aNodes[ i ].m_iParents; ++j ) {
				if( ( iParent = m_mapNodes[ sParser.m_vecvecstrParents[ i ][ j ] ] ) == i ) {
					g_CatSleipnir( ).error( "COntologyOBOImpl::OpenOntology( ) found a loop for node %d: %s has parent %s",
						i, m_aNodes[ i ].m_strID.c_str( ), sParser.m_vecvecstrParents[ i ][ j ].c_str( ) );
					return false; }
				m_aNodes[ i ].m_aiParents[ j ] = iParent;
				vecveciChildren[ m_aNodes[ i ].m_aiParents[ j ] ].push_back( i ); } } }
	for( i = 0; i < m_iNodes; ++i ) {
		if( !vecveciChildren[ i ].size( ) )
			continue;
		m_aNodes[ i ].m_aiChildren = new size_t[ m_aNodes[ i ].m_iChildren =
			vecveciChildren[ i ].size( ) ];
		for( j = 0; j < m_aNodes[ i ].m_iChildren; ++j )
			m_aNodes[ i ].m_aiChildren[ j ] = vecveciChildren[ i ][ j ]; }

	return true; }

bool COntologyOBOImpl::OpenHeader( SParserOBO& sParser ) {

	while( sParser.m_szLine[ 0 ] )
		if( !sParser.GetLine( ) )
			return false;

	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenBlock( SParserOBO& sParser ) {

	if( sParser.IsStart( c_szTerm ) )
		return ( sParser.GetLine( ) && OpenTerm( sParser ) );

	while( sParser.m_szLine[ 0 ] )
		if( !sParser.GetLine( ) )
			return false;

	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenTerm( SParserOBO& sParser ) {
	bool	fRet, fHit;
	SNode	sNode;
	size_t	i;

	sParser.Reset( );
	while( sParser.m_vecvecstrParents.size( ) < ( sParser.m_vecNodes.size( ) + 1 ) )
		sParser.m_vecvecstrParents.push_back( vector<string>( ) );
	sParser.m_vecvecstrParents[ sParser.m_vecNodes.size( ) ].clear( );
	while( sParser.m_szLine[ 0 ] ) {
		fRet = fHit = false;
		switch( sParser.m_szLine[ 0 ] ) {
			case 'a':
				if( sParser.IsStart( c_szAltID ) ) {
					fHit = true;
					if( !( fRet = OpenAltID( sParser ) ) )
						g_CatSleipnir( ).error( "COntologyOBOImpl::OpenTerm( ) failed: %s", c_szAltID ); }
				break;

			case 'i':
				if( sParser.IsStart( c_szID ) ) {
					fHit = true;
					if( !( fRet = OpenID( sParser ) ) )
						g_CatSleipnir( ).error( "COntologyOBOImpl::OpenTerm( ) failed: %s", c_szID ); }
				else if( sParser.IsStart( c_szIsA ) ) {
					fHit = true;
					if( !( fRet = OpenParent( sParser ) ) )
						g_CatSleipnir( ).error( "COntologyOBOImpl::OpenTerm( ) failed: %s", c_szIsA ); }
				else if( sParser.IsStart( c_szIsObsolete ) ) {
					fHit = true;
					if( !( fRet = OpenObsolete( sParser ) ) )
						g_CatSleipnir( ).error( "COntologyOBOImpl::OpenTerm( ) failed: %s",
							c_szIsObsolete ); }
				break;

			case 'n':
				if( sParser.IsStart( c_szName ) ) {
					fHit = true;
					if( !( fRet = OpenName( sParser ) ) )
						g_CatSleipnir( ).error( "COntologyOBOImpl::OpenTerm( ) failed: %s", c_szName ); }
				else if( sParser.IsStart( c_szNamespace ) ) {
					fHit = true;
					if( !( fRet = OpenNamespace( sParser ) ) )
						g_CatSleipnir( ).error( "COntologyOBOImpl::OpenTerm( ) failed: %s",
							c_szNamespace ); }
				break;

			case 'r':
				if( sParser.IsStart( c_szRelationship ) ) {
					fHit = true;
					if( !( fRet = OpenRelationship( sParser ) ) )
						g_CatSleipnir( ).error( "COntologyOBOImpl::OpenTerm( ) failed: %s",
							c_szRelationship ); }
				break; }
		if( !fHit ) {
			g_CatSleipnir( ).info( "COntologyOBOImpl::OpenTerm( ) skipping: %s", sParser.m_szLine );
			fRet = sParser.GetLine( ); }
		if( !fRet ) {
			g_CatSleipnir( ).error( "COntologyOBOImpl::OpenTerm( ) failed: %s", sParser.m_szLine );
			return false; } }

	if( !sParser.m_fObsolete && ( sParser.m_strNamespace == sParser.m_szTarget ) ) {
		sNode.m_strGloss = sParser.m_strGloss;
		sNode.m_strID = sParser.m_vecstrIDs[ 0 ];
		m_mapNodes[ sNode.m_strID ] = sParser.m_vecNodes.size( );
		for( i = 1; i < sParser.m_vecstrIDs.size( ); ++i )
			m_mapNodes[ sParser.m_vecstrIDs[ i ] ] = sParser.m_vecNodes.size( );
		sParser.m_vecNodes.push_back( sNode ); }

	return true; }

bool COntologyOBOImpl::OpenID( SParserOBO& sParser ) {

	sParser.m_vecstrIDs.push_back( sParser.m_szLine + strlen( c_szID ) );
	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenAltID( SParserOBO& sParser ) {

	sParser.m_vecstrIDs.push_back( sParser.m_szLine + strlen( c_szAltID ) );
	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenName( SParserOBO& sParser ) {

	sParser.m_strGloss = sParser.m_szLine + strlen( c_szName );
	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenParent( SParserOBO& sParser ) {
	const char*	szStart;
	const char*	szEnd;

	szStart = sParser.m_szLine + strlen( c_szIsA );
	for( szEnd = szStart; *szEnd && !isspace( *szEnd ); ++szEnd );
	sParser.m_vecvecstrParents[ sParser.m_vecNodes.size( ) ].push_back( string( szStart, szEnd ) );
	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenRelationship( SParserOBO& sParser ) {
	const char*	szStart;
	const char*	szEnd;

	if( strncmp( sParser.m_szLine + strlen( c_szRelationship ), c_szPartOf,
		strlen( c_szPartOf ) ) ) {
		g_CatSleipnir( ).info( "COntologyOBOImpl::OpenRelationship( %s ) unknown relationship",
			sParser.m_szLine );
		return sParser.GetLine( ); }

	szStart = sParser.m_szLine + strlen( c_szRelationship ) + strlen( c_szPartOf );
	for( szEnd = szStart; *szEnd && !isspace( *szEnd ); ++szEnd );
	sParser.m_vecvecstrParents[ sParser.m_vecNodes.size( ) ].push_back( string( szStart, szEnd ) );
	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenNamespace( SParserOBO& sParser ) {

	sParser.m_strNamespace = sParser.m_szLine + strlen( c_szNamespace );
	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenObsolete( SParserOBO& sParser ) {

	sParser.m_fObsolete = true;
	return sParser.GetLine( ); }

bool COntologyOBOImpl::OpenGenes( SParserOBO& sParser ) {
	size_t									i, j;
	SParserOBO::TSetPGene::const_iterator	iterGene;

	g_CatSleipnir( ).info( "COntologyOBOImpl::OpenGenes( )" );
	if( !sParser.GetLine( ) )
		return false;
	if( !sParser.m_szLine[ 0 ] )
		return true;

	sParser.m_vecsetpGenes.resize( m_iNodes );
	while( sParser.m_istm.peek( ) != EOF )
		if( !OpenGene( sParser ) )
			return false;
	if( !OpenGene( sParser ) )
		return false;

	for( i = 0; i < m_iNodes; ++i ) {
		if( sParser.m_vecsetpGenes[ i ].empty( ) )
			continue;
		m_aNodes[ i ].m_apGenes = new const CGene*[ m_aNodes[ i ].m_iGenes =
			sParser.m_vecsetpGenes[ i ].size( ) ];
		for( j = 0,iterGene = sParser.m_vecsetpGenes[ i ].begin( );
			iterGene != sParser.m_vecsetpGenes[ i ].end( ); ++j,++iterGene )
			m_aNodes[ i ].m_apGenes[ j ] = *iterGene; }

	return true; }

bool COntologyOBOImpl::OpenGene( SParserOBO& sParser ) {
	size_t						i;
	string						strID, strName;
	vector<string>				vecstrLine, vecstrNames;
	TMapStrI::const_iterator	iterNode;

	if( sParser.m_szLine[ 0 ] == '!' )
		return sParser.GetLine( );

//	1	DB ID
//	2	Name
//	3	NOT
//	4	GO ID
//	6	Annotation source
//	9	Gloss
//	10	Syns
	CMeta::Tokenize( sParser.m_szLine, vecstrLine );
	if( ( vecstrLine.size( ) < 11 ) || !( strID = vecstrLine[ 4 ] ).length( ) )
		return false;
	if( vecstrLine[ 3 ].length( ) || ( ( iterNode = m_mapNodes.find( strID ) ) ==
		m_mapNodes.end( ) ) )
		return sParser.GetLine( );
	CMeta::Tokenize( vecstrLine[ 10 ].c_str( ), vecstrNames, "|" );

	while( !vecstrNames.empty( ) && vecstrNames[ 0 ].empty( ) )
		vecstrNames.erase( vecstrNames.begin( ) );
	strName = ( sParser.m_fSynonyms || vecstrNames.empty( ) ) ? vecstrLine[ 2 ] : vecstrNames[ 0 ];	
	strName = ( sParser.m_fDBIDs ? vecstrLine[ 1 ] : strName );

	if( strName.empty( ) ) {
		g_CatSleipnir( ).error( "COntologyOBOImpl::OpenGene( ) null name: %s",
			sParser.m_szLine );
		return false; }
	{
		CGene&	Gene	= sParser.m_Genome.AddGene( strName );

		if( sParser.m_fSynonyms )
			sParser.m_Genome.AddSynonym( Gene, vecstrLine[ 2 ] );
		//if( sParser.m_fDBIDs )
		//	sParser.m_Genome.AddSynonym( Gene, vecstrLine[ 1 ] );
		if( vecstrLine[ 2 ].length( ) ) {
			strID = ( ( i = vecstrLine[ 2 ].find( c_szHUMAN ) ) == string::npos ) ? vecstrLine[ 2 ] :
				vecstrLine[ 2 ].substr( 0, i );				
			sParser.m_Genome.AddSynonym( Gene, strID ); }
		for( i = 1; i < vecstrNames.size( ); ++i )
			sParser.m_Genome.AddSynonym( Gene, vecstrNames[ i ] );
		Gene.AddAnnotation( m_pOntology, iterNode->second );
		if( Gene.GetGloss( ).length( ) == 0 )
			Gene.SetGloss( vecstrLine[ 9 ] );
		sParser.m_vecsetpGenes[ iterNode->second ].insert( &Gene );
	}

	return sParser.GetLine( ); }

}
