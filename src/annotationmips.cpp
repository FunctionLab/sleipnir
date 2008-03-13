#include "stdafx.h"
#include "annotation.h"
#include "genome.h"

namespace Sleipnir {

const char	COntologyMIPSImpl::c_szMIPS[]	= "MIPS";

COntologyMIPS::COntologyMIPS( ) {

	m_pOntology = this; }

COntologyMIPSImpl::SParserMIPS::SParserMIPS( istream& istm, CGenome& Genome ) :
	SParser( istm, Genome ) { }

COntologyMIPSImpl::COntologyMIPSImpl( ) : COntologyImpl( c_szMIPS ) { }

/*!
 * \brief
 * Initializes the ontology using the "funcat" scheme (structure) and annotation files.
 * 
 * \param istmOntology
 * Stream from which the catalog structure file is read.
 * 
 * \param istmAnnotations
 * Stream from which the catalog annotation file is read.
 * 
 * \param Genome
 * Genome into which genes are inserted or read during annotation parsing.
 * 
 * \returns
 * True if the ontology was successfully initialized.
 * 
 * COntologyMIPS::Open parses the structure of the MIPS functional catalog from the (organism-independent)
 * funcat scheme file and a set of (organism-dependent) annotations from a data file.  Terms are identified
 * by MIPS IDs (e.g. "16.03.03"), and genes are identified by the name given by the annotation file.  Genes
 * are retrieved from Genome if already present or inserted if not; it is thus important to ensure that the
 * proper primary gene names are used so as to agree with any identifiers already present in Genome.
 */
bool COntologyMIPS::Open( istream& istmOnto, istream& istmGene, CGenome& Genome ) {
	SParserMIPS	sParserOnto( istmOnto, Genome );
	SParserMIPS	sParserGene( istmGene, Genome );

	if( !OpenOntology( sParserOnto ) ) {
		g_CatSleipnir.error( "COntologyMIPS::Open( ) failed on ontology line %d: %s", sParserOnto.m_iLine,
			sParserOnto.m_szLine );
		return false; }
	if( !OpenGenes( sParserGene ) ) {
		g_CatSleipnir.error( "COntologyMIPS::Open( ) failed on genes line %d: %s", sParserGene.m_iLine,
			sParserGene.m_szLine );
		return false; }

	return true; }

bool COntologyMIPSImpl::OpenOntology( SParserMIPS& sParser ) {
	size_t					i, j;
	vector<vector<size_t> >	vecveciChildren;

	g_CatSleipnir.info( "COntologyMIPSImpl::OpenOntology( )" );
	if( !( sParser.GetLine( ) && ( sParser.m_szLine[ 0 ] == '#' ) && sParser.GetLine( ) ) )
		return false;

	while( sParser.m_istm.peek( ) != EOF )
		if( !OpenCategory( sParser ) )
			return false;
	if( !OpenCategory( sParser ) )
		return false;

	m_aNodes = new SNode[ m_iNodes = sParser.m_veciParents.size( ) ];
	vecveciChildren.resize( m_iNodes );
	for( i = 0; i < m_iNodes; ++i ) {
		m_aNodes[ i ].m_strID = sParser.m_vecstrIDs[ i ];
		m_mapNodes[ m_aNodes[ i ].m_strID ] = i;
		m_aNodes[ i ].m_strGloss = sParser.m_vecstrGlosses[ i ];
		if( sParser.m_veciParents[ i ] != -1 ) {
			m_aNodes[ i ].m_aiParents = new size_t[ m_aNodes[ i ].m_iParents = 1 ];
			m_aNodes[ i ].m_aiParents[ 0 ] = sParser.m_veciParents[ i ];
			vecveciChildren[ sParser.m_veciParents[ i ] ].push_back( i ); } }
	for( i = 0; i < m_iNodes; ++i ) {
		if( !vecveciChildren[ i ].size( ) )
			continue;
		m_aNodes[ i ].m_aiChildren = new size_t[ m_aNodes[ i ].m_iChildren =
			vecveciChildren[ i ].size( ) ];
		for( j = 0; j < m_aNodes[ i ].m_iChildren; ++j )
			m_aNodes[ i ].m_aiChildren[ j ] = vecveciChildren[ i ][ j ]; }

	return true; }

bool COntologyMIPSImpl::OpenCategory( SParserMIPS& sParser ) {
	char*	pch;
	size_t	i, iDepth;

	if( !( pch = strchr( sParser.m_szLine, ' ' ) ) )
		return false;

	*(pch++) = 0;
	sParser.m_vecstrIDs.push_back( sParser.m_szLine );
	while( *pch && isspace( *pch ) )
		pch++;
	sParser.m_vecstrGlosses.push_back( pch );
	if( ( iDepth = OpenID( sParser ) ) == -1 )
		return false;
	while( iDepth < sParser.m_stakiHier.size( ) )
		sParser.m_stakiHier.pop( );
	i = sParser.m_veciParents.size( );
	sParser.m_veciParents.push_back( sParser.m_stakiHier.empty( ) ? -1 :
		sParser.m_stakiHier.top( ) );
	if( iDepth >= sParser.m_stakiHier.size( ) )
		sParser.m_stakiHier.push( i );

	return sParser.GetLine( ); }

size_t COntologyMIPSImpl::OpenID( SParserMIPS& sParser ) {
	size_t	iRet;
	char*	pch;

	for( iRet = 0,pch = strchr( sParser.m_szLine, '.' ); pch; ++iRet,
		pch = strchr( ++pch, '.' ) );

	return iRet; }

bool COntologyMIPSImpl::OpenGenes( SParserMIPS& sParser ) {
	size_t	i, j;

	g_CatSleipnir.info( "COntologyMIPSImpl::OpenGenes( )" );
	if( !sParser.GetLine( ) )
		return false;
	if( !sParser.m_szLine[ 0 ] )
		return true;

	sParser.m_vecpGenes.resize( m_iNodes );
	while( sParser.m_istm.peek( ) != EOF )
		if( !OpenGene( sParser ) )
			return false;
	if( !OpenGene( sParser ) )
		return false;

	for( i = 0; i < m_iNodes; ++i ) {
		if( !sParser.m_vecpGenes[ i ].size( ) )
			continue;
		m_aNodes[ i ].m_apGenes = new const CGene*[ m_aNodes[ i ].m_iGenes =
			sParser.m_vecpGenes[ i ].size( ) ];
		for( j = 0; j < m_aNodes[ i ].m_iGenes; ++j )
			m_aNodes[ i ].m_apGenes[ j ] = sParser.m_vecpGenes[ i ][ j ]; }

	return true; }

bool COntologyMIPSImpl::OpenGene( SParserMIPS& sParser ) {
	char*	pchOne;
	char*	pchTwo;
	size_t	iNode;

	if( !( ( pchOne = strchr( sParser.m_szLine, '|' ) ) &&
		( pchTwo = strchr( pchOne + 1, '|' ) ) ) )
		return false;
	*(pchOne++) = *pchTwo = 0;

	iNode = m_mapNodes[ pchOne ];
	{
		CGene&	Gene	= sParser.m_Genome.AddGene( sParser.m_szLine );

		Gene.AddAnnotation( m_pOntology, iNode );
		sParser.m_vecpGenes[ iNode ].push_back( &Gene );
	}

	return sParser.GetLine( ); }

const char	COntologyMIPSPhenotypes::c_szMIPSPhen[]	= "MIPSP";

COntologyMIPSPhenotypes::COntologyMIPSPhenotypes( ) {

	m_pOntology = this;
	m_strID = c_szMIPSPhen; }

}
