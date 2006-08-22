#include "stdafx.h"
#include "annotation.h"
#include "genome.h"

namespace libBioUtils {

const char	COntologyMIPSImpl::c_szMIPS[]	= "MIPS";

COntologyMIPS::COntologyMIPS( ) {

	m_pOntology = this; }

COntologyMIPSImpl::SParserMIPS::SParserMIPS( istream& istm, CGenome& Genome ) :
	SParser( istm, Genome ) { }

COntologyMIPSImpl::COntologyMIPSImpl( ) : COntologyImpl( c_szMIPS ) { }

size_t COntologyMIPS::GetNode( const string& strID ) const {

	return COntologyImpl::GetNode( strID ); }

bool COntologyMIPS::IsAnnotated( size_t iNode, const CGene& Gene, bool fKids ) const {

	return COntologyImpl::IsAnnotated( iNode, Gene, fKids ); }

size_t COntologyMIPS::GetNodes( ) const {

	return COntologyImpl::GetNodes( ); }

const string& COntologyMIPS::GetID( ) const {

	return COntologyImpl::GetID( ); }

const string& COntologyMIPS::GetID( size_t iNode ) const {

	return COntologyImpl::GetID( iNode ); }

const string& COntologyMIPS::GetGloss( size_t iNode ) const {

	return COntologyImpl::GetGloss( iNode ); }

size_t COntologyMIPS::GetParents( size_t iNode ) const {

	return COntologyImpl::GetParents( iNode ); }

size_t COntologyMIPS::GetParent( size_t iNode, size_t iParent ) const {

	return COntologyImpl::GetParent( iNode, iParent ); }

size_t COntologyMIPS::GetChildren( size_t iNode ) const {

	return COntologyImpl::GetChildren( iNode ); }

size_t COntologyMIPS::GetChild( size_t iNode, size_t iChild ) const {

	return COntologyImpl::GetChild( iNode, iChild ); }

size_t COntologyMIPS::GetGenes( size_t iNode, bool fKids ) const {

	return COntologyImpl::GetGenes( iNode, fKids ); }

const CGene& COntologyMIPS::GetGene( size_t iNode, size_t iGene ) const {

	return COntologyImpl::GetGene( iNode, iGene ); }

bool COntologyMIPS::Open( istream& istmOnto, istream& istmGene, CGenome& Genome ) {
	SParserMIPS	sParserOnto( istmOnto, Genome );
	SParserMIPS	sParserGene( istmGene, Genome );

	return ( OpenOntology( sParserOnto ) && OpenGenes( sParserGene ) ); }

bool COntologyMIPSImpl::OpenOntology( SParserMIPS& sParser ) {
	size_t					i, j;
	vector<vector<size_t> >	vecveciChildren;

	g_CatBioUtils.info( "COntologyMIPSImpl::OpenOntology( )" );
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

	g_CatBioUtils.info( "COntologyMIPSImpl::OpenGenes( )" );
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

void COntologyMIPS::GetGeneNames( vector<string>& vecstrGenes ) const {

	return COntologyImpl::GetGeneNames( vecstrGenes ); }

void COntologyMIPS::TermFinder( const CGenes& Genes, vector<STermFound>& vecsTerms,
	bool fBon, bool fKids, bool fBack, const CGenes* pBkg ) const {

	return COntologyImpl::TermFinder( Genes, vecsTerms, fBon, fKids, fBack, pBkg ); }

const char	COntologyMIPSPhenotypes::c_szMIPSPhen[]	= "MIPSP";

COntologyMIPSPhenotypes::COntologyMIPSPhenotypes( ) {

	m_pOntology = this;
	m_strID = c_szMIPSPhen; }

}
