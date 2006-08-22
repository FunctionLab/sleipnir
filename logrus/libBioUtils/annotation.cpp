#include "stdafx.h"
#include "annotation.h"
#include "genome.h"
#include "statistics.h"

namespace libBioUtils {

COntologyImpl::SNode::SNode( ) : m_iParents(0), m_aiParents(NULL), m_iChildren(0),
	m_aiChildren(NULL), m_iGenes(0), m_apGenes(NULL), m_iCacheGenes(-1),
	m_apCacheGenes(NULL) { }

void COntologyImpl::SNode::Reset( ) {

	m_strID.clear( );
	m_iParents = m_iChildren = m_iGenes = 0;
	if( m_aiParents )
		delete[] m_aiParents;
	if( m_aiChildren )
		delete[] m_aiChildren;
	if( m_apGenes )
		delete[] m_apGenes;
	if( m_apCacheGenes )
		delete[] m_apCacheGenes; }

COntologyImpl::SParser::SParser( istream& istm, CGenome& Genome ) : m_istm(istm),
	m_Genome(Genome) {

	m_szLine[ 0 ] = 0; }

bool COntologyImpl::SParser::GetLine( ) {

	m_istm.getline( m_szLine, c_iBuffer - 1 );
	g_CatBioUtils.debug( "COntologyImpl::SParser::GetLine( ) %s", m_szLine );
	return true; }

bool COntologyImpl::SParser::IsStart( const char* szStart ) const {

	return !strncmp( m_szLine, szStart, strlen( szStart ) ); }

COntologyImpl::COntologyImpl( const string& strID ) : m_strID(strID), m_iNodes(0),
	m_aNodes(NULL) { }

COntologyImpl::~COntologyImpl( ) {

	Reset( ); }

size_t COntologyImpl::GetNodes( ) const {

	return m_iNodes; }

size_t COntologyImpl::GetParents( size_t iNode ) const {

	return m_aNodes[ iNode ].m_iParents; }

size_t COntologyImpl::GetParent( size_t iNode, size_t iParent ) const {

	return m_aNodes[ iNode ].m_aiParents[ iParent ]; }

size_t COntologyImpl::GetChildren( size_t iNode ) const {

	return m_aNodes[ iNode ].m_iChildren; }

size_t COntologyImpl::GetChild( size_t iNode, size_t iChild ) const {

	return m_aNodes[ iNode ].m_aiChildren[ iChild ]; }

size_t COntologyImpl::GetGenes( size_t iNode, bool fKids ) const {
	size_t	iRet;

	iRet = m_aNodes[ iNode ].m_iGenes;
	if( fKids ) {
		CollectGenes( iNode );
		iRet += m_aNodes[ iNode ].m_iCacheGenes; }

	return iRet; }

const CGene& COntologyImpl::GetGene( size_t iNode, size_t iGene ) const {
	size_t	i;

	if( iGene < ( i = m_aNodes[ iNode ].m_iGenes ) )
		return *m_aNodes[ iNode ].m_apGenes[ iGene ];

	CollectGenes( iNode );
	return *m_aNodes[ iNode ].m_apCacheGenes[ iGene - i ]; }

const string& COntologyImpl::GetID( ) const {

	return m_strID; }

const string& COntologyImpl::GetID( size_t iNode ) const {

	return m_aNodes[ iNode ].m_strID; }

const string& COntologyImpl::GetGloss( size_t iNode ) const {

	return m_aNodes[ iNode ].m_strGloss; }

void COntologyImpl::Reset( ) {
	size_t	i;

	m_mapNodes.clear( );
	if( !m_aNodes )
		return;

	for( i = 0; i < m_iNodes; ++i )
		m_aNodes[ i ].Reset( );
	m_iNodes = 0;
	delete[] m_aNodes;
	m_aNodes = NULL; }

bool COntologyImpl::IsAnnotated( size_t iNode, const CGene& Gene, bool fKids ) const {
	size_t	i;

	for( i = 0; i < m_aNodes[ iNode ].m_iGenes; ++i )
		if( m_aNodes[ iNode ].m_apGenes[ i ] == &Gene )
			return true;
	if( fKids ) {
		CollectGenes( iNode );
		for( i = 0; i < m_aNodes[ iNode ].m_iCacheGenes; ++i )
			if( m_aNodes[ iNode ].m_apCacheGenes[ i ] == &Gene )
				return true; }

	return false; }

void COntologyImpl::CollectGenes( size_t iNode ) const {
	TSetPGenes	setpGenes;

	if( m_aNodes[ iNode ].m_iCacheGenes == -1 )
		((COntologyImpl*)this)->CollectGenes( iNode, setpGenes ); }

void COntologyImpl::CollectGenes( size_t iNode, TSetPGenes& setpGenes ) {
	size_t						i;
	TSetPGenes					setpKids;
	TSetPGenes::const_iterator	iterGenes;

	if( m_aNodes[ iNode ].m_iCacheGenes != -1 ) {
		for( i = 0; i < m_aNodes[ iNode ].m_iGenes; ++i )
			setpGenes.insert( m_aNodes[ iNode ].m_apGenes[ i ] );
		for( i = 0; i < m_aNodes[ iNode ].m_iCacheGenes; ++i )
			setpGenes.insert( m_aNodes[ iNode ].m_apCacheGenes[ i ] );
		return; }

	for( i = 0; i < m_aNodes[ iNode ].m_iGenes; ++i )
		setpGenes.insert( m_aNodes[ iNode ].m_apGenes[ i ] );
	for( i = 0; i < m_aNodes[ iNode ].m_iChildren; ++i )
		CollectGenes( m_aNodes[ iNode ].m_aiChildren[ i ], setpKids );
	if( m_aNodes[ iNode ].m_iCacheGenes = setpKids.size( ) ) {
		m_aNodes[ iNode ].m_apCacheGenes = new const CGene*[ setpKids.size( ) ];
		for( i = 0,iterGenes = setpKids.begin( ); iterGenes != setpKids.end( );
			++i,++iterGenes ) {
			m_aNodes[ iNode ].m_apCacheGenes[ i ] = *iterGenes;
			setpGenes.insert( *iterGenes ); } } }

size_t COntologyImpl::GetNode( const string& strID ) const {
	TMapStrI::const_iterator	iterNode;

	iterNode = m_mapNodes.find( strID );
	return ( ( iterNode == m_mapNodes.end( ) ) ? -1 : iterNode->second ); }

void COntologyImpl::GetGeneNames( vector<string>& vecstrGenes ) const {
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGene;
	size_t						i, j;

	for( i = 0; i < m_iNodes; ++i )
		for( j = 0; j < m_aNodes[ i ].m_iGenes; ++j )
			setstrGenes.insert( m_aNodes[ i ].m_apGenes[ j ]->GetName( ) );

	for( iterGene = setstrGenes.begin( ); iterGene != setstrGenes.end( ); ++iterGene )
		vecstrGenes.push_back( *iterGene ); }

void CSlimImpl::Reset( const IOntology* pOntology ) {

	m_pOntology = pOntology;
	m_vecstrSlims.clear( );
	m_vecveciTerms.clear( );
	m_vecvecpGenes.clear( ); }

bool CSlim::Open( istream& istm, const IOntology* pOnto ) {
	static const size_t	c_iBuffer	= 1024;
	char								szBuf[ c_iBuffer ];
	size_t								i, j, k, iNode;
	string								str;
	set<const CGene*>					setiGenes;
	set<const CGene*>::const_iterator	iterGene;

	g_CatBioUtils.info( "CSlim::Open( %s )", pOnto->GetID( ).c_str( ) );

	Reset( pOnto );
	while( istm.peek( ) != EOF ) {
		i = m_vecveciTerms.size( );
		m_vecveciTerms.resize( i + 1 );
		istm.getline( szBuf, c_iBuffer - 1 );
		{
			istrstream	issm( szBuf );

			while( issm.peek( ) != EOF ) {
				str = OpenToken( issm );
				if( !str.length( ) )
					break;
				if( m_vecstrSlims.size( ) <= i )
					m_vecstrSlims.push_back( str );
				else {
					if( ( j = m_pOntology->GetNode( str ) ) == -1 ) {
						g_CatBioUtils.error( "CSlim::Open( %s ) unknown node: %s",
							pOnto->GetID( ).c_str( ), str.c_str( ) );
						return false; }
					m_vecveciTerms[ i ].push_back( j ); } }
		} }

	m_vecvecpGenes.resize( m_vecveciTerms.size( ) );
	for( i = 0; i < m_vecveciTerms.size( ); ++i ) {
		setiGenes.clear( );
		for( j = 0; j < m_vecveciTerms[ i ].size( ); ++j ) {
			iNode = m_vecveciTerms[ i ][ j ];
			for( k = 0; k < m_pOntology->GetGenes( iNode, true ); ++k )
				setiGenes.insert( &m_pOntology->GetGene( iNode, k ) ); }
		for( iterGene = setiGenes.begin( ); iterGene != setiGenes.end( ); ++iterGene )
			m_vecvecpGenes[ i ].push_back( *iterGene ); }

	return true; }

const CGene& CSlim::GetGene( size_t iSlim, size_t iGene ) const {

	return *m_vecvecpGenes[ iSlim ][ iGene ]; }

size_t CSlim::GetSlims( ) const {

	return m_vecstrSlims.size( ); }

size_t CSlim::GetGenes( size_t iSlim ) const {

	return m_vecvecpGenes[ iSlim ].size( ); }

void CSlim::GetGeneNames( vector<string>& vecstrGenes ) const {
	set<const CGene*>					setpGenes;
	set<const CGene*>::const_iterator	iterGene;
	size_t								i, j;

	for( i = 0; i < m_vecvecpGenes.size( ); ++i )
		for( j = 0; j < m_vecvecpGenes[ i ].size( ); ++j )
			setpGenes.insert( m_vecvecpGenes[ i ][ j ] );

	for( iterGene = setpGenes.begin( ); iterGene != setpGenes.end( ); ++iterGene )
		vecstrGenes.push_back( (*iterGene)->GetName( ) ); }

const string& CSlim::GetSlim( size_t iSlim ) const {

	return m_vecstrSlims[ iSlim ]; }

void COntologyImpl::TermFinder( const CGenes& Genes, vector<STermFound>& vecsTerms, bool fBon,
	bool fKids, bool fBack, const CGenes* pBkg ) const {
	size_t			i, j, iMult, iBkg, iGenes;
	double			d;
	vector<size_t>	veciAnno;

	iBkg = pBkg ? pBkg->GetGenes( ) :
		( fBack ? Genes.GetGenome( ).GetGenes( ) : Genes.GetGenome( ).CountGenes( m_pOntology ) );
	veciAnno.resize( m_iNodes );
	for( iMult = i = 0; i < veciAnno.size( ); ++i )
		if( veciAnno[ i ] = Genes.CountAnnotations( m_pOntology, i, fKids, pBkg ) )
			iMult++;
	for( i = 0; i < m_iNodes; ++i ) {
		if( !veciAnno[ i ] )
			continue;
		if( pBkg ) {
			for( iGenes = j = 0; j < pBkg->GetGenes( ); ++j )
				if( IsAnnotated( i, pBkg->GetGene( j ), fKids ) )
					iGenes++; }
		else
			iGenes = GetGenes( i, fKids );
		d = CStatistics::HypergeometricCDF( veciAnno[ i ], Genes.GetGenes( ),
			iGenes, iBkg );
/*
//if( m_aNodes[ i ].m_strID == "GO:0050875" )
cerr << veciAnno[ i ] << ", " << Genes.GetGenes( ) << ", " << iGenes << ", " << iBkg << ", " << d << ", " << iMult << endl;
//*/
		if( fBon ) {
			if( ( d *= iMult ) > 1 )
				d = 1; }
		vecsTerms.push_back( STermFound( i, d, veciAnno[ i ], Genes.GetGenes( ), iGenes, iBkg ) ); } }

}
