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
#include "statistics.h"

namespace Sleipnir {

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

COntologyImpl::SParser::SParser( std::istream& istm, CGenome& Genome ) : m_istm(istm),
	m_Genome(Genome), m_iLine(0) {

	m_szLine[ 0 ] = 0; }

bool COntologyImpl::SParser::GetLine( ) {

	m_iLine++;
	m_istm.getline( m_szLine, c_iBuffer - 1 );
	g_CatSleipnir( ).debug( "COntologyImpl::SParser::GetLine( ) %s", m_szLine );
	return true; }

bool COntologyImpl::SParser::IsStart( const char* szStart ) const {

	return !strncmp( m_szLine, szStart, strlen( szStart ) ); }

const CGene& COntologyImpl::GetGene( size_t iNode, size_t iGene ) const {
	size_t	i;

	if( iGene < ( i = m_aNodes[ iNode ].m_iGenes ) )
		return *m_aNodes[ iNode ].m_apGenes[ iGene ];

	CollectGenes( iNode );
	return *m_aNodes[ iNode ].m_apCacheGenes[ iGene - i ]; }

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

void COntologyImpl::CollectGenes( size_t iNode, TSetPGenes& setpGenes ) {
	SNode&						sNode	= m_aNodes[ iNode ];
	size_t						i;
	TSetPGenes					setpSelf, setpKids;
	TSetPGenes::const_iterator	iterGenes;

	if( sNode.m_iCacheGenes != -1 ) {
		for( i = 0; i < sNode.m_iGenes; ++i )
			setpGenes.insert( sNode.m_apGenes[ i ] );
		for( i = 0; i < sNode.m_iCacheGenes; ++i )
			setpGenes.insert( sNode.m_apCacheGenes[ i ] );
		return; }

	for( i = 0; i < sNode.m_iGenes; ++i ) {
		setpSelf.insert( sNode.m_apGenes[ i ] );
		setpGenes.insert( sNode.m_apGenes[ i ] ); }
	for( i = 0; i < sNode.m_iChildren; ++i )
		CollectGenes( sNode.m_aiChildren[ i ], setpKids );
	sNode.m_iCacheGenes = 0;
	for( iterGenes = setpKids.begin( ); iterGenes != setpKids.end( ); ++iterGenes )
		if( setpSelf.find( *iterGenes ) == setpSelf.end( ) )
			sNode.m_iCacheGenes++;
	if( sNode.m_iCacheGenes ) {
		sNode.m_apCacheGenes = new const CGene*[ sNode.m_iCacheGenes ];
		for( i = 0,iterGenes = setpKids.begin( ); iterGenes != setpKids.end( ); ++iterGenes )
			if( setpSelf.find( *iterGenes ) == setpSelf.end( ) ) {
				sNode.m_apCacheGenes[ i++ ] = *iterGenes;
				setpGenes.insert( *iterGenes ); } } }

size_t COntologyImpl::GetNode( const std::string& strID ) const {
	TMapStrI::const_iterator	iterNode;

	iterNode = m_mapNodes.find( strID );
	return ( ( iterNode == m_mapNodes.end( ) ) ? -1 : iterNode->second ); }

void COntologyImpl::GetGeneNames( std::vector<std::string>& vecstrGenes ) const {
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGene;
	size_t						i, j;

	for( i = 0; i < m_iNodes; ++i )
		for( j = 0; j < m_aNodes[ i ].m_iGenes; ++j )
			setstrGenes.insert( m_aNodes[ i ].m_apGenes[ j ]->GetName( ) );

	for( iterGene = setstrGenes.begin( ); iterGene != setstrGenes.end( ); ++iterGene )
		vecstrGenes.push_back( *iterGene ); }

void COntologyImpl::TermFinder( const CGenes& Genes, vector<STermFound>& vecsTerms, bool fBon,
	bool fKids, bool fBack, float dPValue, const CGenes* pBkg ) const {
	size_t			i, j, iMult, iBkg, iGenes, iGiven;
	double			d;
	vector<size_t>	veciAnno;

	for( iGiven = i = 0; i < Genes.GetGenes( ); ++i )
		if( Genes.GetGene( i ).IsAnnotated( m_pOntology ) )
			iGiven++;
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
		d = CStatistics::HypergeometricCDF( veciAnno[ i ], iGiven, iGenes, iBkg );
		if( fBon ) {
			if( ( d *= iMult ) > 1 )
				d = 1; }
		if( d <= dPValue )
			vecsTerms.push_back( STermFound( i, d, veciAnno[ i ], iGiven, iGenes, iBkg ) ); } }

void CSlimImpl::Reset( const IOntology* pOntology ) {

	m_pOntology = pOntology;
	m_vecstrSlims.clear( );
	m_vecveciTerms.clear( );
	m_vecvecpGenes.clear( ); }

/*!
 * \brief
 * Constructs a slim from a text file listing ontology terms.
 * 
 * \param istmSlim
 * Stream containing term ID strings to include in the slim.
 * 
 * \param pOntology
 * Ontology from which terms are drawn.
 * 
 * \returns
 * True if slim construction succeeded (i.e. all terms were found in the ontology).
 * 
 * Constructs a slim from a tab-delimited text file of the form:
 * \code
 * gloss1	id1
 * gloss2	id2
 * ...
 * glossN	idN
 * \endcode
 * The glosses are human-readable names for ontology terms (e.g. "protein targeting to ER"); these are
 * ignored by the parser.  The IDs are the ontology-specific ID strings (e.g. "GO:0009605") used to look up
 * terms in the given ontology.  Each gloss must be separated from its accompanying ID by a tab, and each
 * line should consist of a single gloss/ID pair.
 */
bool CSlim::Open( std::istream& istmSlim, const IOntology* pOntology ) {
	static const size_t	c_iBuffer	= 1024;
	char								szBuf[ c_iBuffer ];
	size_t								i, j, k, iNode;
	string								str;
	set<const CGene*>					setiGenes;
	set<const CGene*>::const_iterator	iterGene;

	g_CatSleipnir( ).info( "CSlim::Open( %s )", pOntology->GetID( ).c_str( ) );

	Reset( pOntology );
	while( istmSlim.peek( ) != EOF ) {
		i = m_vecveciTerms.size( );
		m_vecveciTerms.resize( i + 1 );
		istmSlim.getline( szBuf, c_iBuffer - 1 );
		{
			istringstream	issm( szBuf );

			while( issm.peek( ) != EOF ) {
				str = OpenToken( issm );
				cout << str << endl;
				if( !str.length( ) )
					break;
				if( m_vecstrSlims.size( ) <= i )
					m_vecstrSlims.push_back( str );
				else {
					if( ( j = m_pOntology->GetNode( str ) ) == -1 ) {
						g_CatSleipnir( ).error( "CSlim::Open( %s ) unknown node: %s",
							m_pOntology->GetID( ).c_str( ), str.c_str( ) );
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

/*!
 * \brief
 * Retrieve the names of all genes annotated below terms in this slim.
 * 
 * \param vecstrGenes
 * Output unique gene names annotated below this slim's ontology terms.
 */
void CSlim::GetGeneNames( std::vector<std::string>& vecstrGenes ) const {
	set<const CGene*>					setpGenes;
	set<const CGene*>::const_iterator	iterGene;
	size_t								i, j;

	for( i = 0; i < m_vecvecpGenes.size( ); ++i )
		for( j = 0; j < m_vecvecpGenes[ i ].size( ); ++j )
			setpGenes.insert( m_vecvecpGenes[ i ][ j ] );

	for( iterGene = setpGenes.begin( ); iterGene != setpGenes.end( ); ++iterGene )
		vecstrGenes.push_back( (*iterGene)->GetName( ) ); }

}
