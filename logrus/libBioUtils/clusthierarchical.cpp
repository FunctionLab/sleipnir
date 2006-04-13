#include "stdafx.h"
#include "clusthierarchical.h"
#include "dat.h"
#include "measure.h"

namespace libBioUtils {

CHierarchy::CHierarchy( size_t iID, float dScore, const CHierarchy* pLeft,
	const CHierarchy* pRight ) {

	assert( ( pLeft && pRight ) || !( pLeft || pRight ) );

	m_iID = iID;
	m_dScore = dScore;
	m_pLeft = pLeft;
	m_pRight = pRight;
	m_iWeight = ( m_pLeft && m_pRight ) ? ( m_pLeft->m_iWeight + m_pRight->m_iWeight ) : 1; }

CHierarchyImpl::~CHierarchyImpl( ) {

	if( m_pLeft )
		delete m_pLeft;
	if( m_pRight )
		delete m_pRight; }

void CHierarchy::Destroy( ) {

	delete this; }

float CHierarchy::GetSimilarity( ) const {

	return m_dScore; }

bool CHierarchy::IsGene( ) const {

	return CHierarchyImpl::IsGene( ); }

bool CHierarchyImpl::IsGene( ) const {

	return !( m_pLeft && m_pRight ); }

size_t CHierarchy::GetID( ) const {

	return m_iID; }

const CHierarchy& CHierarchy::Get( bool fRight ) const {

	return *( fRight ? m_pRight : m_pLeft ); }

void CHierarchy::Save( ostream& ostm, size_t iGenes ) const {
	size_t	i;

	for( i = 0; ( i + 1 ) < iGenes; ++i )
		CHierarchyImpl::Save( ostm, i ); }

bool CHierarchyImpl::Save( ostream& ostm, size_t iNode ) const {

	if( IsGene( ) )
		return false;

	if( iNode == m_iID ) {
		ostm << GetSave( ) << '\t' << m_pLeft->GetSave( ) << '\t' <<
			m_pRight->GetSave( ) << '\t' << m_dScore << endl;
		return true; }

	return ( ((const CHierarchyImpl*)m_pLeft)->Save( ostm, iNode ) ||
		((const CHierarchyImpl*)m_pRight)->Save( ostm, iNode ) ); }

string CHierarchyImpl::GetSave( ) const {
	string	strRet;
	char	achBuf[ 16 ];

	strRet = IsGene( ) ? "GENE" : "NODE";
	sprintf_s( achBuf, "%d", m_iID );
	strRet += achBuf;

	return strRet; }

void CHierarchy::GetGenes( vector<size_t>& veciGenes ) const {

	if( IsGene( ) )
		veciGenes.push_back( m_iID );
	else {
		m_pLeft->GetGenes( veciGenes );
		m_pRight->GetGenes( veciGenes ); } }

float CHierarchy::SortChildren( const vector<float>& vecdPCL ) {
	float				dLeft, dRight, dRet;
	const CHierarchy*	pTemp;

	if( IsGene( ) )
		return vecdPCL[ m_iID ];

	dLeft = ((CHierarchy*)m_pLeft)->SortChildren( vecdPCL );
	dRight = ((CHierarchy*)m_pRight)->SortChildren( vecdPCL );
	dRet = ( ( dLeft * m_pLeft->m_iWeight ) + ( dRight * m_pRight->m_iWeight ) ) /
		( m_pRight->m_iWeight + m_pLeft->m_iWeight );
	if( dLeft < dRight ) {
		pTemp = m_pLeft;
		m_pLeft = m_pRight;
		m_pRight = pTemp; }

	return dRet; }

// Implementation courtesy of TIGR MeV

CHierarchy* CClustHierarchical::Cluster( const CDistanceMatrix& Dist ) {
	CDistanceMatrix	Sim;
	size_t			i, j, k, m, iP, iAssigned, iParentless;
	float			d, dTotal, dMin;
	vector<float>	vecdHeight, vecdMax;
	vector<size_t>	veciChild1, veciChild2, veciChildren, veciMax, veciOwner;

	Sim.Initialize( Dist.GetSize( ) );
	for( i = 0; i < Sim.GetSize( ); ++i )
		Sim.Set( i, Dist.Get( i ) );
	iAssigned = iParentless = Sim.GetSize( );
	dTotal = FLT_MAX;

	vecdHeight.resize( Sim.GetSize( ) );
	veciChild1.resize( Sim.GetSize( ) );
	veciChild2.resize( Sim.GetSize( ) );
	veciChildren.resize( Sim.GetSize( ) * 2 );
	for( i = 0; i < veciChild1.size( ); ++i ) {
		veciChild1[ i ] = veciChild2[ i ] = -1;
		veciChildren[ i ] = 1; }

	dMin = FLT_MAX;
	vecdMax.resize( Sim.GetSize( ) );
	veciMax.resize( Sim.GetSize( ) );
	for( i = 0; ( i + 1 ) < Sim.GetSize( ); ++i ) {
		vecdMax[ i ] = -FLT_MAX;
		for( j = ( i + 1 ); j < Sim.GetSize( ); ++j ) {
			if( ( d = Sim.Get( i, j ) ) > vecdMax[ i ] ) {
				vecdMax[ i ] = d;
				veciMax[ i ] = j; }
			if( d < dMin )
				dMin = d; } }

	veciOwner.resize( Sim.GetSize( ) );
	for( i = 0; i < veciOwner.size( ); ++i )
		veciOwner[ i ] = i;
	while( iParentless > 1 ) {
		float	dHeight;
		size_t	iOne, iTwo;

		if( !( iParentless % 500 ) )
			g_CatBioUtils.notice( "CClustHierarchical::Cluster( ) %d/%d nodes remaining", iParentless,
				Dist.GetSize( ) );
		dHeight = -FLT_MAX;
		for( i = 0; i < Sim.GetSize( ); ++i )
			if( ( veciOwner[ i ] != -1 ) && ( vecdMax[ i ] > dHeight ) ) {
				dHeight = vecdMax[ i ];
				j = i;
				iOne = veciMax[ i ]; }
		i = iOne;

		if( ( vecdHeight[ ( k = iAssigned++ ) - Sim.GetSize( ) ] = dHeight ) < dTotal )
			dTotal = dHeight;
		AssertParentage( veciChildren, veciChild1, veciChild2, veciOwner[ i ], k );
		AssertParentage( veciChildren, veciChild1, veciChild2, veciOwner[ j ], k );
		iParentless--;

		iOne = veciChildren[ veciOwner[ i ] ];
		iTwo = veciChildren[ veciOwner[ j ] ];
		veciOwner[ i ] = k;
		veciOwner[ j ] = -1;
		for( iP = 0; iP < Sim.GetSize( ); ++iP )
			if( ( iP != i ) && ( veciOwner[ iP ] != -1 ) )
				Sim.Set( i, iP, ( ( Sim.Get( i, iP ) * iOne ) + ( Sim.Get( j, iP ) * iTwo ) ) /
					( iOne + iTwo ) );

		veciMax[ i ] = i;
		for( iP = 0; iP < Sim.GetSize( ); ++iP )
			if( ( veciOwner[ iP ] != -1 ) && ( ( veciMax[ iP ] == i ) || ( veciMax[ iP ] == j ) ) ) {
				if( ( vecdMax[ iP ] == dMin ) && ( iP != i ) ) {
					veciMax[ iP ] = i;
					continue; }
				vecdMax[ iP ] = -FLT_MAX;
				for( m = ( iP + 1 ); m < Sim.GetSize( ); ++m )
					if( ( veciOwner[ m ] != -1 ) && ( ( d = Sim.Get( iP, m ) ) > vecdMax[ iP ] ) ) {
						vecdMax[ iP ] = d;
						veciMax[ iP ] = m; } } }

	return ConstructHierarchy( veciChild1, veciChild2, vecdHeight, ( 2 * Sim.GetSize( ) ) - 2 ); }

void CClustHierarchicalImpl::AssertParentage( vector<size_t>& veciChildren, vector<size_t>& veciChild1,
	vector<size_t>& veciChild2, size_t iChild, size_t iParent ) {

	veciChildren[ iParent ] += veciChildren[ iChild ];
	iParent -= veciChild1.size( );
	veciChild2[ iParent ] = veciChild1[ iParent ];
	veciChild1[ iParent ] = iChild; }

CHierarchy* CClustHierarchicalImpl::ConstructHierarchy( const vector<size_t>& veciChild1,
	const vector<size_t>& veciChild2, const vector<float>& vecdHeight, size_t iID ) {
	bool	fNode;

	if( fNode = ( iID >= veciChild1.size( ) ) )
		iID -= veciChild1.size( );
	return new CHierarchy( iID, vecdHeight[ iID ],
			fNode ? ConstructHierarchy( veciChild1, veciChild2, vecdHeight, veciChild1[ iID ] ) : NULL,
			fNode ? ConstructHierarchy( veciChild1, veciChild2, vecdHeight, veciChild2[ iID ] ) : NULL ); }

}
