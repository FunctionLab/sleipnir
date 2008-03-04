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

CHierarchy* CClustHierarchical::Cluster( const CDistanceMatrix& Dist, const vector<bool>& vecfGenes ) {

	return CClustHierarchicalImpl::Cluster( Dist, &vecfGenes ); }

CHierarchy* CClustHierarchical::Cluster( const CDistanceMatrix& Dist ) {

	return CClustHierarchicalImpl::Cluster( Dist ); }

// Implementation inspired by TIGR MeV

CHierarchy* CClustHierarchicalImpl::Cluster( const CDistanceMatrix& Dist, const vector<bool>* pvecfGenes ) {
	CDistanceMatrix	Sim;
	size_t			i, j, k, iAssigned, iParentless, iOne, iTwo;
	float			d, dTotal, dMin;
	vector<float>	vecdHeight, vecdMax;
	vector<size_t>	veciChild1, veciChild2, veciChildren, veciMax, veciOwner;

	if( pvecfGenes ) {
		for( i = j = 0; i < pvecfGenes->size( ); ++i )
			if( (*pvecfGenes)[ i ] )
				j++;
		Sim.Initialize( j );
		for( iOne = i = 0; i < Dist.GetSize( ); ++i )
			if( (*pvecfGenes)[ i ] ) {
				for( iTwo = 0,j = 0; j < Dist.GetSize( ); ++j )
					if( (*pvecfGenes)[ j ] )
						Sim.Set( iOne, iTwo++, Dist.Get( i, j ) );
				iOne++; } }
	else {
		Sim.Initialize( Dist.GetSize( ) );
		for( i = 0; i < Sim.GetSize( ); ++i )
			Sim.Set( i, Dist.Get( i ) ); }
	iAssigned = iParentless = Sim.GetSize( );
	dTotal = FLT_MAX;

	vecdHeight.resize( Sim.GetSize( ) );
	veciChild1.resize( Sim.GetSize( ) );
	veciChild2.resize( Sim.GetSize( ) );
	veciChildren.resize( Sim.GetSize( ) * 2 );
	for( i = 0; i < veciChild1.size( ); ++i ) {
		veciChild1[ i ] = veciChild2[ i ] = -1;
		veciChildren[ i ] = 1; }

	vecdMax.resize( Sim.GetSize( ) );
	veciMax.resize( Sim.GetSize( ) );
	vecdMax[ 0 ] = -FLT_MAX;
	veciMax[ 0 ] = -1;
	for( i = 1; i < Sim.GetSize( ); ++i ) {
		size_t	iMin;

		iMin = 0;
		dMin = Sim.Get( 0, i );
		for( j = 1; j < i; ++j )
			if( ( d = Sim.Get( j, i ) ) > dMin ) {
				dMin = d;
				iMin = j; }
		vecdMax[ i ] = dMin;
		veciMax[ i ] = iMin; }

	veciOwner.resize( Sim.GetSize( ) );
	for( i = 0; i < veciOwner.size( ); ++i )
		veciOwner[ i ] = i;
	while( iParentless > 1 ) {
		float	dHeight;

		if( !( iParentless % 500 ) )
			g_CatBioUtils.notice( "CClustHierarchical::Cluster( ) %d/%d nodes remaining", iParentless,
				Sim.GetSize( ) );
		dHeight = -FLT_MAX;
		for( k = 0; k < Sim.GetSize( ); ++k )
			if( vecdMax[ k ] > dHeight )
				dHeight = vecdMax[ i = k ];
		j = veciMax[ i ];

		if( ( vecdHeight[ ( k = iAssigned++ ) - Sim.GetSize( ) ] = dHeight ) < dTotal )
			dTotal = dHeight;
		iParentless--;

		UpdateDistances( i, j, Sim, veciChildren[ veciOwner[ i ] ], veciChildren[ veciOwner[ j ] ], vecdMax,
			veciMax );
		AssertParentage( veciChildren, veciChild1, veciChild2, veciOwner[ i ], k );
		AssertParentage( veciChildren, veciChild1, veciChild2, veciOwner[ j ], k );
		veciOwner[ i ] = k;
		veciOwner[ j ] = -1; }

	return ConstructHierarchy( veciChild1, veciChild2, vecdHeight, ( 2 * Sim.GetSize( ) ) - 2 ); }

void CClustHierarchicalImpl::AssertParentage( vector<size_t>& veciChildren, vector<size_t>& veciChild1,
	vector<size_t>& veciChild2, size_t iChild, size_t iParent ) {

	veciChildren[ iParent ] += veciChildren[ iChild ];
	iParent -= veciChild1.size( );
	veciChild2[ iParent ] = veciChild1[ iParent ];
	veciChild1[ iParent ] = iChild; }

void CClustHierarchicalImpl::UpdateDistances( size_t iOne, size_t iTwo, CDistanceMatrix& Sim, size_t iWOne,
	size_t iWTwo, vector<float>& vecdMax, vector<size_t>& veciMax ) {
	float	d, dOne, dTwo;
	size_t	i, j;

	vecdMax[ iOne ] = -FLT_MAX;
	veciMax[ iOne ] = -1;
	// Update row iOne across
	for( i = 0; i < iOne; ++i )
		if( !CMeta::IsNaN( dOne = Sim.Get( i, iOne ) ) && !CMeta::IsNaN( dTwo = Sim.Get( i, iTwo ) ) ) {
			Sim.Set( i, iOne, d = ( ( iWOne * dOne ) + ( iWTwo * dTwo ) ) / ( iWOne + iWTwo ) );
			if( d > vecdMax[ iOne ] ) {
				vecdMax[ iOne ] = d;
				veciMax[ iOne ] = i; } }
	// Update row iOne down
	for( i = ( iOne + 1 ); i < Sim.GetSize( ); ++i )
		if( !CMeta::IsNaN( dOne = Sim.Get( iOne, i ) ) && !CMeta::IsNaN( dTwo = Sim.Get( iTwo, i ) ) ) {
			Sim.Set( iOne, i, d = ( ( iWOne * dOne ) + ( iWTwo * dTwo ) ) / ( iWOne + iWTwo ) );
			if( veciMax[ i ] == iOne ) {
				vecdMax[ i ] = -FLT_MAX;
				veciMax[ i ] = 0;
				for( j = 0; j < i; ++j )
					if( !CMeta::IsNaN( d = Sim.Get( j, i ) ) && ( d > vecdMax[ i ] ) ) {
						vecdMax[ i ] = d;
						veciMax[ i ] = j; } }
			else if( d > vecdMax[ i ] ) {
				vecdMax[ i ] = d;
				veciMax[ i ] = iOne; } }
	// Delete row iTwo across
	for( i = 0; i < iTwo; ++i )
		Sim.Set( i, iTwo, CMeta::GetNaN( ) );
	vecdMax[ iTwo ] = -FLT_MAX;
	veciMax[ iTwo ] = -1;
	// Delete row iTwo down
	for( i = ( iTwo + 1 ); i < Sim.GetSize( ); ++i ) {
		Sim.Set( iTwo, i, CMeta::GetNaN( ) );
		if( veciMax[ i ] == iTwo ) {
			vecdMax[ i ] = -FLT_MAX;
			veciMax[ i ] = 0;
			for( j = 0; j < i; ++j )
				if( !CMeta::IsNaN( d = Sim.Get( j, i ) ) && ( d > vecdMax[ i ] ) ) {
					vecdMax[ i ] = d;
					veciMax[ i ] = j; } } } }

CHierarchy* CClustHierarchicalImpl::ConstructHierarchy( const vector<size_t>& veciChild1,
	const vector<size_t>& veciChild2, const vector<float>& vecdHeight, size_t iID ) {
	bool	fNode;

	if( fNode = ( iID >= veciChild1.size( ) ) )
		iID -= veciChild1.size( );
	return new CHierarchy( iID, vecdHeight[ iID ],
			fNode ? ConstructHierarchy( veciChild1, veciChild2, vecdHeight, veciChild1[ iID ] ) : NULL,
			fNode ? ConstructHierarchy( veciChild1, veciChild2, vecdHeight, veciChild2[ iID ] ) : NULL ); }

}
