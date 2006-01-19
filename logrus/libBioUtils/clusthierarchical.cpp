#include "stdafx.h"
#include "clusthierarchical.h"
#include "dat.h"
#include "measure.h"

namespace libBioUtils {

CHierarchy::CHierarchy( size_t iID ) {

	m_iID = iID;
	m_dScore = 0;
	m_iWeight = 1;
	m_pLeft = m_pRight = NULL; }

CHierarchy::CHierarchy( size_t iID, float dScore, const CHierarchy* pLeft,
	const CHierarchy* pRight ) {

	m_iID = iID;
	m_dScore = dScore;
	m_pLeft = pLeft;
	m_pRight = pRight;
	m_iWeight = m_pLeft->GetWeight( ) + m_pRight->GetWeight( ); }

CHierarchyImpl::~CHierarchyImpl( ) {

	if( m_pLeft )
		delete m_pLeft;
	if( m_pRight )
		delete m_pRight; }

void CHierarchy::Destroy( ) {

	delete this; }

bool CHierarchy::IsGene( ) const {

	return CHierarchyImpl::IsGene( ); }

bool CHierarchyImpl::IsGene( ) const {

	return !( m_pLeft && m_pRight ); }

size_t CHierarchy::GetID( ) const {

	return m_iID; }

const CHierarchy& CHierarchy::Get( bool fRight ) const {

	return *( fRight ? m_pRight : m_pLeft ); }

size_t CHierarchy::GetWeight( ) const {

	return m_iWeight; }

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
	dRet = ( ( dLeft * m_pLeft->GetWeight( ) ) + ( dRight * m_pRight->GetWeight( ) ) ) /
		( m_pRight->GetWeight( ) + m_pLeft->GetWeight( ) );
	if( dLeft < dRight ) {
		pTemp = m_pLeft;
		m_pLeft = m_pRight;
		m_pRight = pTemp; }

	return dRet; }

CHierarchy* CClustHierarchical::Cluster( const CDistanceMatrix& Dist ) {
	size_t				i, iNodes;
	vector<SHierarchy>	vecsGenes, vecsNodes;
	TPrHier				prHier;
	TVecVecD			vecvecdDist;
	CDistanceMatrix		DistNodes;
	CHierarchy*			pHier;
	CHierarchy*			pRight;
	float				dScore;

	iNodes = 0;
	DistNodes.Initialize( Dist.GetSize( ) );
	vecvecdDist.resize( Dist.GetSize( ) );
	vecsGenes.resize( Dist.GetSize( ) );
	vecsNodes.resize( Dist.GetSize( ) );
	vecvecdDist.resize( Dist.GetSize( ) );
	for( i = 0; i < Dist.GetSize( ); ++i ) {
		vecsGenes[ i ].m_fUsed = vecsNodes[ i ].m_fUsed = false;
		vecsNodes[ i ].m_pHier = NULL;
		vecsGenes[ i ].m_pHier = new CHierarchy( i ); }

	for( i = 1; i < Dist.GetSize( ); ++i ) {
		g_CatBioUtils.notice( "CClustHierarchical::Cluster( ) processing node %d/%d", i,
			Dist.GetSize( ) );
		FindHighestSimilarity( Dist, vecsGenes, vecsNodes, vecvecdDist, DistNodes, iNodes,
			dScore, pHier, pRight );
		vecsNodes[ iNodes++ ].m_pHier = pHier = new CHierarchy( iNodes, dScore, pHier, pRight );
		CleanDistance( pHier->Get( false ), vecsGenes, vecsNodes, vecvecdDist );
		CleanDistance( pHier->Get( true ), vecsGenes, vecsNodes, vecvecdDist );
		RefreshDistances( Dist, vecsGenes, vecsNodes, vecvecdDist, DistNodes, iNodes - 1 ); }

	return pHier; }

void CClustHierarchicalImpl::FindHighestSimilarity( const CDistanceMatrix& Dist,
	const vector<SHierarchy>& vecsGenes, const vector<SHierarchy>& vecsNodes,
	const TVecVecD& vecvecdDist, const CDistanceMatrix& DistNodes, size_t iNodes,
	float& dMax, CHierarchy*& pLeft, CHierarchy*& pRight ) {
	size_t	i, j;

	dMax = -1;
	for( i = 0; i < Dist.GetSize( ); ++i ) {
		if( vecsGenes[ i ].m_fUsed )
			continue;
		for( j = ( i + 1 ); j < Dist.GetSize( ); ++j )
			if( !vecsGenes[ j ].m_fUsed && ( Dist.Get( i, j ) > dMax ) ) {
				dMax = Dist.Get( i, j );
				pLeft = vecsGenes[ i ].m_pHier;
				pRight = vecsGenes[ j ].m_pHier; } }

	for( i = 0; i < iNodes; ++i ) {
		if( vecsNodes[ i ].m_fUsed )
			continue;
		for( j = 0; j < vecvecdDist[ i ].size( ); ++j )
			if( !vecsGenes[ j ].m_fUsed && ( vecvecdDist[ i ][ j ] > dMax ) ) {
				dMax = vecvecdDist[ i ][ j ];
				pLeft = vecsNodes[ i ].m_pHier;
				pRight = vecsGenes[ j ].m_pHier; }
		for( j = ( i + 1 ); j < iNodes; ++j )
			if( !vecsNodes[ j ].m_fUsed && ( DistNodes.Get( i, j ) > dMax ) ) {
				dMax = DistNodes.Get( i, j );
				pLeft = vecsNodes[ i ].m_pHier;
				pRight = vecsNodes[ j ].m_pHier; } }

	if( dMax >= 0 )
		return;

	dMax = 0;
	for( i = 0; i < Dist.GetSize( ); ++i ) {
		if( vecsGenes[ i ].m_fUsed )
			continue;
		for( j = ( i + 1 ); j < Dist.GetSize( ); ++j )
			if( !vecsGenes[ j ].m_fUsed ) {
				pLeft = vecsGenes[ i ].m_pHier;
				pRight = vecsGenes[ j ].m_pHier;
				return; } }
	for( i = 0; i < iNodes; ++i ) {
		if( vecsNodes[ i ].m_fUsed )
			continue;
		for( j = 0; j < vecvecdDist[ i ].size( ); ++j )
			if( !vecsGenes[ j ].m_fUsed ) {
				pLeft = vecsNodes[ i ].m_pHier;
				pRight = vecsGenes[ j ].m_pHier;
				return; }
		for( j = ( i + 1 ); j < iNodes; ++j )
			if( !vecsNodes[ j ].m_fUsed ) {
				pLeft = vecsNodes[ i ].m_pHier;
				pRight = vecsNodes[ j ].m_pHier;
				return; } }

	g_CatBioUtils.error( "CClustHierarchicalImpl::FindHighestSimilarity( ) found no unused nodes!" ); }

void CClustHierarchicalImpl::CleanDistance( const CHierarchy& Hier,
	vector<SHierarchy>& vecsGenes, vector<SHierarchy>& vecsNodes, TVecVecD& vecvecdDist ) {

	if( Hier.IsGene( ) ) {
		vecsGenes[ Hier.GetID( ) ].m_fUsed = true;
		return; }

	vecsNodes[ Hier.GetID( ) ].m_fUsed = true;
	if( !Hier.Get( false ).IsGene( ) )
		vecvecdDist[ Hier.Get( false ).GetID( ) ].clear( );
	if( !Hier.Get( true ).IsGene( ) )
		vecvecdDist[ Hier.Get( true ).GetID( ) ].clear( ); }

void CClustHierarchicalImpl::RefreshDistances( const CDistanceMatrix& Dist,
	const vector<SHierarchy>& vecsGenes, const vector<SHierarchy>& vecsNodes,
	TVecVecD& vecvecdDist, CDistanceMatrix& DistNodes, size_t iNode ) {
	size_t				i;
	vector<float>&		vecdDist	= vecvecdDist[ iNode ];
	const CHierarchy&	Hier		= *vecsNodes[ iNode ].m_pHier;

	vecdDist.resize( Dist.GetSize( ) );
	for( i = 0; i < vecdDist.size( ); ++i )
		if( !vecsGenes[ i ].m_fUsed )
			vecdDist[ i ] = CalculateDistance( Hier, *vecsGenes[ i ].m_pHier, Dist, vecvecdDist,
				DistNodes );
	for( i = 0; i < iNode; ++i )
		if( !vecsNodes[ i ].m_fUsed )
			DistNodes.Set( i, iNode, CalculateDistance( Hier, *vecsNodes[ i ].m_pHier, Dist,
				vecvecdDist, DistNodes ) ); }

float CClustHierarchicalImpl::CalculateDistance( const CHierarchy& HierNew,
	const CHierarchy& HierOld, const CDistanceMatrix& Dist, const TVecVecD& vecvecdDist,
	const CDistanceMatrix& DistNodes ) {
	size_t	iLeftWeight, iRightWeight;
	float	dLeftDist, dRightDist;

	if( HierNew.IsGene( ) ) {
		if( HierOld.IsGene( ) ) {
			iLeftWeight = 1;
			dLeftDist = Dist.Get( HierNew.GetID( ), HierOld.GetID( ) );
			iRightWeight = 0; }
		else {
			CalculateDistance( HierOld.Get( false ), HierNew, Dist, vecvecdDist, DistNodes,
				iLeftWeight, dLeftDist );
			CalculateDistance( HierOld.Get( true ), HierNew, Dist, vecvecdDist, DistNodes,
				iRightWeight, dRightDist ); } }
	else {
		CalculateDistance( HierNew.Get( false ), HierOld, Dist, vecvecdDist, DistNodes,
			iLeftWeight, dLeftDist );
		CalculateDistance( HierNew.Get( true ), HierOld, Dist, vecvecdDist, DistNodes,
			iRightWeight, dRightDist ); }

	return ( ( ( iLeftWeight * dLeftDist ) + ( iRightWeight * dRightDist ) ) /
		( iLeftWeight + iRightWeight ) ); }

void CClustHierarchicalImpl::CalculateDistance( const CHierarchy& HierOne,
	const CHierarchy& HierTwo, const CDistanceMatrix& Dist, const TVecVecD& vecvecdDist,
	const CDistanceMatrix& DistNodes, size_t& iWeight, float& dDist ) {

	if( HierOne.IsGene( ) ) {
		iWeight = 1;
		if( HierTwo.IsGene( ) )
			dDist = Dist.Get( HierOne.GetID( ), HierTwo.GetID( ) );
		else
			dDist = vecvecdDist[ HierTwo.GetID( ) ][ HierOne.GetID( ) ]; }
	else {
		iWeight = HierOne.GetWeight( );
		if( HierTwo.IsGene( ) )
			dDist = vecvecdDist[ HierOne.GetID( ) ][ HierTwo.GetID( ) ];
		else
			dDist = DistNodes.Get( HierOne.GetID( ), HierTwo.GetID( ) ); } }

}
