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
#include "clusthierarchical.h"
#include "dat.h"
#include "measure.h"

namespace Sleipnir {

/*!
 * \brief
 * Constructs a new hierarchy node.
 * 
 * \param iID
 * Unique ID of the new node.
 * 
 * \param dSimilarity
 * Height of the new node within the tree.
 * 
 * \param pLeft
 * First child of the new node, possibly null.
 * 
 * \param pRight
 * Second child of the new node, possibly null.
 * 
 * \remarks
 * If either of pLeft or pRight is null, they should both be null (i.e. a node should have exactly zero or
 * exactly two children).
 * 
 * \see
 * CClustHierarchical::Cluster
 */
CHierarchy::CHierarchy( size_t iID, float dSimilarity, const CHierarchy* pLeft,
	const CHierarchy* pRight ) {

	assert( ( pLeft && pRight ) || !( pLeft || pRight ) );

	m_iID = iID;
	m_dScore = dSimilarity;
	m_pLeft = pLeft;
	m_pRight = pRight;
	m_iWeight = ( m_pLeft && m_pRight ) ? ( m_pLeft->m_iWeight + m_pRight->m_iWeight ) : 1; }

CHierarchyImpl::~CHierarchyImpl( ) {

	if( m_pLeft )
		delete m_pLeft;
	if( m_pRight )
		delete m_pRight; }

bool CHierarchyImpl::Save( std::ostream& ostm, size_t iNode, const vector<string>* pvecstrGenes ) const {

	if( IsGene( ) )
		return false;

	if( iNode == m_iID ) {
		ostm << GetSave( ) << '\t' << m_pLeft->GetSave( pvecstrGenes ) << '\t' <<
			m_pRight->GetSave( pvecstrGenes ) << '\t' << m_dScore << endl;
		return true; }

	return ( ((const CHierarchyImpl*)m_pLeft)->Save( ostm, iNode, pvecstrGenes ) ||
		((const CHierarchyImpl*)m_pRight)->Save( ostm, iNode, pvecstrGenes ) ); }

string CHierarchyImpl::GetSave( const vector<string>* pvecstrGenes ) const {
	string	strRet;
	char	achBuf[ 16 ];

	if( !( pvecstrGenes && IsGene( ) ) ) {
		strRet = IsGene( ) ? "GENE" : "NODE";
		sprintf_s( achBuf, "%d", m_iID );
		strRet += achBuf; }
	else
		strRet = (*pvecstrGenes)[ m_iID ];

	return strRet; }

/*!
 * \brief
 * Retrieve the IDs of the hierarchy's genes by inorder traversal.
 * 
 * \param veciGenes
 * Output vector into which gene IDs are placed.
 * 
 * \remarks
 * Genes are leaf nodes with IDs generally corresponding to their original PCL indices before clustering.
 * This method will return the PCL indices as they are currently ordered in the hierarchy.
 * 
 * \see
 * SortChildren
 */
void CHierarchy::GetGenes( vector<size_t>& veciGenes ) const {

	if( IsGene( ) )
		veciGenes.push_back( m_iID );
	else {
		m_pLeft->GetGenes( veciGenes );
		m_pRight->GetGenes( veciGenes ); } }

/*!
 * \brief
 * Performs node-flipping in the hierarchy according to a given set of leaf node scores.
 * 
 * \param vecdScores
 * Scores for the hierarchy's leaf nodes indexed by ID.
 * 
 * \returns
 * Score of the current node.
 * 
 * SortChildren can be used to node-flip a hierarchy such that an inorder traversal of the leaf nodes results
 * in a strictly increasing value for some precomputed score.  Since optimal node ordering is NP-hard, this
 * is often used to heuristically order microarray vectors, e.g. from most green to most red.
 * 
 * \remarks
 * vecdScores must be of size equal to the total number of leaves in the hierarchy, and elements of the
 * vector are indexed by the IDs of the leaf nodes (generally the original index of the leaf genes within
 * the pre-clustered PCL file).
 */
float CHierarchy::SortChildren( const vector<float>& vecdScores ) {
	float				dLeft, dRight, dRet;
	const CHierarchy*	pTemp;

	if( IsGene( ) )
		return vecdScores[ m_iID ];

	dLeft = ((CHierarchy*)m_pLeft)->SortChildren( vecdScores );
	dRight = ((CHierarchy*)m_pRight)->SortChildren( vecdScores );
	dRet = ( ( dLeft * m_pLeft->m_iWeight ) + ( dRight * m_pRight->m_iWeight ) ) /
		( m_pRight->m_iWeight + m_pLeft->m_iWeight );
	if( dLeft < dRight ) {
		pTemp = m_pLeft;
		m_pLeft = m_pRight;
		m_pRight = pTemp; }

	return dRet; }

/*!
 * \brief
 * Hierarchically cluster a subset of elements with the given pairwise similarity scores.
 * 
 * \param MatSimilarities
 * Matrix of similarity scores between each pair of elements.
 * 
 * \param vecfIncluded
 * Vector indicating which elements should be clustered (true to be included, false to be excluded).
 * 
 * \returns
 * Resulting clustered hierarchy; null on failure.
 * 
 * Efficiently performs hierarchical clustering on a subset of elements using a precalculated pairwise
 * similarity matrix.  This similarity matrix is often generated by applying some similarity measure to the
 * genes in a PCL file.  The leaves of the resulting hierarchy represent elements of the input matrix and are
 * identified by their indices in this matrix.  Clusters are joined using UPGMA (i.e. average linkage); the
 * pair with the greatest similarity is joined first, then the next greatest, and so forth.  Indices for which
 * vecfIncluded is false are ignored.
 * 
 * \remarks
 * vecfIncluded must be of the same size as MatSimilarities.
 * 
 * \see
 * CClustKMeans::Cluster
 */
CHierarchy* CClustHierarchical::Cluster( const CDistanceMatrix& MatSimilarities,
	const vector<bool>& vecfIncluded ) {

	return CClustHierarchicalImpl::Cluster( MatSimilarities, &vecfIncluded ); }

/*!
 * \brief
 * Hierarchically cluster a set of elements with the given pairwise similarity scores.
 * 
 * \param MatSimilarities
 * Matrix of similarity scores between each pair of elements.
 * 
 * \returns
 * Resulting clustered hierarchy; null on failure.
 * 
 * Efficiently performs hierarchical clustering on a set of elements using a precalculated pairwise similarity
 * matrix.  This similarity matrix is often generated by applying some similarity measure to the genes in a
 * PCL file.  The leaves of the resulting hierarchy represent elements of the input matrix and are identified
 * by their indices in this matrix.  Clusters are joined using UPGMA (i.e. average linkage); the pair with
 * the greatest similarity is joined first, then the next greatest, and so forth.
 * 
 * \see
 * CClustKMeans::Cluster
 */
CHierarchy* CClustHierarchical::Cluster( const CDistanceMatrix& MatSimilarities ) {

	return CClustHierarchicalImpl::Cluster( MatSimilarities ); }

// Implementation inspired by TIGR MeV

CHierarchy* CClustHierarchicalImpl::Cluster( const CDistanceMatrix& Dist, const vector<bool>* pvecfGenes ) {
	CDistanceMatrix	Sim;
	size_t			i, j, k, iAssigned, iParentless, iOne, iTwo;
	float			d, dTotal, dMin;
	vector<float>	vecdHeight, vecdMax;
	vector<size_t>	veciChild1, veciChild2, veciChildren, veciMax, veciOwner;

	if( !Dist.GetSize( ) )
		return NULL;

	dMin = FLT_MAX;
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
	for( i = 0; i < Sim.GetSize( ); ++i )
		for( j = ( i + 1 ); j < Sim.GetSize( ); ++j )
			if( ( ( d = Sim.Get( i, j ) ) == -FLT_MAX ) || CMeta::IsNaN( d ) ) {
				g_CatSleipnir.error( "CClustHierarchicalImpl::Cluster( ) illegal input value at %d, %d", i,
					j );
				return NULL; }
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
			g_CatSleipnir.notice( "CClustHierarchical::Cluster( ) %d/%d nodes remaining", iParentless,
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
