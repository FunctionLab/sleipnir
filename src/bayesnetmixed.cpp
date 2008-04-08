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
#ifdef PNL_ENABLED
#pragma warning (disable: 4244 4267)
#include <pnl_dll.hpp>
#pragma warning (default: 4244 4267)
#include "bayesnet.h"

namespace Sleipnir {

bool CBayesNetSmileImpl::ConvertGraph( CBayesNetPNL& BNPNL ) const {
	CNodeType*						aNodeTypes;
	DSL_intArray					IntArrayNodes, IntArrayAdj;
	int								i, j, iSize;
	int**							aaiNeighborCounts;
	int*							aiNeighborCounts;
	int*							aiNodeAssociation;
	map<int,int>					mapSmilePNL, mapNodeTypes;
	ENeighborType**					aaeNeighborTypes;
	vector<int>						veciNeighborCounts;
	vector<ENeighborType>			veceNeighborTypes;
	map<int,int>::const_iterator	iterNodeTypes;
	DSL_node*						pNode;
	CGraph*							pGraph;

	m_SmileNet.GetAllNodes( IntArrayNodes );
	for( i = 0; i < IntArrayNodes.NumItems( ); ++i )
		mapSmilePNL[ IntArrayNodes[ i ] ] = i;

	aiNodeAssociation = new int[ IntArrayNodes.NumItems( ) ];
	for( i = 0; i < IntArrayNodes.NumItems( ); ++i ) {
		pNode = m_SmileNet.GetNode( IntArrayNodes[ i ] );
		if( IsGaussian( m_SmileNet ) )
			iSize = -1;
		else
			iSize = pNode->Definition( )->GetNumberOfOutcomes( );
		if( ( iterNodeTypes = mapNodeTypes.find( iSize ) ) == mapNodeTypes.end( ) ) {
			j = (int)mapNodeTypes.size( );
			mapNodeTypes[ iSize ] = j;
			iterNodeTypes = mapNodeTypes.find( iSize ); }
		aiNodeAssociation[ i ] = iterNodeTypes->second; }
	aNodeTypes = new CNodeType[ mapNodeTypes.size( ) ];
	for( iterNodeTypes = mapNodeTypes.begin( ); iterNodeTypes != mapNodeTypes.end( );
		++iterNodeTypes )
		aNodeTypes[ iterNodeTypes->second ].SetType( ( iterNodeTypes->first > 0 ),
			abs( iterNodeTypes->first ) );

	aiNeighborCounts = new int[ IntArrayNodes.NumItems( ) ];
	aaiNeighborCounts = new int*[ IntArrayNodes.NumItems( ) ];
	aaeNeighborTypes = new ENeighborType*[ IntArrayNodes.NumItems( ) ];
	for( i = 0; i < IntArrayNodes.NumItems( ); ++i ) {
		IntArrayAdj = m_SmileNet.GetParents( IntArrayNodes[ i ] );
		veciNeighborCounts.resize( iSize = aiNeighborCounts[ i ] = IntArrayAdj.NumItems( ) );
		veceNeighborTypes.resize( iSize );
		for( j = 0; j < IntArrayAdj.NumItems( ); ++j ) {
			veciNeighborCounts[ j ] = mapSmilePNL[ IntArrayAdj[ j ] ];
			veceNeighborTypes[ j ] = ntParent; }

		IntArrayAdj = m_SmileNet.GetChildren( IntArrayNodes[ i ] );
		veciNeighborCounts.resize( aiNeighborCounts[ i ] += IntArrayAdj.NumItems( ) );
		veceNeighborTypes.resize( aiNeighborCounts[ i ] );
		for( j = 0; j < IntArrayAdj.NumItems( ); ++j ) {
			veciNeighborCounts[ iSize + j ] = mapSmilePNL[ IntArrayAdj[ j ] ];
			veceNeighborTypes[ iSize + j ] = ntChild; }

		aaiNeighborCounts[ i ] = new int[ veciNeighborCounts.size( ) ];
		aaeNeighborTypes[ i ] = new ENeighborType[ veceNeighborTypes.size( ) ];
		for( j = 0; j < veciNeighborCounts.size( ); ++j ) {
			aaiNeighborCounts[ i ][ j ] = veciNeighborCounts[ j ];
			aaeNeighborTypes[ i ][ j ] = veceNeighborTypes[ j ]; } }

	if( BNPNL.m_pPNLNet )
		delete BNPNL.m_pPNLNet;
	pGraph = CGraph::Create( IntArrayNodes.NumItems( ), aiNeighborCounts,
		aaiNeighborCounts, aaeNeighborTypes );
	BNPNL.m_pPNLNet = CBNet::Create( IntArrayNodes.NumItems( ), (int)mapNodeTypes.size( ),
		aNodeTypes, aiNodeAssociation, pGraph );

	for( i = 0; i < IntArrayNodes.NumItems( ); ++i ) {
		delete[] aaiNeighborCounts[ i ];
		delete[] aaeNeighborTypes[ i ]; }
	delete[] aaiNeighborCounts;
	delete[] aaeNeighborTypes;
	delete[] aiNeighborCounts;
	delete[] aiNodeAssociation;
	delete[] aNodeTypes;

	return true; }

bool CBayesNetSmileImpl::ConvertCPTs( CBayesNetPNL& BNPNL ) const {
	int							i, j, k, iCPT;
	floatVector					vecdCPT;
	DSL_Dmatrix*				pCPT;
	DSL_intArray				IntArrayNodes, IntArrayCoords;
	DSL_node*					pNode;
	CNumericDenseMatrix<float>*	pMatrix;
	int							aiDims[ 2 ]	= { 1, 1 };
	float						dCPT;

	BNPNL.m_pPNLNet->AllocFactors( );
	m_SmileNet.GetAllNodes( IntArrayNodes );
	for( i = 0; i < IntArrayNodes.NumItems( ); ++i ) {
		pNode = m_SmileNet.GetNode( IntArrayNodes[ i ] );
		pCPT = pNode->Definition( )->GetMatrix( );
		if( IsGaussian( m_SmileNet ) ) {
			BNPNL.m_pPNLNet->AllocFactor( i );

			dCPT = (float)(*pCPT)[ iCPT = 0 ];
			pMatrix = CNumericDenseMatrix<float>::Create( 2, aiDims, &dCPT );
			BNPNL.m_pPNLNet->GetFactor( i )->AttachMatrix( pMatrix, matMean );

			dCPT = (float)(*pCPT)[ ++iCPT ];
			pMatrix = CNumericDenseMatrix<float>::Create( 2, aiDims, &dCPT );
			BNPNL.m_pPNLNet->GetFactor( i )->AttachMatrix( pMatrix, matCovariance );

			j = m_SmileNet.GetParents( IntArrayNodes[ i ] ).NumItems( );
			for( k = 0; k < j; ++k ) {
				dCPT = (float)(*pCPT)[ ++iCPT ];
				pMatrix = CNumericDenseMatrix<float>::Create( 2, aiDims, &dCPT );
				BNPNL.m_pPNLNet->GetFactor( i )->AttachMatrix( pMatrix, matWeights, k ); } }
		else {
			vecdCPT.resize( pCPT->GetSize( ) );
			for( j = 0; j < pCPT->GetSize( ); ++j )
				vecdCPT[ j ] = (float)(*pCPT)[ j ];
			BNPNL.m_pPNLNet->CreateTabularCPD( i, vecdCPT ); } }

	return true; }

}

#endif // PNL_ENABLED
