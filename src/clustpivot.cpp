#include "stdafx.h"
#include "clustpivot.h"
#include "meta.h"

namespace Sleipnir {

uint16_t CClustPivot::Cluster( const CDistanceMatrix& Dist, float dCutoff,
	vector<uint16_t>& vecsClusters ) {
	size_t			i, j, iRand, iTmp, iPivot;
	uint16_t		sRet;
	vector<size_t>	veciPerm;
	float			d;

	if( vecsClusters.size( ) != Dist.GetSize( ) )
		return -1;

	veciPerm.resize( Dist.GetSize( ) );
	// Pick a random permutation of the genes
	for( i = 0; i < veciPerm.size( ); ++i )
		veciPerm[ i ] = i;
	for( i = 0; i < Dist.GetSize( ); ++i ) {
		iRand = rand( ) % ( veciPerm.size( ) - i );
		iTmp = veciPerm[ i ];
		veciPerm[ i ] = veciPerm[ i + iRand ];
		veciPerm[ i + iRand ] = iTmp; }

	// reset the cluster data
	for( i = 0; i < vecsClusters.size( ); ++i )
		vecsClusters[ i ] = -1;

	for( sRet = i = 0; i < Dist.GetSize( ); ++i ) {
		iPivot = veciPerm[ i ];
		// If gene was already clustered (or excluded), continue
		if( vecsClusters[ iPivot ] != -1 )
			continue;

		vecsClusters[ iPivot ] = sRet++;
		for( j = 0; j < Dist.GetSize( ); ++j ) {
			// check if already clustered (or thrown away)
			if( vecsClusters[ j ] != -1 )
				continue;

		if( !CMeta::IsNaN( d = Dist.Get( iPivot, j ) ) && ( d > dCutoff ) )
			vecsClusters[ j ] = vecsClusters[ iPivot ]; } }

	return sRet; }

}
