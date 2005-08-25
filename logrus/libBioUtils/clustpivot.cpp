#include "stdafx.h"
#include "clustpivot.h"
#include "meta.h"

namespace libBioUtils {

size_t CClustPivot::Cluster( const CDistanceMatrix& Dist, float dCutoff,
	vector<size_t>& veciClusters ) {
	size_t			i, j, iRand, iTmp, iRet, iPivot;
	vector<size_t>	veciPerm;
	float			d;

	if( veciClusters.size( ) != Dist.GetSize( ) )
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
	for( i = 0; i < veciClusters.size( ); ++i )
		veciClusters[ i ] = -1;

	for( iRet = i = 0; i < Dist.GetSize( ); ++i ) {
		iPivot = veciPerm[ i ];
		// If gene was already clustered (or excluded), continue
		if( veciClusters[ iPivot ] != -1 )
			continue;

		veciClusters[ iPivot ] = iRet++;
		for( j = 0; j < Dist.GetSize( ); ++j ) {
			// check if already clustered (or thrown away)
			if( veciClusters[ j ] != -1 )
				continue;

		if( !CMeta::IsNaN( d = Dist.Get( iPivot, j ) ) && ( d > dCutoff ) )
			veciClusters[ j ] = veciClusters[ iPivot ]; } }

	return iRet; }

}
