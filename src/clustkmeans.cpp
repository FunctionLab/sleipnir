#include "stdafx.h"
#include "clustkmeans.h"
#include "measure.h"
#include "meta.h"

namespace libBioUtils {

bool CClustKMeans::Cluster( const CDataMatrix& MatData, const IMeasure* pMeasure, size_t iK,
	vector<uint16_t>& vecsClusters, const CDataMatrix* pMatWeights ) {
	size_t			i, j, iIteration;
	float			d;
	CDataMatrix		MatMeans;
	bool			fDone;
	vector<size_t>	veciCounts;
	uint16_t		sMax;

	if( !pMeasure || ( MatData.GetRows( ) < iK ) )
		return false;

	MatMeans.Initialize( iK, MatData.GetColumns( ) );
	for( i = 0; i < MatMeans.GetRows( ); ++i )
		Randomize( MatMeans, i, MatData );

	vecsClusters.resize( MatData.GetRows( ) );
	fill( vecsClusters.begin( ), vecsClusters.end( ), -1 );
	veciCounts.resize( iK );
	for( iIteration = 0,fDone = false; !fDone; ++iIteration ) {
		if( !( iIteration % 10 ) )
			g_CatBioUtils.info( "CClustKMeans::Cluster( %d ) iteration %d", iK, iIteration );
		fill( veciCounts.begin( ), veciCounts.end( ), 0 );
		fDone = true;
		for( i = 0; i < vecsClusters.size( ); ++i ) {
			float	dMax;

			dMax = -FLT_MAX;
			for( sMax = j = 0; j < MatMeans.GetRows( ); ++j ) {
				d = (float)pMeasure->Measure( MatData.Get( i ), MatData.GetColumns( ), MatMeans.Get( j ),
					MatMeans.GetColumns( ), IMeasure::EMapCenter, pMatWeights ? pMatWeights->Get( i ) : NULL );
				if( CMeta::IsNaN( d ) )
					return false;
				if( d > dMax ) {
					dMax = d;
					sMax = j; } }
			veciCounts[ sMax ]++;
			if( vecsClusters[ i ] != sMax ) {
				fDone = false;
				vecsClusters[ i ] = sMax; } }

		MatMeans.Clear( );
		for( i = 0; i < vecsClusters.size( ); ++i )
			for( j = 0; j < MatMeans.GetColumns( ); ++j )
				MatMeans.Get( vecsClusters[ i ], j ) += MatData.Get( i, j );
		for( i = 0; i < MatMeans.GetRows( ); ++i )
			if( veciCounts[ i ] )
				for( j = 0; j < MatMeans.GetColumns( ); ++j )
					MatMeans.Get( i, j ) /= veciCounts[ i ];
			else
				Randomize( MatMeans, i, MatData ); }

	return true; }

void CClustKMeansImpl::Randomize( CDataMatrix& MatMeans, size_t iRow, const CDataMatrix& MatData ) {

	MatMeans.Set( iRow, MatData.Get( rand( ) % MatData.GetRows( ) ) ); }

}
