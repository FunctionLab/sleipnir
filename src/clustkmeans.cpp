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
#include "clustkmeans.h"
#include "measure.h"
#include "meta.h"

namespace Sleipnir {

/*!
 * \brief
 * Cluster a set of elements into k groups using the given data and pairwise similarity score.
 * 
 * \param MatData
 * Data vectors for each element, generally microarray values from a PCL file.
 * 
 * \param pMeasure
 * Similarity measure to use for clustering.
 * 
 * \param iK
 * Number of clusters to generate.
 * 
 * \param vecsClusters
 * Output cluster IDs for each gene.
 * 
 * \param pMatWeights
 * If non-null, weights to use for each gene/condition value.  These can be used to up/downweight aneuploidies
 * present under only certain conditions, for example.  Default assumes all ones.
 * 
 * \returns
 * True if clustering succeeded.
 * 
 * Performs k-means clustering on the given data using the specified similarity measure and number of
 * clusters.  The indices of each element's final cluster are indicated in the output vector.  If given,
 * individual gene/condition scores can be weighted (e.g. to up/downweight aneuploidies present only under
 * certain conditions).  During k-means clustering, K centers are initially chosen at random.  Each gene is
 * assigned to the center most similar to it, and the centers are moved to the mean of their assigned
 * genes.  This process is iterated until no gene assignments change.  This places each gene in exactly one
 * cluster.
 * 
 * \remarks
 * The size of MatData must be at least iK; on successful return, the size of vecsClusters will be equal to
 * the size of MatData.
 * 
 * \see
 * CClustHierarchical::Cluster
 */
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
			g_CatSleipnir.notice( "CClustKMeans::Cluster( %d ) iteration %d", iK, iIteration );
		fill( veciCounts.begin( ), veciCounts.end( ), 0 );
		fDone = true;
		for( i = 0; i < vecsClusters.size( ); ++i ) {
			float	dMax;

			dMax = -FLT_MAX;
			for( sMax = j = 0; j < MatMeans.GetRows( ); ++j ) {
				d = (float)pMeasure->Measure( MatData.Get( i ), MatData.GetColumns( ), MatMeans.Get( j ),
					MatMeans.GetColumns( ), IMeasure::EMapNone, pMatWeights ? pMatWeights->Get( i ) : NULL );
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

bool CClustKMeans::Cluster( const CDistanceMatrix& MatSimilarities, size_t iK,
	vector<uint16_t>& vecsClusters ) {
	size_t			i, j, iOne, iIteration, iChanged, iState;
	float			d, dMax, dMin;
	CDataMatrix		MatPrev, MatNext;
	vector<size_t>	veciPrev, veciNext;
	uint16_t		sMax;
	set<size_t>		setiStates;

	if( MatSimilarities.GetSize( ) < iK )
		return false;

	dMax = -( dMin = FLT_MAX );
	for( i = 0; i < MatSimilarities.GetSize( ); ++i )
		for( j = ( i + 1 ); j < MatSimilarities.GetSize( ); ++j )
			if( !CMeta::IsNaN( d = MatSimilarities.Get( i, j ) ) ) {
				if( d > dMax )
					dMax = d;
				if( d < dMin )
					dMin = d; }
	if( dMin == FLT_MAX )
		return false;
	dMax++;
	dMin--;
	MatPrev.Initialize( MatSimilarities.GetSize( ), iK );
	for( i = 0; i < MatPrev.GetColumns( ); ++i ) {
		iOne = rand( ) % MatSimilarities.GetSize( );
		for( j = 0; j < MatPrev.GetRows( ); ++j )
			MatPrev.Set( j, i, GetClean( iOne, j, dMin, dMax, MatSimilarities ) ); }
	MatNext.Initialize( MatPrev.GetRows( ), MatPrev.GetColumns( ) );
	MatNext.Clear( );

	vecsClusters.resize( MatSimilarities.GetSize( ) );
	fill( vecsClusters.begin( ), vecsClusters.end( ), iK );
	veciPrev.resize( iK );
	fill( veciPrev.begin( ), veciPrev.end( ), 1 );
	veciNext.resize( veciPrev.size( ) );
	for( iChanged = MatSimilarities.GetSize( ),iIteration = 0; iChanged > 2; ++iIteration ) {
		g_CatSleipnir.log( ( iIteration % 10 ) ? Priority::DEBUG : Priority::NOTICE,
			"CClustKMeans::Cluster( %d ) iteration %d", iK, iIteration );
		for( iChanged = i = 0; i < vecsClusters.size( ); ++i ) {
			float	dMax;

			dMax = -FLT_MAX;
			for( sMax = j = 0; j < MatPrev.GetColumns( ); ++j ) {
				if( CMeta::IsNaN( d = MatPrev.Get( i, j ) ) ) {
					g_CatSleipnir.error( "CClustKMeans::Cluster( %d ) failed on gene %d, cluster %d", i, j );
					return false; }
				d /= veciPrev[ j ];
				if( d > dMax ) {
					dMax = d;
					sMax = j; } }
			if( vecsClusters[ i ] != sMax ) {
				iChanged++;
				if( vecsClusters[ i ] != iK )
					veciNext[ vecsClusters[ i ] ]--;
				veciNext[ sMax ]++;
				for( j = 0; j < MatSimilarities.GetSize( ); ++j ) {
					d = GetClean( i, j, dMin, dMax, MatSimilarities );
					if( vecsClusters[ i ] != iK )
						MatNext.Get( j, vecsClusters[ i ] ) -= d;
					MatNext.Get( j, sMax ) += d; }
				vecsClusters[ i ] = sMax; } }

		for( i = 0; i < veciNext.size( ); ++i )
			if( !veciNext[ i ] ) {
				do
					iOne = rand( ) % vecsClusters.size( );
				while( veciNext[ vecsClusters[ iOne ] ] < 2 );
				g_CatSleipnir.info( "CClustKMeans::Cluster( %d ) moving gene %d into empty cluster %d", iK,
					iOne, i );
				veciNext[ vecsClusters[ iOne ] ]--;
				for( j = 0; j < MatNext.GetRows( ); ++j ) {
					d = GetClean( iOne, j, dMin, dMax, MatSimilarities );
					MatNext.Get( j, vecsClusters[ iOne ] ) -= d;
					MatNext.Set( j, i, d ); }
				veciNext[ i ]++; }

// This calculates a simple hash for the current cluster assignments
		for( iState = i = 0; i < vecsClusters.size( ); ++i ) {
			j = iState & (size_t)-32768;
			iState <<= 15;
			iState ^= j >> ( ( sizeof(size_t) * 8 ) - 15 );
			iState ^= vecsClusters[ i ] + 1; }
		if( setiStates.find( iState ) != setiStates.end( ) ) {
			g_CatSleipnir.info( "CClustKMeans::Cluster( %d ) found redundant state, terminating", iK );
			break; }
		setiStates.insert( iState );

		g_CatSleipnir.notice( "CClustKMeans::Cluster( %d ) updated %d genes", iK, iChanged );
		if( g_CatSleipnir.isDebugEnabled( ) )
			for( i = 0; i < MatPrev.GetRows( ); ++i ) {
				ostringstream	ossm;

				for( j = 0; j < MatPrev.GetColumns( ); ++j )
					ossm << ( j ? "\t" : "" ) << MatPrev.Get( i, j );
				g_CatSleipnir.debug( "CClustKMeans::Cluster( %d ) object %d:	%s", iK, i,
					ossm.str( ).c_str( ) ); }
		copy( veciNext.begin( ), veciNext.end( ), veciPrev.begin( ) );
		for( i = 0; i < MatPrev.GetRows( ); ++i )
			memcpy( MatPrev.Get( i ), MatNext.Get( i ), MatPrev.GetColumns( ) * sizeof(*MatPrev.Get( i )) ); }

	return true; }

void CClustKMeansImpl::Randomize( CDataMatrix& MatMeans, size_t iRow, const CDataMatrix& MatData ) {

	MatMeans.Set( iRow, MatData.Get( rand( ) % MatData.GetRows( ) ) ); }

}
