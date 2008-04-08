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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
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
			g_CatSleipnir.info( "CClustKMeans::Cluster( %d ) iteration %d", iK, iIteration );
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
