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
#include "clustpivot.h"
#include "meta.h"

namespace Sleipnir {

/*!
 * \brief
 * Implements a randomized pivot clustering algorithm due to Ailon and Charikar.
 * 
 * \param MatSimilarities
 * Matrix of similarity scores between each pair of elements.
 * 
 * \param dCutoff
 * Description of parameter dCutoff.
 * 
 * \param vecsClusters
 * Output cluster IDs for each gene.
 * 
 * \returns
 * Total number of clusters.
 * 
 * Performs a hard clustering of a set of elements using the randomized pivot algorithm due to Ailon and
 * Charikar.  This places each gene in exactly one cluster.  Briefly, an unclustered gene is chosen at
 * random to be the current pivot.  Each gene more similar than the given cutoff is placed in that pivot's
 * cluster.  This process is iterated until all genes are clustered.  This algorithm has some nice
 * theoretical properties, but in practice, it doesn't do so hot at uncovering useful biological information.
 * 
 * \remarks
 * On successful return, the size of vecsClusters will be equal to the size of MatSimilarities.
 * 
 * \see
 * CClustKMeans::Cluster
 */
uint16_t CClustPivot::Cluster( const CDistanceMatrix& MatSimilarities, float dCutoff,
	vector<uint16_t>& vecsClusters ) {
	size_t			i, j, iRand, iTmp, iPivot;
	uint16_t		sRet;
	vector<size_t>	veciPerm;
	float			d;

	vecsClusters.resize( MatSimilarities.GetSize( ) );
	veciPerm.resize( MatSimilarities.GetSize( ) );
	// Pick a random permutation of the genes
	for( i = 0; i < veciPerm.size( ); ++i )
		veciPerm[ i ] = i;
	for( i = 0; i < MatSimilarities.GetSize( ); ++i ) {
		iRand = rand( ) % ( veciPerm.size( ) - i );
		iTmp = veciPerm[ i ];
		veciPerm[ i ] = veciPerm[ i + iRand ];
		veciPerm[ i + iRand ] = iTmp; }

	// reset the cluster data
	for( i = 0; i < vecsClusters.size( ); ++i )
		vecsClusters[ i ] = -1;

	for( sRet = i = 0; i < MatSimilarities.GetSize( ); ++i ) {
		iPivot = veciPerm[ i ];
		// If gene was already clustered (or excluded), continue
		if( vecsClusters[ iPivot ] != (uint16_t)-1 )
			continue;

		vecsClusters[ iPivot ] = sRet++;
		for( j = 0; j < MatSimilarities.GetSize( ); ++j ) {
			// check if already clustered (or thrown away)
			if( vecsClusters[ j ] != (uint16_t)-1 )
				continue;

		if( !CMeta::IsNaN( d = MatSimilarities.Get( iPivot, j ) ) && ( d > dCutoff ) )
			vecsClusters[ j ] = vecsClusters[ iPivot ]; } }

	return sRet; }

}
