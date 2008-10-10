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
#include "clustqtc.h"
#include "pcl.h"
#include "measure.h"
#include "meta.h"

namespace Sleipnir {

/*!
 * \brief
 * Cluster a set of elements with the quality threshold algorithm using the given data and pairwise
 * similarity score.
 * 
 * \param MatData
 * Data vectors for each element, generally microarray values from a PCL file.
 * 
 * \param pMeasure
 * Similarity measure to use for clustering.
 * 
 * \param dDiameter
 * Maximum cluster diameter.
 * 
 * \param iSize
 * Minimum cluster size.
 * 
 * \param vecsClusters
 * Output cluster IDs for each gene; unclustered genes are grouped in the last cluster.
 * 
 * \param pMatWeights
 * If non-null, weights to use for each gene/condition value.  These can be used to up/downweight aneuploidies
 * present under only certain conditions, for example.  Default assumes all ones.
 * 
 * \returns
 * Total number of clusters.
 * 
 * Clusters elements using the quality threshold algorithm due to Heyer et al.  Each gene is assigned to at
 * most one cluster.  Briefly, the most similar gene pair is grouped together, and each other gene within
 * the given diameter of that group is added to the cluster.  These genes are then removed from the pool and
 * the process is repeated.  Gene groups that cannot reach the minimum cluster size are discarded.
 * 
 * \remarks
 * If N clusters are generated, unclustered genes will be assigned ID N in the output vector.
 * 
 * \see
 * CClustKMeans::Cluster
 */
uint16_t CClustQTC::Cluster( const CDataMatrix& MatData, const IMeasure* pMeasure,
	float dDiameter, size_t iSize, vector<uint16_t>& vecsClusters, const CDataMatrix* pMatWeights ) {
	CDistanceMatrix	Dist;

	InitializeDistances( MatData, pMeasure, Dist, pMatWeights );
	return QualityThresholdAll( MatData.GetRows( ), dDiameter, iSize, Dist, vecsClusters ); }

/*!
 * \brief
 * Cluster a set of elements with the quality threshold algorithm using the given similarity scores.
 * 
 * \param MatSimilarities
 * Precalculated similarity scores, generally using microarray values from a PCL file.
 * 
 * \param dDiameter
 * Maximum cluster diameter.
 * 
 * \param iSize
 * Minimum cluster size.
 * 
 * \param vecsClusters
 * Output cluster IDs for each gene; unclustered genes are grouped in the last cluster.
 * 
 * \returns
 * Total number of clusters.
 * 
 * Clusters elements using the quality threshold algorithm due to Heyer et al.  Each gene is assigned to at
 * most one cluster.  Briefly, the most similar gene pair is grouped together, and each other gene within
 * the given diameter of that group is added to the cluster.  These genes are then removed from the pool and
 * the process is repeated.  Gene groups that cannot reach the minimum cluster size are discarded.
 * 
 * \remarks
 * If N clusters are generated, unclustered genes will be assigned ID N in the output vector.
 * 
 * \see
 * CClustKMeans::Cluster
 */
uint16_t CClustQTC::Cluster( const CDistanceMatrix& MatSimilarities, float dDiameter, size_t iSize,
	vector<uint16_t>& vecsClusters ) {

	return QualityThresholdAll( MatSimilarities.GetSize( ), dDiameter, iSize, MatSimilarities,
		vecsClusters ); }

/*!
 * \brief
 * Record the smallest cluster diameter within some range at which each gene pair clusters.
 * 
 * \param MatData
 * Data vectors for each element, generally microarray values from a PCL file.
 * 
 * \param pMeasure
 * Similarity measure to use for clustering.
 * 
 * \param dMinDiameter
 * Minimum cluster diameter at which to attempt clustering.
 * 
 * \param dMaxDiameter
 * Maximum cluster diameter at which to attempt clustering.
 * 
 * \param dDeltaDiameter
 * Increment of cluster diameters to scan between minimum and maximum.
 * 
 * \param iSize
 * Minimum cluster size.
 * 
 * \param MatResults
 * Output matrix recording the smallest diameter at which each gene pair coclustered, or NaN if the pair
 * did not cocluster within the given diameter range.
 * 
 * \param pMatWeights
 * If non-null, weights to use for each gene/condition value.  These can be used to up/downweight aneuploidies
 * present under only certain conditions, for example.  Default assumes all ones.
 * 
 * This clustering method incrementally attempts to quality threshold cluster the given elements at each
 * cluster diameter between the given minimum and maximum, in steps of the requested delta.  For each gene
 * pair, the smallest diameter at which they coclustered (appeared in some cluster together) is recorded.
 * This can be used to rapidly scan through a range of cluster "sizes" to find the strictest diameter
 * cutoff at which gene pairs cocluster.
 * 
 * \remarks
 * MatResults must be pre-initialized to the same size as MatData.
 */
void CClustQTC::Cluster( const CDataMatrix& MatData, const IMeasure* pMeasure,
	float dMinDiameter, float dMaxDiameter, float dDeltaDiameter, size_t iSize, CDistanceMatrix& MatResults,
	const CDataMatrix* pMatWeights ) {
	CDistanceMatrix	Dist;

	InitializeDistances( MatData, pMeasure, Dist, pMatWeights );
	return Cluster( Dist, dMinDiameter, dMaxDiameter, dDeltaDiameter, iSize, MatResults ); }

/*!
 * \brief
 * Record the smallest cluster diameter within some range at which each gene pair clusters.
 * 
 * \param MatSimilarities
 * Precalculated similarity scores for every pair of entities, generally genes from a microarray PCL.
 * 
 * \param dMinDiameter
 * Minimum cluster diameter at which to attempt clustering.
 * 
 * \param dMaxDiameter
 * Maximum cluster diameter at which to attempt clustering.
 * 
 * \param dDeltaDiameter
 * Increment of cluster diameters to scan between minimum and maximum.
 * 
 * \param iSize
 * Minimum cluster size.
 * 
 * \param MatResults
 * Output matrix recording the smallest diameter at which each gene pair coclustered, or NaN if the pair
 * did not cocluster within the given diameter range.
 * 
 * This clustering method incrementally attempts to quality threshold cluster the given elements at each
 * cluster diameter between the given minimum and maximum, in steps of the requested delta.  For each gene
 * pair, the smallest diameter at which they coclustered (appeared in some cluster together) is recorded.
 * This can be used to rapidly scan through a range of cluster "sizes" to find the strictest diameter
 * cutoff at which gene pairs cocluster.
 * 
 * \remarks
 * MatResults must be pre-initialized to the same size as MatSimilarities.
 */
void CClustQTC::Cluster( const CDistanceMatrix& MatSimilarities, float dMinDiameter, float dMaxDiameter,
	float dDeltaDiameter, size_t iSize, CDistanceMatrix& MatResults ) {
	float				dDiameter;
	vector<uint16_t>	vecsClusters;
	uint16_t			sClusters;
	size_t				i, j;

	for( dDiameter = dMinDiameter; dDiameter <= dMaxDiameter; dDiameter += dDeltaDiameter ) {
		g_CatSleipnir.notice( "CClustQTC::Cluster( %g, %g, %g, %d, %d ) processing diameter %g", dMinDiameter,
			dMaxDiameter, dDeltaDiameter, iSize, dDiameter );
		sClusters = QualityThresholdAll( MatSimilarities.GetSize( ), dDiameter, iSize, MatSimilarities,
			vecsClusters );
		for( i = 0; i < vecsClusters.size( ); ++i ) {
			if( ( vecsClusters[ i ] + 1 ) == sClusters )
				continue;
			for( j = ( i + 1 ); j < vecsClusters.size( ); ++j )
				if( ( vecsClusters[ j ] == vecsClusters[ i ] ) && CMeta::IsNaN( MatResults.Get( i, j ) ) )
					MatResults.Set( i, j, 1 - dDiameter ); } } }

uint16_t CClustQTCImpl::QualityThresholdAll( size_t iGenes, float dDiameter, size_t iSize,
	const CDistanceMatrix& Dist, vector<uint16_t>& vecsClusters ) {
	size_t				i, iAssigned;
	uint16_t			sCluster;
	vector<uint16_t>	vecsCur;
	vector<bool>		vecfAssigned;

	vecfAssigned.resize( iGenes );
	for( i = 0; i < vecfAssigned.size( ); ++i )
		vecfAssigned[ i ] = false;

	vecsClusters.resize( iGenes );
	for( iAssigned = sCluster = 0; ; ++sCluster ) {
		g_CatSleipnir.notice( "CClustQTCImpl::QualityThresholdAll( ) cluster %d, assigned %d/%d genes",
			sCluster + 1, iAssigned, iGenes );
		QualityThresholdLargest( iGenes, dDiameter, Dist, vecfAssigned, vecsCur );
		for( i = 0; i < vecsCur.size( ); ++i )
			vecsClusters[ vecsCur[ i ] ] = sCluster;

		if( vecsCur.size( ) < iSize ) {
			for( i = 0; i < vecfAssigned.size( ); ++i )
				if( !vecfAssigned[ i ] )
					vecsClusters[ i ] = sCluster;
			break; }

		if( ( iAssigned += vecsCur.size( ) ) >= iGenes ) {
			sCluster++;
			break; }
		for( i = 0; i < vecsCur.size( ); ++i )
			vecfAssigned[ vecsCur[ i ] ] = true; }

	return ( sCluster + 1 ); }

void CClustQTCImpl::QualityThresholdLargest( size_t iGenes, float dDiameter, const CDistanceMatrix& Dist,
	const vector<bool>& vecfAssigned,
	vector<uint16_t>& vecsCluster ) {
	vector<bool>		vecfClone;
	size_t				iGene;
	vector<uint16_t>	vecsCur;
	vector<float>		vecdDiameter;

	vecsCluster.clear( );
	vecfClone.resize( vecfAssigned.size( ) );
	for( iGene = 0; iGene < vecfAssigned.size( ); ++iGene ) {
		if( !( iGene % 1000 ) )
			g_CatSleipnir.notice( "CClustQTCImpl::QualityThresholdLargest( %g ) processing gene %d/%d",
				dDiameter, iGene, vecfAssigned.size( ) );
		if( vecfAssigned[ iGene ] )
			continue;
		copy( vecfAssigned.begin( ), vecfAssigned.end( ), vecfClone.begin( ) );
		vecfClone[ iGene ] = true;
		QualityThresholdGene( iGene, iGenes, dDiameter, Dist, vecfClone, vecdDiameter, vecsCur );
		if( vecsCur.size( ) > vecsCluster.size( ) ) {
			vecsCluster.resize( vecsCur.size( ) );
			copy( vecsCur.begin( ), vecsCur.end( ), vecsCluster.begin( ) ); } } }

void CClustQTCImpl::QualityThresholdGene( size_t iGene, size_t iGenes, float dDiameter,
	const CDistanceMatrix& Dist, vector<bool>& vecfAssigned, vector<float>& vecdDiameter,
	vector<uint16_t>& vecsCluster ) {
	size_t	iAdded, iLocal, iBest;
	float	dBest;

	vecsCluster.resize( 1 );
	vecsCluster[ 0 ] = iAdded = iGene;
	vecdDiameter.resize( iGenes );
	for( iLocal = 0; iLocal < vecfAssigned.size( ); ++iLocal )
		if( !vecfAssigned[ iLocal ] )
			vecdDiameter[ iLocal ] = -(float)HUGE_VAL;

	while( true ) {
		iBest = -1;
		dBest = (float)HUGE_VAL;

		for( iLocal = 0; iLocal < vecfAssigned.size( ); ++iLocal ) {
			if( vecfAssigned[ iLocal ] )
				continue;
			vecdDiameter[ iLocal ] = max( vecdDiameter[ iLocal ], Dist.Get( iLocal, iAdded ) );
			if( vecdDiameter[ iLocal ] > dDiameter ) {
				vecfAssigned[ iLocal ] = true;
				continue; }
			if( vecdDiameter[ iLocal ] < dBest ) {
				dBest = vecdDiameter[ iLocal ];
				iBest = iLocal; } }

		if( iBest == -1 )
			break;
		vecfAssigned[ iAdded = iBest ] = true;
		vecsCluster.push_back( iAdded ); } }

void CClustQTCImpl::InitializeDistances( const CDataMatrix& Data, const IMeasure* pMeasure,
	CDistanceMatrix& Dist, const CDataMatrix* pWeights ) {
	size_t	i, j;
	float*	adA;
	float*	adB;
	float*	adWA;
	float*	adWB;

	Dist.Initialize( Data.GetRows( ) );
	adA = new float[ Data.GetColumns( ) - 1 ];
	adB = new float[ Data.GetColumns( ) - 1 ];
	if( pWeights ) {
		adWA = new float[ Data.GetColumns( ) - 1 ];
		adWB = new float[ Data.GetColumns( ) - 1 ]; }
	else
		adWA = adWB = NULL;
	for( i = 0; i < Data.GetRows( ); ++i ) {
		if( !( i % 10 ) )
			g_CatSleipnir.notice( "CClustQTCImpl::InitializeDistances( %d ) initializing %d/%d genes",
				i, Data.GetRows( ) );
		for( j = ( i + 1 ); j < Data.GetRows( ); ++j )
			Dist.Set( i, j, (float)GetJackDistance( Data.Get( i ), Data.Get( j ),
				Data.GetColumns( ), adA, adB, pMeasure, pWeights ? pWeights->Get( i ) : NULL,
				pWeights ? pWeights->Get( j ) : NULL, adWA, adWB ) ); }
	if( pWeights ) {
		delete[] adWA;
		delete[] adWB; }
	delete[] adA;
	delete[] adB; }

double CClustQTCImpl::GetJackDistance( const float* adX, const float* adY, size_t iN, float* adA, float* adB,
	const IMeasure* pMeasure, const float* adWX, const float* adWY,
	float* adWA, float* adWB ) {
	size_t	i;
	double	dRet, dCur;

	dRet = pMeasure->Measure( adX, iN, adY, iN, IMeasure::EMapCenter, adWX, adWY );
	dRet = 1 - dRet;
	for( i = 0; i < iN; ++i ) {
		memcpy( adA, adX, i * sizeof(*adX) );
		memcpy( adA + i, adX + i + 1, ( iN - 1 - i ) * sizeof(*adX) );
		memcpy( adB, adY, i * sizeof(*adY) );
		memcpy( adB + i, adY + i + 1, ( iN - 1 - i ) * sizeof(*adY) );
		if( adWX && adWY && adWA && adWB ) {
			memcpy( adWA, adWX, i * sizeof(*adWX) );
			memcpy( adWA + i, adWX + i + 1, ( iN - 1 - i ) * sizeof(*adWX) );
			memcpy( adWB, adWY, i * sizeof(*adWY) );
			memcpy( adWB + i, adWY + i + 1, ( iN - 1 - i ) * sizeof(*adWY) ); }
		dCur = pMeasure->Measure( adA, iN - 1, adB, iN - 1, IMeasure::EMapCenter, adWA, adWB );
		if( ( dCur = ( 1 - dCur ) ) > dRet )
			dRet = dCur; }

	return dRet; }

}
