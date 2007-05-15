#include "stdafx.h"
#include "clustqtc.h"
#include "pcl.h"
#include "measure.h"
#include "meta.h"

namespace libBioUtils {

uint16_t CClustQTC::Cluster( const CDataMatrix& Data, const IMeasure* pMeasure,
	float dDiameter, size_t iSize, bool fAutoc, vector<uint16_t>& vecsClusters, const CDataMatrix* pWeights ) {
	CDistanceMatrix	Dist;

	InitializeDistances( Data, pMeasure, fAutoc, Dist, pWeights );
	return QualityThresholdAll( Data, dDiameter, iSize, Dist, vecsClusters ); }

void CClustQTC::Cluster( const CDataMatrix& Data, const IMeasure* pMeasure,
	float dMin, float dMax, float dDelta, size_t iSize, bool fAutoc, CDistanceMatrix& DatOut,
	const CDataMatrix* pWeights ) {
	CDistanceMatrix		Dist;
	float				dDiameter;
	vector<uint16_t>	vecsClusters;
	uint16_t			sClusters;
	size_t				i, j;

	InitializeDistances( Data, pMeasure, fAutoc, Dist, pWeights );
	for( dDiameter = dMin; dDiameter <= dMax; dDiameter += dDelta ) {
		g_CatBioUtils.notice( "CClustQTC::Cluster( %g, %g, %g, %d, %d ) processing diameter %g", dMin, dMax,
			dDelta, iSize, fAutoc, dDiameter );
		sClusters = QualityThresholdAll( Data, dDiameter, iSize, Dist, vecsClusters );
		for( i = 0; i < vecsClusters.size( ); ++i ) {
			if( ( vecsClusters[ i ] + 1 ) == sClusters )
				continue;
			for( j = ( i + 1 ); j < vecsClusters.size( ); ++j )
				if( ( vecsClusters[ j ] == vecsClusters[ i ] ) && CMeta::IsNaN( DatOut.Get( i, j ) ) )
					DatOut.Set( i, j, 1 - dDiameter ); } } }

uint16_t CClustQTCImpl::QualityThresholdAll( const CDataMatrix& Data, float dDiameter,
	size_t iSize, const CDistanceMatrix& Dist, vector<uint16_t>& vecsClusters ) {
	size_t				i, iAssigned;
	uint16_t			sCluster;
	vector<uint16_t>	vecsCur;
	vector<bool>		vecfAssigned;

	vecfAssigned.resize( Data.GetRows( ) );
	for( i = 0; i < vecfAssigned.size( ); ++i )
		vecfAssigned[ i ] = false;

	vecsClusters.resize( Data.GetRows( ) );
	for( iAssigned = sCluster = 0; ; ++sCluster ) {
		g_CatBioUtils.notice( "CClustQTCImpl::QualityThresholdAll( ) cluster %d, assigned %d/%d genes",
			sCluster + 1, iAssigned, Data.GetRows( ) );
		QualityThresholdLargest( Data, dDiameter, Dist, vecfAssigned, vecsCur );
		for( i = 0; i < vecsCur.size( ); ++i )
			vecsClusters[ vecsCur[ i ] ] = sCluster;

		if( vecsCur.size( ) < iSize ) {
			for( i = 0; i < vecfAssigned.size( ); ++i )
				if( !vecfAssigned[ i ] )
					vecsClusters[ i ] = sCluster;
			break; }

		if( ( iAssigned += vecsCur.size( ) ) >= Data.GetRows( ) ) {
			sCluster++;
			break; }
		for( i = 0; i < vecsCur.size( ); ++i )
			vecfAssigned[ vecsCur[ i ] ] = true; }

	return ( sCluster + 1 ); }

void CClustQTCImpl::QualityThresholdLargest( const CDataMatrix& Data, float dDiameter,
	const CDistanceMatrix& Dist, const vector<bool>& vecfAssigned,
	vector<uint16_t>& vecsCluster ) {
	vector<bool>		vecfClone;
	size_t				iGene;
	vector<uint16_t>	vecsCur;
	vector<float>		vecdDiameter;

	vecsCluster.clear( );
	vecfClone.resize( vecfAssigned.size( ) );
	for( iGene = 0; iGene < vecfAssigned.size( ); ++iGene ) {
		if( !( iGene % 1000 ) )
			g_CatBioUtils.notice( "CClustQTCImpl::QualityThresholdLargest( %g ) processing gene %d/%d",
				dDiameter, iGene, vecfAssigned.size( ) );
		if( vecfAssigned[ iGene ] )
			continue;
		copy( vecfAssigned.begin( ), vecfAssigned.end( ), vecfClone.begin( ) );
		vecfClone[ iGene ] = true;
		QualityThresholdGene( iGene, Data, dDiameter, Dist, vecfClone, vecdDiameter, vecsCur );
		if( vecsCur.size( ) > vecsCluster.size( ) ) {
			vecsCluster.resize( vecsCur.size( ) );
			copy( vecsCur.begin( ), vecsCur.end( ), vecsCluster.begin( ) ); } } }

void CClustQTCImpl::QualityThresholdGene( size_t iGene, const CDataMatrix& Data,
	float dDiameter, const CDistanceMatrix& Dist, vector<bool>& vecfAssigned,
	vector<float>& vecdDiameter, vector<uint16_t>& vecsCluster ) {
	size_t	iAdded, iLocal, iBest;
	float	dBest;

	vecsCluster.resize( 1 );
	vecsCluster[ 0 ] = iAdded = iGene;
	vecdDiameter.resize( Data.GetRows( ) );
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
	bool fAutoc, CDistanceMatrix& Dist, const CDataMatrix* pWeights ) {
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
			g_CatBioUtils.notice( "CClustQTCImpl::InitializeDistances( %d ) initializing %d/%d genes",
				fAutoc, i, Data.GetRows( ) );
		for( j = ( i + 1 ); j < Data.GetRows( ); ++j )
			Dist.Set( i, j, (float)GetJackDistance( Data.Get( i ), Data.Get( j ),
				Data.GetColumns( ), fAutoc, adA, adB, pMeasure, pWeights ? pWeights->Get( i ) : NULL,
				pWeights ? pWeights->Get( j ) : NULL, adWA, adWB ) ); }
	if( pWeights ) {
		delete[] adWA;
		delete[] adWB; }
	delete[] adA;
	delete[] adB; }

double CClustQTCImpl::GetJackDistance( const float* adX, const float* adY, size_t iN,
	bool fAutoc, float* adA, float* adB, const IMeasure* pMeasure, const float* adWX, const float* adWY,
	float* adWA, float* adWB ) {
	size_t					i;
	double					dRet, dCur;
	CMeasureAutocorrelate	Autocorrelate( pMeasure, false );
	const IMeasure*			pUse	= fAutoc ? &Autocorrelate : pMeasure;

	dRet = pUse->Measure( adX, iN, adY, iN, IMeasure::EMapCenter, adWX, adWY );
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
		dCur = pUse->Measure( adA, iN - 1, adB, iN - 1, IMeasure::EMapCenter, adWA, adWB );
		if( ( dCur = ( 1 - dCur ) ) > dRet )
			dRet = dCur; }

	return dRet; }

}
