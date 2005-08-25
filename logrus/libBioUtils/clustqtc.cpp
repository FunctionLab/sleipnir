#include "stdafx.h"
#include "clustqtc.h"
#include "pcl.h"
#include "measure.h"

namespace libBioUtils {

size_t CClustQTC::Cluster( const CDataMatrix& Data, const IMeasure* pMeasure,
	float dDiameter, size_t iSize, bool fAutoc, vector<size_t>& veciClusters ) {
	CDistanceMatrix	Dist;

	InitializeDistances( Data, pMeasure, fAutoc, Dist );
	return QualityThresholdAll( Data, 2 * dDiameter, iSize, Dist, veciClusters ); }

size_t CClustQTCImpl::QualityThresholdAll( const CDataMatrix& Data, float dDiameter,
	size_t iSize, const CDistanceMatrix& Dist, vector<size_t>& veciClusters ) {
	size_t			i, iCluster, iAssigned;
	vector<size_t>	veciCur;
	vector<bool>	vecfAssigned;

	vecfAssigned.resize( Data.GetRows( ) );
	for( i = 0; i < vecfAssigned.size( ); ++i )
		vecfAssigned[ i ] = false;

	veciClusters.resize( Data.GetRows( ) );
	for( iAssigned = iCluster = 0; ; ++iCluster ) {
		g_CatBioUtils.notice( "CClustQTCImpl::QualityThresholdAll( ) cluster %d, assigned %d/%d genes",
			iCluster + 1, iAssigned, Data.GetRows( ) );
		QualityThresholdLargest( Data, dDiameter, Dist, vecfAssigned, veciCur );
		for( i = 0; i < veciCur.size( ); ++i )
			veciClusters[ veciCur[ i ] ] = iCluster;

		if( veciCur.size( ) < iSize ) {
			for( i = 0; i < vecfAssigned.size( ); ++i )
				if( !vecfAssigned[ i ] )
					veciClusters[ i ] = iCluster;
			break; }

		iAssigned += veciCur.size( );
		for( i = 0; i < veciCur.size( ); ++i )
			vecfAssigned[ veciCur[ i ] ] = true;

		for( i = 0; i < vecfAssigned.size( ); ++i )
			if( vecfAssigned[ i ] )
				break;
		if( i >= vecfAssigned.size( ) ) {
			iCluster++;
			break; } }

	return ( iCluster + 1 ); }

void CClustQTCImpl::QualityThresholdLargest( const CDataMatrix& Data, float dDiameter,
	const CDistanceMatrix& Dist, const vector<bool>& vecfAssigned,
	vector<size_t>& veciCluster ) {
	vector<bool>	vecfClone;
	size_t			i, iGene;
	vector<size_t>	veciCur;
	vector<float>	vecdDiameter;

	veciCluster.clear( );
	vecfClone.resize( vecfAssigned.size( ) );
	for( iGene = 0; iGene < vecfAssigned.size( ); ++iGene ) {
		if( vecfAssigned[ iGene ] )
			continue;
		for( i = 0; i < vecfClone.size( ); ++i )
			vecfClone[ i ] = vecfAssigned[ i ];
		vecfClone[ iGene ] = true;

		QualityThresholdGene( iGene, Data, dDiameter, Dist, vecfClone, vecdDiameter,
			veciCur );

		if( veciCur.size( ) > veciCluster.size( ) ) {
			veciCluster.resize( veciCur.size( ) );
			for( i = 0; i < veciCluster.size( ); ++i )
				veciCluster[ i ] = veciCur[ i ]; } } }

void CClustQTCImpl::QualityThresholdGene( size_t iGene, const CDataMatrix& Data,
	float dDiameter, const CDistanceMatrix& Dist, vector<bool>& vecfAssigned,
	vector<float>& vecdDiameter,
	vector<size_t>& veciCluster ) {
	size_t	iAdded, iLocal, iBest;
	float	dBest;

	veciCluster.resize( 1 );
	veciCluster[ 0 ] = iAdded = iGene;
	vecdDiameter.resize( Data.GetRows( ) );
	for( iLocal = 0; iLocal < vecfAssigned.size( ); ++iLocal ) {
		if( vecfAssigned[ iLocal ] )
			continue;
		vecdDiameter[ iLocal ] = -(float)HUGE_VAL; }

	while( true ) {
		iBest = -1;
		dBest = (float)HUGE_VAL;

		for( iLocal = 0; iLocal < vecfAssigned.size( ); ++iLocal ) {
			if( vecfAssigned[ iLocal ] )
				continue;

			vecdDiameter[ iLocal ] = max( vecdDiameter[ iLocal ],
				Dist.Get( iLocal, iAdded ) );
			if( vecdDiameter[ iLocal ] > dDiameter ) {
				vecfAssigned[ iLocal ] = true;
				continue; }
			if( vecdDiameter[ iLocal ] < dBest ) {
				dBest = vecdDiameter[ iLocal ];
				iBest = iLocal; } }

		if( iBest == -1 )
			break;
		vecfAssigned[ iAdded = iBest ] = true;
		veciCluster.push_back( iAdded ); } }

void CClustQTCImpl::InitializeDistances( const CDataMatrix& Data, const IMeasure* pMeasure,
	bool fAutoc, CDistanceMatrix& Dist ) {
	size_t	i, j;
	float*	adA;
	float*	adB;

	Dist.Initialize( Data.GetRows( ) );
	adA = new float[ Data.GetColumns( ) - 1 ];
	adB = new float[ Data.GetColumns( ) - 1 ];
	for( i = 0; i < Data.GetRows( ); ++i )
		for( j = ( i + 1 ); j < Data.GetRows( ); ++j )
			Dist.Set( i, j, (float)GetJackDistance( Data.Get( i ), Data.Get( j ),
				Data.GetColumns( ), fAutoc, adA, adB, pMeasure ) );
	delete[] adA;
	delete[] adB; }

double CClustQTCImpl::GetJackDistance( const float* adX, const float* adY, size_t iN,
	bool fAutoc, float* adA, float* adB, const IMeasure* pMeasure ) {
	size_t					i, j, k;
	double					dRet, dCur;
	CMeasureAutocorrelate	Autocorrelate( pMeasure );
	const IMeasure*			pAuto	= &Autocorrelate;

	dRet = fAutoc ? pAuto->Measure( adX, iN, adY, iN ) :
		pMeasure->Measure( adX, iN, adY, iN );
	dRet = 1 - dRet;

	for( i = 0; i < iN; ++i ) {
		for( j = k = 0; j < iN; ++j ) {
			if( j == i )
				continue;
			adA[ k ] = adX[ j ];
			adB[ k++ ] = adY[ j ]; }
		dCur = fAutoc ? pAuto->Measure( adA, iN - 1, adB, iN - 1 ) :
			pMeasure->Measure( adA, iN - 1, adB, iN - 1 );
		if( ( dCur = ( 1 - dCur ) ) > dRet )
			dRet = dCur; }

	return dRet; }

}
