#include "stdafx.h"
#include "measure.h"
#include "meta.h"
#include "statistics.h"

namespace libBioUtils {

static float GetWeight( const float* adW, size_t iW ) {

	return ( adW ? adW[ iW ] : 1 ); }

CMeasureImpl::CMeasureImpl( const IMeasure* pMeasure, bool fMemory ) : m_pMeasure((IMeasure*)pMeasure),
	m_fMemory(fMemory) { }

CMeasureImpl::~CMeasureImpl( ) {

	if( m_fMemory && m_pMeasure )
		delete m_pMeasure; }

bool CMeasureImpl::IsNaN( const float* adX, size_t iX ) {
	size_t	i;

	for( i = 0; i < iX; ++i )
		if( CMeta::IsNaN( adX[ i ] ) )
			return true;

	return false; }

double CMeasureImpl::MeasureTrim( const IMeasure* pMeasure, const float* adX, size_t iM, const float* adY,
	size_t iN, const IMeasure::EMap eMap, const float* adWX, const float* adWY ) {
	float*	adA;
	float*	adB;
	float*	adWA;
	float*	adWB;
	size_t	i, j, iA, iB;
	double	dRet;

	for( iA = i = 0; i < iM; ++i )
		if( CMeta::IsNaN( adX[ i ] ) )
			iA++;
	for( iB = i = 0; i < iN; ++i )
		if( CMeta::IsNaN( adY[ i ] ) )
			iB++;
	iA = iM - iA;
	iB = iN - iB;

	adA = new float[ iA ];
	adWA = adWX ? new float[ iA ] : NULL;
	for( i = j = 0; i < iM; ++i )
		if( !CMeta::IsNaN( adX[ i ] ) ) {
			if( adWA )
				adWA[ j ] = adWX[ i ];
			adA[ j++ ] = adX[ i ]; }
	adB = new float[ iB ];
	adWB = adWY ? new float[ iB ] : NULL;
	for( i = j = 0; i < iN; ++i )
		if( !CMeta::IsNaN( adY[ i ] ) ) {
			if( adWB )
				adWB[ j ] = adWY[ i ];
			adB[ j++ ] = adY[ i ]; }

	dRet = pMeasure->Measure( adA, iA, adB, iB, eMap, adWA, adWB );
	delete[] adA;
	delete[] adB;
	if( adWA )
		delete[] adWA;
	if( adWB )
		delete[] adWB;

	return dRet; }

const char* CMeasureKolmogorovSmirnov::GetName( ) const {

	return "kolm-smir"; }

bool CMeasureKolmogorovSmirnov::IsRank( ) const {

	return false; }

IMeasure* CMeasureKolmogorovSmirnov::Clone( ) const {

	return new CMeasureKolmogorovSmirnov( ); }

double CMeasureKolmogorovSmirnov::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	double			dCur, dMax, dRet;
	size_t			i, iX, iY;
	vector<float>	vecdX, vecdY, vecdZ;

	if( adWX || adWY )
		return CMeta::GetNaN( );
	if( CMeasureImpl::IsNaN( adX, iM ) || CMeasureImpl::IsNaN( adY, iN ) )
		return CMeasureImpl::MeasureTrim( this, adX, iM, adY, iN, eMap, adWX, adWY );
	if( iM > iN )
		return Measure( adY, iN, adX, iM, eMap, adWY, adWX );

	vecdX.resize( iM );
	vecdY.resize( iN );
	vecdZ.resize( iM + iN );
	for( i = 0; i < iM; ++i )
		vecdZ[ i ] = vecdX[ i ] = adX[ i ];
	for( i = 0; i < iN; ++i )
		vecdZ[ iM + i ] = vecdY[ i ] = adY[ i ];
	sort( vecdX.begin( ), vecdX.end( ) );
	sort( vecdY.begin( ), vecdY.end( ) );
	sort( vecdZ.begin( ), vecdZ.end( ) );

	for( dMax = iX = iY = i = 0; i < vecdZ.size( ); ++i ) {
		while( ( iX < iM ) && ( vecdX[ iX ] <= vecdZ[ i ] ) )
			iX++;
		while( ( iY < iN ) && ( vecdY[ iY ] <= vecdZ[ i ] ) )
			iY++;
		if( ( dCur = fabs( ( (double)iX / iM ) - ( (double)iY / iN ) ) ) > dMax )
			dMax = dCur; }

	dCur = sqrt( iM * iN / (double)( iM + iN ) );
	dMax *= dCur + 0.12 + ( 0.11 / dCur );
	dMax = -2 * pow( dMax, 2 );
	for( dRet = 0,i = 1; i < 100; ++i ) {
		dCur = exp( i * i * dMax );
		if( !( i % 2 ) )
			dCur *= -1;
		dRet += dCur; }
	if( dRet != 1 )
		dRet *= 2;

	return dRet; }

const char* CMeasureEuclidean::GetName( ) const {

	return "euclidean"; }

bool CMeasureEuclidean::IsRank( ) const {

	return false; }

IMeasure* CMeasureEuclidean::Clone( ) const {

	return new CMeasureEuclidean( ); }

double CMeasureEuclidean::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	size_t	i;
	double	dRet, d;

	if( iM != iN )
		return CMeta::GetNaN( );

	dRet = 0;
	for( i = 0; i < iN; ++i )
		if( ( adX[ i ] || adY[ i ] ) && !( CMeta::IsNaN( adX[ i ] ) || CMeta::IsNaN( adY[ i ] ) ) ) {
			d = adX[ i ] - adY[ i ];
			dRet += d * d * GetWeight( adWX, i ) * GetWeight( adWY, i ); }

	return sqrt( dRet ); }

const char* CMeasurePearson::GetName( ) const {

	return "pearson"; }

bool CMeasurePearson::IsRank( ) const {

	return false; }

IMeasure* CMeasurePearson::Clone( ) const {

	return new CMeasurePearson( ); }

double CMeasurePearson::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {

	return CMeasurePearson::Pearson( adX, iM, adY, iN, eMap, adWX, adWY ); }

double CMeasurePearson::Pearson( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) {
	double	dMX, dMY, dRet, dDX, dDY, dX, dY;
	size_t	i;

	if( iM != iN )
		return CMeta::GetNaN( );

	dMX = dMY = dX = dY = 0;
	for( i = 0; i < iN; ++i ) {
		if( CMeta::IsNaN( adX[ i ] ) || CMeta::IsNaN( adY[ i ] ) )
			continue;
		dX += GetWeight( adWX, i );
		dY += GetWeight( adWY, i );
		dMX += adX[ i ] * GetWeight( adWX, i );
		dMY += adY[ i ] * GetWeight( adWY, i ); }
	dMX /= dX;
	dMY /= dY;

	dRet = dDX = dDY = 0;
	for( i = 0; i < iN; ++i ) {
		if( CMeta::IsNaN( adX[ i ] ) || CMeta::IsNaN( adY[ i ] ) )
			continue;
		dX = adX[ i ] - dMX;
		dY = adY[ i ] - dMY;
		dRet += dX * dY * GetWeight( adWX, i ) * GetWeight( adWY, i );
		dDX += dX * dX * GetWeight( adWX, i );
		dDY += dY * dY * GetWeight( adWY, i ); }
	if( !( dDX || dDY ) )
		dRet = 1;
	else {
		if( dDX )
			dRet /= sqrt( dDX );
		if( dDY )
			dRet /= sqrt( dDY ); }

	switch( eMap ) {
		case EMapCenter:
			dRet = ( 1 + dRet ) / 2;
			break;

		case EMapAbs:
			dRet = fabs( dRet );
			break; }

	return dRet; }

const char* CMeasureQuickPearson::GetName( ) const {

	return "quickpear"; }

bool CMeasureQuickPearson::IsRank( ) const {

	return false; }

IMeasure* CMeasureQuickPearson::Clone( ) const {

	return new CMeasureQuickPearson( ); }

double CMeasureQuickPearson::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	double	dMX, dMY, dRet, dDX, dDY, dX, dY;
	size_t	i;

	dMX = dMY = 0;
	for( i = 0; i < iN; ++i ) {
		dMX += adX[ i ];
		dMY += adY[ i ]; }
	dMX /= iN;
	dMY /= iN;

	dRet = dDX = dDY = 0;
	for( i = 0; i < iN; ++i ) {
		dX = adX[ i ] - dMX;
		dY = adY[ i ] - dMY;
		dRet += dX * dY;
		dDX += dX * dX;
		dDY += dY * dY; }
	if( !( dDX || dDY ) )
		dRet = CMeta::GetNaN( );
	else {
		if( dDX )
			dRet /= sqrt( dDX );
		if( dDY )
			dRet /= sqrt( dDY ); }

	switch( eMap ) {
		case EMapCenter:
			dRet = ( 1 + dRet ) / 2;
			break;

		case EMapAbs:
			dRet = fabs( dRet );
			break; }

	return dRet; }

const char* CMeasureKendallsTau::GetName( ) const {

	return "kendalls"; }

bool CMeasureKendallsTau::IsRank( ) const {

	return true; }

IMeasure* CMeasureKendallsTau::Clone( ) const {

	return new CMeasureKendallsTau( ); }

double CMeasureKendallsTau::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	double	dRet;

	if( iM != iN )
		return CMeta::GetNaN( );
	if( CMeasureImpl::IsNaN( adX, iM ) || CMeasureImpl::IsNaN( adY, iN ) )
		return CMeasureImpl::MeasureTrim( this, adX, iM, adY, iN, eMap, adWX, adWY );

	dRet = ( adWX || adWY ) ? CMeasureKendallsTauImpl::MeasureWeighted( adX, adY, iN, adWX,
		adWY ) : CMeasureKendallsTauImpl::MeasureUnweighted( adX, adY, iN );
	if( dRet < -1 )
		dRet = -1;
	else if( dRet > 1 )
		dRet = 1;
	switch( eMap ) {
		case EMapCenter:
			dRet = ( 1 + dRet ) / 2;
			break;

		case EMapAbs:
			dRet = fabs( dRet );
			break; }

	return dRet; }

double CMeasureKendallsTauImpl::MeasureWeighted( const float* adX, const float* adY,
	size_t iN, const float* adWX, const float* adWY ) {
	size_t	i, j;
	double	dA1, dA2, dWX, dWY, dW, dN1, dN2, dS, dAA;

	dN1 = dN2 = dS = 0;
	for( i = 0; ( i + 1 ) < iN; ++i )
		for( j = ( i + 1 ); j < iN; ++j ) {
			dA1 = adX[ i ] - adX[ j ];
			dA2 = adY[ i ] - adY[ j ];
			dWX = GetWeight( adWX, i ) * GetWeight( adWX, j );
			dWY = GetWeight( adWY, i ) * GetWeight( adWY, j );
			dW = sqrt( dWX * dWY );
			if( dAA = ( dA1 * dA2 ) ) {
				dN1 += dWX;
				dN2 += dWY;
				dS += ( dAA > 0 ) ? dW : -dW; }
			else if( dA1 )
				dN1 += dWX;
			else
				dN2 += dWY; }

	return ( dS / ( sqrt( dN1 ) * sqrt( dN2 ) ) ); }

/*
double CMeasureKendallsTauImpl::MeasureUnweighted( const float* adX, const float* adY,
	size_t iN ) {
	size_t	i, j, iN1, iN2;
	float	dA1, dA2, dAA;
	int		iS;

	for( iN1 = iN2 = iS = i = 0; ( i + 1 ) < iN; ++i )
		for( j = ( i + 1 ); j < iN; ++j ) {
			dA1 = adX[ i ] - adX[ j ];
			dA2 = adY[ i ] - adY[ j ];
			if( dAA = ( dA1 * dA2 ) ) {
				iN1++;
				iN2++;
				iS += ( dAA > 0 ) ? 1 : -1; }
			else if( dA1 )
				iN1++;
			else
				iN2++; }

	return ( ( iN1 && iN2 ) ? ( iS / ( sqrt( (float)iN1 ) * sqrt( (float)iN2 ) ) ) : 1 ); }
*/

double CMeasureKendallsTauImpl::MeasureUnweighted( const float* adX, const float* adY,
	size_t iN ) {
	static const size_t	c_iCache	= 1024;
	static size_t		l_aiPerm[ c_iCache ];
	static size_t		l_aiTemp[ c_iCache ];
	size_t*	aiPerm;
	size_t*	aiTemp;
	size_t	i, iFirst, iT, iU, iV, iExchanges;
	int		iTotal;
	double	dBottom;

	aiTemp = ( iN > c_iCache ) ? new size_t[ iN ] : l_aiTemp;
	aiPerm = ( iN > c_iCache ) ? new size_t[ iN ] : l_aiPerm;
	for( i = 0; i < iN; ++i )
		aiPerm[ i ] = i;

	// First of all we first by the first ordering.
	sort( aiPerm, aiPerm + iN, SKendallsFirst( adX, adY ) );
	iFirst = iT = 0;
	// Next, we compute the number of joint ties.
	for( i = 1; i < iN; ++i )
		if( ( adX[ aiPerm[ iFirst ] ] != adX[ aiPerm[ i ] ] ) ||
			( adY[ aiPerm[ iFirst ] ] != adY[ aiPerm[ i ] ] ) ) {
			iT += ( ( i - iFirst ) * ( i - iFirst - 1 ) ) / 2;
			iFirst = i; }
	iT += ( ( i - iFirst ) * ( i - iFirst - 1 ) ) / 2;

	// Now we compute the number of ties.
	iFirst = iU = 0;
	for( i = 1; i < iN; ++i )
		if( adX[ aiPerm[ iFirst ] ] != adX[ aiPerm[ i ] ] ) {
			iU += ( ( i - iFirst ) * ( i - iFirst - 1 ) ) / 2;
			iFirst = i; }
	iU += ( ( i - iFirst ) * ( i - iFirst - 1 ) ) / 2;

	// Now we use an exchange counter to order by the second ordering and count the number
	// of exchanges (i.e., discordances).
	memset( aiTemp, 0, iN * sizeof(*aiTemp) );
	iExchanges = CountExchanges( aiPerm, iN, aiTemp, SKendallsSecond( adX, adY ) );

	// Now we compute the number of ties.
	iFirst = iV = 0;
	for( i = 1; i < iN; ++i )
		if( adY[ aiPerm[ iFirst ] ] != adY[ aiPerm[ i ] ] ) {
			iV += ( ( i - iFirst ) * ( i - iFirst - 1 ) ) / 2;
			iFirst = i; }
	iV += ( ( i - iFirst ) * ( i - iFirst - 1 ) ) / 2;

	if( iN > c_iCache ) {
		delete[] aiPerm;
		delete[] aiTemp; }

	iTotal = ( iN * ( iN - 1 ) ) / 2;
	dBottom = sqrt( (float)( iTotal - iU ) ) * sqrt( (float)( iTotal - iV ) );
	iTotal = ( iTotal - ( iV + iU - iT ) ) - ( 2 * iExchanges );
	return ( dBottom ? ( iTotal / dBottom ) : 1 ); }

size_t CMeasureKendallsTauImpl::CountExchanges( size_t* aiPerm, size_t iN,
	size_t* aiTemp, const SKendallsSecond& sCompare, size_t iOffset ) {
	size_t	iExchanges, iT, iL0, iL1, iMiddle, i, j, k;
	int		iD;

	if( iN == 1 )
		return 0;
	if( iN == 2 ) {
		if( sCompare( aiPerm[ iOffset ], aiPerm[ iOffset + 1 ] ) <= 0 )
			return 0;
		iT = aiPerm[ iOffset ];
		aiPerm[ iOffset ] = aiPerm[ iOffset + 1 ];
		aiPerm[ iOffset + 1 ] = iT;
		return 1; }

	iL1 = iN - ( iL0 = iN / 2 );
	iMiddle = iOffset + iL0;
	iExchanges = CountExchanges( aiPerm, iL0, aiTemp, sCompare, iOffset ) +
		CountExchanges( aiPerm, iL1, aiTemp, sCompare, iMiddle );
	
	// If the last element of the first subarray is smaller than the first element of 
	// the second subarray, there is nothing to do and we can return the exchanges got so far.
	if( sCompare( aiPerm[ iMiddle - 1 ], aiPerm[ iMiddle ] ) < 0 )
		return iExchanges;

	// We merge the lists into temp, adding the number of forward moves to exchanges.
	for( i = j = k = 0; ( j < iL0 ) || ( k < iL1 ); ++i ) {
		if( ( k >= iL1 ) || ( ( j < iL0 ) &&
			( sCompare( aiPerm[ iOffset + j ], aiPerm[ iMiddle + k ] ) <= 0 ) ) ) {
			aiTemp[ i ] = aiPerm[ iOffset + j ];
			iD = i - j++; }
		else {
			aiTemp[ i ] = aiPerm[ iMiddle + k ];
			iD = ( iOffset + i ) - ( iMiddle + k++ ); }

		if( iD > 0 )
			iExchanges += iD; }

	memcpy( aiPerm + iOffset, aiTemp, iN * sizeof(*aiPerm) );
	return iExchanges; }

CMeasureNegate::CMeasureNegate( const IMeasure* pMeasure, bool fMemory ) : CMeasureImpl( pMeasure, fMemory ) { }

const char* CMeasureNegate::GetName( ) const {

	return m_pMeasure->GetName( ); }

bool CMeasureNegate::IsRank( ) const {

	return m_pMeasure->IsRank( ); }

IMeasure* CMeasureNegate::Clone( ) const {

	return new CMeasureNegate( m_pMeasure->Clone( ), true ); }

double CMeasureNegate::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {

	return -m_pMeasure->Measure( adX, iM, adY, iN, eMap, adWX, adWY ); }

CMeasureSigmoid::CMeasureSigmoid( const IMeasure* pMeasure, bool fMemory, float dMult ) :
	CMeasureSigmoidImpl( pMeasure, fMemory, dMult ) { }

CMeasureSigmoidImpl::CMeasureSigmoidImpl( const IMeasure* pMeasure, bool fMemory, float dMult ) :
	m_dMult(dMult), CMeasureImpl( pMeasure, fMemory ) { }

const char* CMeasureSigmoid::GetName( ) const {

	return m_pMeasure->GetName( ); }

bool CMeasureSigmoid::IsRank( ) const {

	return m_pMeasure->IsRank( ); }

IMeasure* CMeasureSigmoid::Clone( ) const {

	return new CMeasureSigmoid( m_pMeasure->Clone( ), true, m_dMult ); }

double CMeasureSigmoid::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	double	dRet;

	dRet = m_pMeasure->Measure( adX, iM, adY, iN, eMap, adWX, adWY );
	return ( 2 * ( 1 - ( 1 / ( 1 + exp( -dRet * m_dMult ) ) ) ) ); }

CMeasureAutocorrelate::CMeasureAutocorrelate( const IMeasure* pMeasure, bool fMemory ) :
	CMeasureImpl( pMeasure, fMemory ) { }

const char* CMeasureAutocorrelate::GetName( ) const {

	return m_pMeasure->GetName( ); }

bool CMeasureAutocorrelate::IsRank( ) const {

	return m_pMeasure->IsRank( ); }

IMeasure* CMeasureAutocorrelate::Clone( ) const {

	return new CMeasureAutocorrelate( m_pMeasure->Clone( ), true ); }

double CMeasureAutocorrelate::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	size_t	i, j;
	double	dCur, dMax;
	float*	adZ;
	float*	adWZ;

	if( iM != iN )
		return CMeta::GetNaN( );

	dMax = m_pMeasure->Measure( adX, iM, adY, iN, eMap, adWX, adWY );
	adZ = new float[ iN ];
	adWZ = adWY ? new float[ iN ] : NULL;
	for( i = 1; i < iN; ++i ) {
		for( j = 0; j < iN; ++j ) {
			adZ[ j ] = adY[ ( j + i ) % iN ];
			if( adWZ )
				adWZ[ j ] = adWY[ ( j + i ) % iN ]; }
		if( ( dCur = m_pMeasure->Measure( adX, iM, adZ, iN, eMap, adWX, adWZ ) ) > dMax )
			dMax = dCur; }
	delete[] adZ;

	return dMax; }

CMeasureSpearman::CMeasureSpearman( bool fTransformed ) :
	CMeasureSpearmanImpl( fTransformed ) { }

CMeasureSpearmanImpl::CMeasureSpearmanImpl( bool fTransformed ) :
	m_fTransformed(fTransformed) { }

const char* CMeasureSpearman::GetName( ) const {

	return "spearman"; }

bool CMeasureSpearman::IsRank( ) const {

	return true; }

IMeasure* CMeasureSpearman::Clone( ) const {

	return new CMeasureSpearman( m_fTransformed ); }

double CMeasureSpearman::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	static const size_t	c_iCache	= 1024;
	static size_t		l_aiX[ c_iCache ];
	static size_t		l_aiY[ c_iCache ];
	size_t*	aiX;
	size_t*	aiY;
	size_t	i, j, iSum;
	double	dRet, d, dSum;

	if( ( iM != iN ) || adWX || adWY )
		return CMeta::GetNaN( );
	if( CMeasureImpl::IsNaN( adX, iM ) || CMeasureImpl::IsNaN( adY, iN ) )
		return CMeasureImpl::MeasureTrim( this, adX, iM, adY, iN, eMap, adWX, adWY );

	if( m_fTransformed ) {
		dSum = 0;
		for( i = 0; i < iN; ++i ) {
			d = adX[ i ] - adY[ i ];
			dSum += d * d; }
		dRet = dSum ? 1 - ( 6 * dSum / iN / ( ( iN * iN ) - 1 ) ) : 1; }
	else {
		if( iN > c_iCache ) {
			aiX = new size_t[ iM ];
			aiY = new size_t[ iN ]; }
		else {
			aiX = l_aiX;
			aiY = l_aiY; }
		memset( aiX, 0, iM * sizeof(*aiX) );
		memset( aiY, 0, iN * sizeof(*aiY) );
		for( i = 0; i < iN; ++i )
			for( j = 0; j < iN; ++j ) {
				if( i == j )
					continue;
				if( adX[ j ] < adX[ i ] )
					aiX[ i ]++;
				if( adY[ j ] < adY[ i ] )
					aiY[ i ]++; }

		for( iSum = i = 0; i < iN; ++i ) {
			j = aiX[ i ] - aiY[ i ];
			iSum += j * j; }
		if( aiX != l_aiX )
			delete[] aiX;
		if( aiY != l_aiY )
			delete[] aiY;
		dRet = iSum ? 1 - ( 6.0 * iSum / iN / ( ( iN * iN ) - 1 ) ) : 1; }

	switch( eMap ) {
		case EMapCenter:
			dRet = ( 1 + dRet ) / 2;
			break;

		case EMapAbs:
			dRet = fabs( dRet );
			break; }

	return dRet; }

CMeasurePearNorm::CMeasurePearNorm( ) : CMeasurePearNormImpl( HUGE_VAL, HUGE_VAL ) { }

CMeasurePearNorm::CMeasurePearNorm( double dAve, double dStd ) : CMeasurePearNormImpl( dAve, dStd ) { }

CMeasurePearNormImpl::CMeasurePearNormImpl( double dAve, double dStd ) : m_dAverage(dAve), m_dStdDev(dStd) { }

const char* CMeasurePearNorm::GetName( ) const {

	return "pearnorm"; }

bool CMeasurePearNorm::IsRank( ) const {

	return false; }

IMeasure* CMeasurePearNorm::Clone( ) const {

	return new CMeasurePearNorm( m_dAverage, m_dStdDev ); }

double CMeasurePearNorm::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	static const float	c_dBound	= 0.99f;
	double	dP;

	dP = CMeasurePearson::Pearson( adX, iM, adY, iN, EMapNone, adWX, adWY );
	if( fabs( dP ) >= c_dBound )
		dP *= c_dBound;
	dP = log( ( 1 + dP ) / ( 1 - dP ) ) / 2;
	if( m_dAverage != HUGE_VAL )
		dP = ( dP - m_dAverage ) / m_dStdDev;
	return dP; }

const char* CMeasureHypergeometric::GetName( ) const {

	return "hypergeom"; }

bool CMeasureHypergeometric::IsRank( ) const {

	return false; }

IMeasure* CMeasureHypergeometric::Clone( ) const {

	return new CMeasureHypergeometric( ); }

double CMeasureHypergeometric::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	size_t	i, iOne, iTwo, iBoth;

	if( iM != iN )
		return CMeta::GetNaN( );

	iOne = iTwo = iBoth = 0;
	for( i = 0; i < iN; ++i ) {
		if( CMeta::IsNaN( adX[ i ] ) || CMeta::IsNaN( adY[ i ] ) )
			continue;
		if( adX[ i ] )
			iOne++;
		if( adY[ i ] ) {
			iTwo++;
			if( adX[ i ] )
				iBoth++; } }

	return ( 1 - CStatistics::HypergeometricCDF( iBoth, iOne, iTwo, iN ) ); }

const char* CMeasureInnerProduct::GetName( ) const {

	return "innerprod"; }

bool CMeasureInnerProduct::IsRank( ) const {

	return false; }

IMeasure* CMeasureInnerProduct::Clone( ) const {

	return new CMeasureInnerProduct( ); }

double CMeasureInnerProduct::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	size_t	i;
	double	dRet;

	if( iM != iN )
		return CMeta::GetNaN( );

	dRet = 0;
	for( i = 0; i < iN; ++i )
		if( ( adX[ i ] || adY[ i ] ) && !( CMeta::IsNaN( adX[ i ] ) || CMeta::IsNaN( adY[ i ] ) ) )
			dRet += adX[ i ] * adY[ i ] * GetWeight( adWX, i ) * GetWeight( adWY, i );

	return dRet; }

const char* CMeasureBinaryInnerProduct::GetName( ) const {

	return "bininnerprod"; }

bool CMeasureBinaryInnerProduct::IsRank( ) const {

	return false; }

IMeasure* CMeasureBinaryInnerProduct::Clone( ) const {

	return new CMeasureBinaryInnerProduct( ); }

double CMeasureBinaryInnerProduct::Measure( const float* adX, size_t iM, const float* adY,
	size_t iN, EMap eMap, const float* adWX, const float* adWY ) const {
	size_t	i;
	double	dRet;

	if( iM != iN )
		return CMeta::GetNaN( );

	dRet = 0;
	for( i = 0; i < iN; ++i ) {
		if( CMeta::IsNaN( adX[ i ] ) || CMeta::IsNaN( adY[ i ] ) || !adX[ i ] || !adY[ i ] )
			continue;
		dRet += GetWeight( adWX, i ) * GetWeight( adWY, i ); }

	return dRet; }

}
