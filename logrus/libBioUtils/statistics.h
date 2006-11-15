#ifndef STATISTICS_H
#define STATISTICS_H

#include "statisticsi.h"

namespace libBioUtils {

class CBinaryMatrix;
class CDat;

class CStatistics : CStatisticsImpl {
public:
	// Simple stuff
	static double FScore( size_t, size_t, size_t, size_t, double = 1 );
	static double Precision( size_t, size_t, size_t, size_t );
	static double Recall( size_t, size_t, size_t, size_t );

	template<class tType>
	static double Variance( tType Begin, tType End, double dMean ) {
		tType	Cur;
		size_t	iN;
		double	d, dRet;

		for( dRet = iN = 0,Cur = Begin; Cur != End; ++Cur )
			if( !CMeta::IsNaN( *Cur ) ) {
				iN++;
				d = *Cur - dMean;
				dRet += d * d; }

		return ( ( iN < 2 ) ? 0 : ( dRet / ( iN - 1 ) ) ); }

	template<class tType>
	static double Variance( const std::vector<tType>& vecValues, double dMean ) {

		return Variance( vecValues.begin( ), vecValues.end( ), dMean ); }

	template<class tType>
	static double Variance( const std::vector<tType>& vecValues ) {

		return Variance( vecValues, Average( vecValues ) ); }

	template<class tType>
	static double Average( const std::vector<tType>& vecValues ) {

		return Average( vecValues.begin( ), vecValues.end( ) ); }

	template<class tType>
	static double Average( tType Begin, tType End ) {
		tType	Cur;
		double	dRet;
		size_t	iN;

		for( dRet = iN = 0,Cur = Begin; Cur != End; ++Cur )
			if( !CMeta::IsNaN( *Cur ) ) {
				iN++;
				dRet += *Cur; }

		return ( dRet / iN ); }

	// P-value tests
	static double LjungBox( const float*, size_t, size_t );
	static double LjungBox( const float*, size_t );
	static double PValueLognormal( double, double, double );
	static double TTest( double, double, size_t, double, double, size_t );

private:
	template<class tType>
	struct SCompareRank {
		const std::vector<tType>&	m_vecData;

		SCompareRank( const std::vector<tType>& vecData ) : m_vecData(vecData) { }

		bool operator()( size_t iOne, size_t iTwo ) const {

			return ( m_vecData[ iOne ] < m_vecData[ iTwo ] ); }
	};

public:
	// Evaluation statistics
	static double WilcoxonRankSum( const CDat& Data, const CDat& Answers, const std::vector<bool>& vecfGenes,
		const CBinaryMatrix& MatPairs, bool fInvert = false );
	static double SumSquaredError( const CDat& Data, const CDat& Answers, const std::vector<bool>& vecfGenes,
		const CBinaryMatrix& MatPairs, bool fInvert = false );

	// Probability distributions
	static double BinomialCDF( size_t, size_t, double );
	static double LognormalCDF( double, double, double );
	static double NormalCDF( double, double, double );
	static double HypergeometricPDF( size_t, size_t, size_t, size_t );
	static double HypergeometricCDF( size_t, size_t, size_t, size_t );
	static double SampleChi2( size_t );
	static double SampleGamma( double, double );
	static double SampleGammaStandard( double );
	static double SampleGammaLogStandard( double );
	static double SampleNormalStandard( );
	static double SampleExponentialStandard( );

	static double BetaPDF( double dX, double dMin, double dMax, double dAlpha, double dBeta ) {
		double	dFunc, dLocation, dScale;

		dLocation = dMin;
		dScale = dMax - dMin;
		dX = ( dX - dLocation ) / dScale;
		dFunc = exp( SampleGammaLogStandard( dAlpha ) + SampleGammaLogStandard( dBeta ) -
			SampleGammaLogStandard( dAlpha + dBeta ) );
		return ( pow( dX, dAlpha - 1 ) * pow( 1 - dX, dBeta - 1 ) / dFunc / dScale ); }

	static double NormalPDF( double dX, double dMu, double dSigma ) {
		static const double	c_dS2P	= sqrt( 2 * 3.1415926535898 );
		double	d;

		d = dX - dMu;
		return ( exp( -( d * d ) / ( 2 * dSigma * dSigma ) ) / ( dSigma * c_dS2P ) ); }
};

}

#endif // STATISTICS_H
