#ifndef STATISTICS_H
#define STATISTICS_H

#include "statisticsi.h"

namespace libBioUtils {

class CStatistics : CStatisticsImpl {
public:
	// Simple stuff
	static double FScore( size_t, size_t, size_t, size_t, double = 1 );
	static double Precision( size_t, size_t, size_t, size_t );
	static double Recall( size_t, size_t, size_t, size_t );

	template<class tType>
	static double Variance( const std::vector<tType>& vecValues, double dMean ) {
		size_t	i;
		double	d, dRet;

		if( vecValues.size( ) < 2 )
			return 0;

		dRet = 0;
		for( i = 0; i < vecValues.size( ); ++i ) {
			d = vecValues[ i ] - dMean;
			dRet += d * d; }

		return ( dRet / ( i - 1 ) ); }

	template<class tType>
	static double Variance( const std::vector<tType>& vecValues ) {

		return Variance( vecValues, Average( vecValues ) ); }

	template<class tType>
	static double Average( const std::vector<tType>& vecValues ) {
		size_t	i;
		double	dRet;

		dRet = 0;
		for( i = 0; i < vecValues.size( ); ++i )
			dRet += vecValues[ i ];

		return ( dRet / i ); }

	// P-value tests
	static double LjungBox( const float*, size_t, size_t );
	static double LjungBox( const float*, size_t );
	static double PValueLognormal( double, double, double );
	static double TTest( double, double, size_t, double, double, size_t );

	// Probability distributions
	static double BinomialCDF( size_t, size_t, double );
	static double LognormalCDF( double, double, double );
	static double NormalCDF( double, double, double );
	static double HypergeometricPDF( size_t, size_t, size_t, size_t );
	static double HypergeometricCDF( size_t, size_t, size_t, size_t );
	static double SampleChi2( size_t );
	static double SampleGamma( double, double );
	static double SampleGammaStandard( double );
	static double SampleNormalStandard( );
	static double SampleExponentialStandard( );
};

}

#endif // STATISTICS_H
