#ifndef STATISTICS_H
#define STATISTICS_H

#include <algorithm>

#include "statisticsi.h"

#include "dat.h"
#include "halfmatrix.h"
#include "meta.h"

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

private:
	template<class tType>
	struct SCompareRank {
		const std::vector<tType>&	m_vecData;

		SCompareRank( const std::vector<tType>& vecData ) : m_vecData(vecData) { }

		bool operator()( size_t iOne, size_t iTwo ) const {

			return ( m_vecData[ iOne ] < m_vecData[ iTwo ] ); }
	};

public:
	static double WilcoxonRankSum( const CDat& Data, const CDat& Answers, const std::vector<bool>& vecfGenes,
		const CBinaryMatrix& MatPairs, bool fInvert = false ) {
		size_t				i, j, k, iOne, iTwo, iSum, iPos, iNeg;
		float				d, dAnswer;
		std::vector<size_t>	veciGenes, veciRanks;
		std::vector<float>	vecdValues;

		veciGenes.resize( Answers.GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
#pragma warning( disable : 4267 )
			veciGenes[ i ] = Data.GetGene( Answers.GetGene( i ) );
#pragma warning( default : 4267 )

		for( i = 0; i < Answers.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
				if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
					CMeta::IsNaN( d = Data.Get( iOne, iTwo ) ) ||
					CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) )
					continue;
				if( !vecfGenes.empty( ) && !MatPairs.Get( i, j ) &&
					( ( dAnswer && !( vecfGenes[ i ] && vecfGenes[ j ] ) ) ||
					( !dAnswer && !( vecfGenes[ i ] || vecfGenes[ j ] ) ) ) )
					continue;
				if( fInvert )
					d = 1 - d;
				vecdValues.push_back( d ); } }
		{
			std::vector<size_t>	veciIndices;

			veciIndices.resize( vecdValues.size( ) );
			for( i = 0; i < vecdValues.size( ); ++i )
#pragma warning( disable : 4267 )
				veciIndices[ i ] = i;
#pragma warning( default : 4267 )
			std::sort( veciIndices.begin( ), veciIndices.end( ), SCompareRank<float>( vecdValues ) );
			veciRanks.resize( veciIndices.size( ) );
			for( i = 0; i < veciRanks.size( ); ++i )
#pragma warning( disable : 4267 )
				veciRanks[ veciIndices[ i ] ] = i;
#pragma warning( default : 4267 )
		}

		for( iPos = iNeg = iSum = i = k = 0; i < Answers.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
				if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
					CMeta::IsNaN( Data.Get( iOne, iTwo ) ) ||
					CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) ||
					( !vecfGenes.empty( ) && !MatPairs.Get( i, j ) &&
					( ( dAnswer && !( vecfGenes[ i ] && vecfGenes[ j ] ) ) ||
					( !dAnswer && !( vecfGenes[ i ] || vecfGenes[ j ] ) ) ) ) )
					continue;
				if( dAnswer ) {
					iPos++;
					iSum += veciRanks[ k ]; }
				else
					iNeg++;
				k++; } }
		iSum -= ( iPos * ( iPos - 1 ) ) / 2;

		return ( (double)iSum / ( iPos * iNeg ) ); }

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
