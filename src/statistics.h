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
#ifndef STATISTICS_H
#define STATISTICS_H

#include "meta.h"

#include "statisticsi.h"

namespace Sleipnir {

class CDat;

/*!
 * \brief
 * Utility class containing static statistics functions.
 * 
 * CStatistics contains a collection of static utility functions, mainly various statistical distributions.
 * Credit for many of these implementations goes to Press WH, Teukolsky SA, Vetterling WT, Flannery BP.
 * Numerical Recipes in C, 1992, Cambridge University Press.
 */
class CStatistics : CStatisticsImpl {
public:
	// Simple stuff
	static double FScore( size_t iTruePositives, size_t iFalsePositives, size_t iTrueNegatives,
		size_t iFalseNegatives, double dBeta = 1 );

	/*!
	 * \brief
	 * Calculate the precision of a predicted positive set.
	 * 
	 * \param iTruePositives
	 * Number of true positives.
	 * 
	 * \param iFalsePositives
	 * Number of false positives.
	 * 
	 * \param iTrueNegatives
	 * Number of true negatives.
	 * 
	 * \param iFalseNegatives
	 * Number of false negatives.
	 * 
	 * \returns
	 * iTruePositives / (iTruePositives + iFalsePositives)
	 * 
	 * \remarks
	 * True and false negatives are ignored and included only for consistency.
	 * 
	 * \see
	 * Recall | FScore
	 */
	static double Precision( size_t iTruePositives, size_t iFalsePositives, size_t iTrueNegatives,
		size_t iFalseNegatives ) {
		UNUSED_PARAMETER(iTrueNegatives);
		UNUSED_PARAMETER(iFalseNegatives);

		return ( (double)iTruePositives / ( iTruePositives + iFalsePositives ) ); }

	/*!
	 * \brief
	 * Calculate the recall of predictions over some positive set.
	 * 
	 * \param iTruePositives
	 * Number of true positives.
	 * 
	 * \param iFalsePositives
	 * Number of false positives.
	 * 
	 * \param iTrueNegatives
	 * Number of true negatives.
	 * 
	 * \param iFalseNegatives
	 * Number of false negatives.
	 * 
	 * \returns
	 * iTruePositives / (iTruePositives + iFalseNegatives)
	 * 
	 * \remarks
	 * False positives and true negatives are ignored and included only for consistency.
	 * 
	 * \see
	 * Precision | FScore
	 */
	static double Recall( size_t iTruePositives, size_t iFalsePositives, size_t iTrueNegatives,
		size_t iFalseNegatives ) {
		UNUSED_PARAMETER(iFalsePositives);
		UNUSED_PARAMETER(iTrueNegatives);

		return ( (double)iTruePositives / ( iTruePositives + iFalseNegatives ) ); }

	/*!
	 * \brief
	 * Calculate the sample variance of a given sequence.
	 * 
	 * \param Begin
	 * Beginning of the sequence over which variance is estimated.
	 * 
	 * \param End
	 * End of the sequence over which variance is estimated.
	 * 
	 * \param dMean
	 * Mean of the sample or population.
	 * 
	 * \returns
	 * sum((x - dMean)^2) / (n - 1)
	 * 
	 * \see
	 * Average
	 */
	template<class tType>
	static double Variance( const tType Begin, const tType End, double dMean ) {
		tType	Cur;
		size_t	iN;
		double	d, dRet;

		for( dRet = iN = 0,Cur = Begin; Cur != End; ++Cur )
			if( !CMeta::IsNaN( *Cur ) ) {
				iN++;
				d = *Cur - dMean;
				dRet += d * d; }

		return ( ( iN < 2 ) ? 0 : ( dRet / ( iN - 1 ) ) ); }

	/*!
	 * \brief
	 * Calculate the sample variance of a given vector.
	 * 
	 * \param vecValues
	 * Vector over which variance is estimated.
	 * 
	 * \param dMean
	 * Mean of the sample or population.
	 * 
	 * \returns
	 * sum((x - dMean)^2) / (n - 1)
	 * 
	 * \see
	 * Average
	 */
	template<class tType>
	static double Variance( const std::vector<tType>& vecValues, double dMean ) {

		return Variance( vecValues.begin( ), vecValues.end( ), dMean ); }

	/*!
	 * \brief
	 * Calculate the sample variance of a given vector.
	 * 
	 * \param vecValues
	 * Vector over which variance is estimated.
	 * 
	 * \returns
	 * sum((x - xbar)^2) / (n - 1)
	 * 
	 * \see
	 * Average
	 */
	template<class tType>
	static double Variance( const std::vector<tType>& vecValues ) {

		return Variance( vecValues, Average( vecValues ) ); }

	/*!
	 * \brief
	 * Calculate the mean of a given vector.
	 * 
	 * \param vecValues
	 * Vector over which average is calculated.
	 * 
	 * \returns
	 * sum(x) / n
	 * 
	 * \see
	 * Variance
	 */
	template<class tType>
	static double Average( const std::vector<tType>& vecValues ) {

		return Average( vecValues.begin( ), vecValues.end( ) ); }

	/*!
	 * \brief
	 * Calculate the mean of a given sequence.
	 * 
	 * \param Begin
	 * Beginning of the sequence over which average is calculated.
	 * 
	 * \param End
	 * End of the sequence over which average is calculated.
	 * 
	 * \returns
	 * sum(x) / n
	 * 
	 * \see
	 * Variance
	 */
	template<class tType>
	static double Average( const tType Begin, const tType End ) {
		tType	Cur;
		double	dRet;
		size_t	iN;

		for( dRet = iN = 0,Cur = Begin; Cur != End; ++Cur )
			if( !CMeta::IsNaN( *Cur ) ) {
				iN++;
				dRet += *Cur; }

		return ( dRet / iN ); }

	/*!
	 * \brief
	 * Returns the value at the requested percentile of the given sequence.
	 * 
	 * \param pBegin
	 * Pointer/iterator to the beginning of the sequence.
	 * 
	 * \param pEnd
	 * Pointer/iterator to the end of the sequence.
	 * 
	 * \param dPercentile
	 * Value between 0 and 1 indicating the percentile of the value to be retrieved.
	 * 
	 * \returns
	 * Value at the requested percentile.
	 * 
	 * \remarks
	 * The 0th percentile is the first (smallest) element, 1st percentile is the last (largest), and 0.5th is
	 * the median.  Note that the input sequence is modified (sorted) by this function.
	 * 
	 * \see
	 * Median
	 */
	template<class tType>
	static double Percentile( tType pBegin, tType pEnd, double dPercentile ) {
		size_t	iOne, iTwo, iSize;
		double	d, dFrac;

		iSize = pEnd - pBegin;
		if( !iSize )
			return CMeta::GetNaN( );
		sort( pBegin, pEnd );
		d = ( iSize - 1 ) * dPercentile;
		dFrac = d - (size_t)d;
		iOne = (size_t)d;
		iTwo = (size_t)( d + 1 );

		return ( ( iTwo >= iSize ) ? pBegin[ iOne ] :
			( ( pBegin[ iOne ] * ( 1 - dPercentile ) ) + ( pBegin[ iTwo ] * dPercentile ) ) ); }

	/*!
	 * \brief
	 * Returns the median value of the given vector.
	 * 
	 * \param vecData
	 * Vector from which median value is returned.
	 * 
	 * \returns
	 * Median value of the given vector.
	 * 
	 * \remarks
	 * Note that the given vector is modified (sorted) by this function.
	 * 
	 * \see
	 * Percentile
	 */
	template<class tType>
	static double Median( std::vector<tType>& vecData ) {

		return Percentile( vecData.begin( ), vecData.end( ), 0.5 ); }

	/*!
	 * \brief
	 * Winsorizes the requested number of items at each end of the given vector.
	 * 
	 * \param vecValues
	 * Vector to be sorted and Winsorized.
	 * 
	 * \param iCount
	 * Number of items at each end of the vector to be Winsorized.
	 * 
	 * \returns
	 * True if Winsorization was necessary, false if there were too few values.
	 * 
	 * Winsorization is the process of replacing the n largest and smallest elements of a sequence with
	 * copies of the n-1st largest and n-1st smallest, respectively.  For example, Winsorizing the sequence
	 * [1, 2, 3, 4, 5, 6, 7, 8] by two would return [3, 3, 3, 4, 5, 6, 6, 6].  This is a standard method for
	 * finding a robust average in the presence of outliers.
	 * 
	 * \remarks
	 * Note that the vector is modified (sorted) by this function.
	 */
	template<class tType>
	static bool Winsorize( std::vector<tType>& vecValues, size_t iCount = 1 ) {

		return Winsorize( vecValues.begin( ), vecValues.end( ), iCount ); }

	/*!
	 * \brief
	 * Winsorizes the requested number of items at each end of the given sequence.
	 * 
	 * \param pBegin
	 * Pointer/iterator to the beginning of the sequence.
	 * 
	 * \param pEnd
	 * Pointer/iterator to the end of the sequence.
	 * 
	 * \param iCount
	 * Number of items at each end of the sequence to be Winsorized.
	 * 
	 * \returns
	 * True if Winsorization was necessary, false if there were too few values.
	 * 
	 * Winsorization is the process of replacing the n largest and smallest elements of a sequence with
	 * copies of the n-1st largest and n-1st smallest, respectively.  For example, Winsorizing the sequence
	 * [1, 2, 3, 4, 5, 6, 7, 8] by two would return [3, 3, 3, 4, 5, 6, 6, 6].  This is a standard method for
	 * finding a robust average in the presence of outliers.
	 * 
	 * \remarks
	 * Note that the sequence is modified (sorted) by this function.
	 */
	template<class tType>
	static bool Winsorize( tType pBegin, tType pEnd, size_t iCount = 1 ) {
		size_t	i, iLength;

		iLength = pEnd - pBegin;
		if( iLength < ( ( 2 * iCount ) + 1 ) )
			return false;
		std::sort( pBegin, pEnd );
		for( i = 0; i < iCount; ++i ) {
			pBegin[ i ] = pBegin[ i + iCount ];
			pEnd[ -1 - i ] = pEnd[ -1 - iCount ]; }

		return true; }

	// P-value tests
	static double LjungBox( const float* adX, size_t iN, size_t iH );

	/*!
	 * \brief
	 * Calculate the p-value of a Ljung-Box portmanteau test for autocorrelation randomness.
	 * 
	 * \param adX
	 * Array of values to test.
	 * 
	 * \param iN
	 * Number of values in array.
	 * 
	 * \returns
	 * Chi-squared of Q = iN * (iN + 2) * sum(autocorrelation(adX, i)^2 / (iN - i), i = 1..(iN - 1)) with
	 * iN - 1 degrees of freedom.
	 */
	static double LjungBox( const float* adX, size_t iN ) {

		return LjungBox( adX, iN, iN - 1 ); }

	/*!
	 * \brief
	 * Return the two-tailed p-value of a lognormal distribution.
	 * 
	 * \param dX
	 * Sample point.
	 * 
	 * \param dMean
	 * Mean of lognormal.
	 * 
	 * \param dStdev
	 * Standard deviation of lognormal.
	 * 
	 * \returns
	 * For dCDF = CStatistics::LognormalCDF (dX, dMean, dVariance), 2 * ((dX > exp(dMean)) ?
	 * (1 - dCDF) : dCDF).
	 */
	static double PValueLognormal( double dX, double dMean, double dStdev ) {
		double	dCDF;

		dCDF = LognormalCDF( dX, dMean, dStdev );
		if( dX > exp( dMean ) )
			dCDF = 1 - dCDF;

		return ( 2 * dCDF ); }

	/*!
	 * \brief
	 * Return the two-tailed p-value of a Pearson correlation.
	 * 
	 * \param dR
	 * Pearson correlation.
	 * 
	 * \param iN
	 * Length of correlated vectors.
	 * 
	 * \returns
	 * P-value corresponding to the given correlation and array size.
	 * 
	 * \see
	 * CMeasurePearson
	 */
	static double PValuePearson( double dR, size_t iN ) {
		double	dT, dF;

		dF = iN - 2;
		dT = dR * sqrt( dF / ( 1 - ( dR * dR ) ) );
		return IncompleteBeta( dF / 2, 0.5, dF / ( dF + ( dT * dT ) ) ); }

	static double PValueKolmogorovSmirnov( double dD ) {
		static const float	c_dEpsilon1	= 0.001f;
		static const float	c_dEpsilon2	= 1e-8f;
		static const float	c_dEpsilon3	= 0.475f;
		double	d, dRet, dCur, dPrev;
		size_t	i, iIterations;

		if( !dD )
			return 1;

		dD = -2 * pow( dD, 2 );
		iIterations = max( (size_t)250, (size_t)( 1 / dD ) );
		for( dRet = dPrev = 0,i = 1; i < iIterations; ++i ) {
			dCur = exp( i * i * dD );
			if( !( i % 2 ) )
				dCur *= -1;
			dRet += dCur;
			d = fabs( dCur );
			if( ( ( ( d / dRet ) < c_dEpsilon1 ) && ( dRet > c_dEpsilon3 ) ) ||
				( ( d / dPrev ) < c_dEpsilon1 ) || ( ( d / dRet ) < c_dEpsilon2 ) )
				break;
			dPrev = d; }
		if( dRet != 1 )
			dRet *= 2;

		return dRet; }

	/*!
	 * \brief
	 * Return the p-value of a t-test between the two given array statistics assuming equal variance.
	 * 
	 * \param dMeanOne
	 * Mean of the first sample.
	 * 
	 * \param dVarianceOne
	 * Variance of the first sample.
	 * 
	 * \param iNOne
	 * Number of elements in the first sample.
	 * 
	 * \param dMeanTwo
	 * Mean of the second sample.
	 * 
	 * \param dVarianceTwo
	 * Variance of the second sample.
	 * 
	 * \param iNTwo
	 * Number of elements in the second sample.
	 * 
	 * \returns
	 * P-value of T = (dMeanOne - dMeanTwo) / sqrt(((((iNOne - 1) * dVarianceOne) + ((iNTwo - 1) *
	 * dVarianceTwo)) / (iNOne + iNTwo - 2)) * ((1 / iNOne) + (1 / iNTwo)))
	 */
	static double TTestStudent( double dMeanOne, double dVarianceOne, size_t iNOne, double dMeanTwo,
		double dVarianceTwo, size_t iNTwo ) {
		size_t	iDegFree;
		double	dPoolVar, dT;

		iDegFree = iNOne + iNTwo - 2;
		dPoolVar = ( ( ( iNOne - 1 ) * dVarianceOne ) + ( ( iNTwo - 1 ) * dVarianceTwo ) ) / iDegFree;
		dT = ( dMeanOne - dMeanTwo ) / sqrt( dPoolVar * ( ( 1.0 / iNOne ) + ( 1.0 / iNTwo ) ) );

		return IncompleteBeta( 0.5 * iDegFree, 0.5, iDegFree / ( iDegFree + ( dT * dT ) ) ); }

	/*!
	 * \brief
	 * Return the p-value of a t-test between the two given array statistics without assuming equal variance.
	 * 
	 * \param dMeanOne
	 * Mean of the first sample.
	 * 
	 * \param dVarianceOne
	 * Variance of the first sample.
	 * 
	 * \param iNOne
	 * Number of elements in the first sample.
	 * 
	 * \param dMeanTwo
	 * Mean of the second sample.
	 * 
	 * \param dVarianceTwo
	 * Variance of the second sample.
	 * 
	 * \param iNTwo
	 * Number of elements in the second sample.
	 * 
	 * \returns
	 * P-value of T = (dMeanOne - dMeanTwo) / sqrt(((((iNOne - 1) * dVarianceOne) + ((iNTwo - 1) *
	 * dVarianceTwo)) / (iNOne + iNTwo - 2)) * ((1 / iNOne) + (1 / iNTwo)))
	 */
	static double TTestWelch( double dMeanOne, double dVarianceOne, size_t iNOne, double dMeanTwo,
		double dVarianceTwo, size_t iNTwo ) {
		double	dDegFree, dT;

		dDegFree = ( dVarianceOne / iNOne ) + ( dVarianceTwo / iNTwo );
		dDegFree = ( dDegFree * dDegFree ) / ( ( ( dVarianceOne * dVarianceOne ) / iNOne / iNOne / ( iNOne - 1 ) ) +
			( ( dVarianceTwo * dVarianceTwo ) / iNTwo / iNTwo / ( iNTwo - 1 ) ) );
		dT = ( dMeanOne - dMeanTwo ) / sqrt( ( dVarianceOne / iNOne ) + ( dVarianceTwo / iNTwo ) );

		return IncompleteBeta( 0.5 * dDegFree, 0.5, dDegFree / ( dDegFree + ( dT * dT ) ) ); }

	template<class tType>
	static double KSTest( const tType* aOne, const tType* aTwo, size_t iN ) {
		size_t	i;
		tType	SumOne, SumTwo, CumOne, CumTwo;
		float	d, dMax;

		if( !iN )
			return 1;
		for( SumOne = SumTwo = i = 0; i < iN; ++i ) {
			SumOne += aOne[ i ];
			SumTwo += aTwo[ i ]; }
		for( dMax = 0,CumOne = CumTwo = i = 0; i < iN; ++i ) {
			CumOne += aOne[ i ];
			CumTwo += aTwo[ i ];
			if( ( d = fabs( ( (float)CumOne / SumOne ) - ( (float)CumTwo / SumTwo ) ) ) > dMax )
				dMax = d; }

		d = sqrt( (float)( SumOne * SumTwo ) / ( SumOne + SumTwo ) );
		return PValueKolmogorovSmirnov( dMax * ( d + 0.12 + ( 0.11 / d ) ) ); }

	// Evaluation statistics
	static double WilcoxonRankSum( const CDat& DatData, const CDat& DatAnswers,
		const std::vector<bool>& vecfGenesOfInterest, bool fInvert = false );

	// Probability distributions
	static double HypergeometricCDF( size_t iHitsOne, size_t iSizeOne, size_t iHitsTwo, size_t iSizeTwo );
	static double SampleGammaStandard( double dShape );
	static double SampleGammaLogStandard( double dXX );
	static double SampleNormalStandard( );
	static double SampleExponentialStandard( );

	/*!
	 * \brief
	 * Return the binomial p-value of obtaining at least the given sample with the given number of
	 * observations and probability of success.
	 * 
	 * \param iObservations
	 * Count parameter of the binomial.
	 * 
	 * \param iSample
	 * Point at which to sample the CDF.
	 * 
	 * \param dProbability
	 * Probability parameter of the binomial.
	 * 
	 * \returns
	 * 1 - CStatistics::NormalCDF ((iObservations - (iSample * dProbability)) / sqrt(iSample *
	 * dProbability * (1 - dProbability), 0, 1)
	 * 
	 * \remarks
	 * Implementation courtesy of Press WH, Teukolsky SA, Vetterling WT, Flannery BP.  Numerical Recipes in C,
	 * 1992, Cambridge University Press.
	 */
	static double BinomialCDF( size_t iObservations, size_t iSample, double dProbability ) {
		double	d;

		d = ( iObservations - ( iSample * dProbability ) ) / sqrt( iSample * dProbability *
			( 1 - dProbability ) );
		return ( 1 - NormalCDF( d, 0, 1 ) ); }

	/*!
	 * \brief
	 * Calculate the hypergeometric probability distribution given the sizes and overlap of two sets.
	 * 
	 * \param iHitsOne
	 * Number of hits in the first (query) set.
	 * 
	 * \param iSizeOne
	 * Size of the first (query) set.
	 * 
	 * \param iHitsTwo
	 * Number of hits in the second (background) set.
	 * 
	 * \param iSizeTwo
	 * Size of the second (background) set.
	 * 
	 * \returns
	 * choose(iHitsTwo, iHitsOne) * choose(iSizeTwo - iHitsTwo, iSizeOne - iHitsOne) /
	 * choose(iSizeTwo, iSizeOne)
	 * 
	 * \remarks
	 * Calculated using the exponential of CStatistics::LogFact results for increased speed and precision.
	 */
	static double HypergeometricPDF( size_t iHitsOne, size_t iSizeOne, size_t iHitsTwo,
		size_t iSizeTwo ) {

		return exp( LogFact( iSizeTwo - iHitsTwo ) + LogFact( iHitsTwo ) + LogFact( iSizeOne ) +
			LogFact( iSizeTwo - iSizeOne ) - LogFact( iHitsOne ) -
			LogFact( iHitsTwo - iHitsOne ) - LogFact( iSizeTwo - iHitsTwo + iHitsOne - iSizeOne ) -
			LogFact( iSizeOne - iHitsOne ) - LogFact( iSizeTwo ) ); }

	/*!
	 * \brief
	 * Calculate a p-value for the given T and degrees of freedom.
	 * 
	 * \param dT
	 * T value at which to sample the t-distribution.
	 * 
	 * \param iDF
	 * Degrees of freedom of the desired t-distribution.
	 * 
	 * \returns
	 * p-value of the given T and degrees of freedom.
	 */
	static double TCDF( double dT, size_t iDF ) {

		return IncompleteBeta( 0.5 * iDF, 0.5, iDF / ( iDF + ( dT * dT ) ) ); }

	/*!
	 * \brief
	 * CDF of a lognormal distribution with the given parameters at the given point.
	 * 
	 * \param dX
	 * Sample point.
	 * 
	 * \param dMean
	 * Mean of lognormal.
	 * 
	 * \param dStdev
	 * Standard deviation of lognormal.
	 * 
	 * \returns
	 * CStatistics::NormalCDF (log(dX), dMean, dVariance)
	 */
	static double LognormalCDF( double dX, double dMean, double dStdev ) {

		return ( ( dX > 0 ) ? NormalCDF( log( dX ), dMean, dStdev ) : 0 ); }

	/*!
	 * \brief
	 * CDF of an inverse Gaussian distribution with the given parameters at the given point.
	 * 
	 * \param dX
	 * Sample point.
	 * 
	 * \param dMean
	 * Mean of inverse Gaussian.
	 * 
	 * \param dLambda
	 * Lambda parameter of inverse Gaussian.
	 * 
	 * \returns
	 * CStatistics::Normal01CDF ( sqrt( dLambda / dX ) * ( ( dX / dMean ) - 1 ) ) +
	 * 		exp( ( 2 * dLambda / dMean  ) + log ( Normal01CDF( -sqrt( dLambda / dX ) *
	 * 		( ( dX / dMean ) + 1 ) ) )
	 */
	static double InverseGaussianCDF( double dX, double dMean, double dLambda ) {

		return ( Normal01CDF( sqrt( dLambda / dX ) * ( ( dX / dMean ) - 1 ) ) +
			exp( ( 2 * dLambda / dMean  ) + log ( Normal01CDF( -sqrt( dLambda / dX ) *
			( ( dX / dMean ) + 1 ) ) ) ) ); }

	/*!
	 * \brief
	 * CDF of a normal distribution with the given parameters at the given point.
	 * 
	 * \param dX
	 * Sample point.
	 * 
	 * \param dMean
	 * Mean of normal.
	 * 
	 * \param dStdev
	 * Standard deviation of normal.
	 * 
	 * \returns
	 * NCDF((dX - dMean) / dStdev), for NCDF a normal CDF with mean 0, variance 1.
	 */
	static double NormalCDF( double dX, double dMean, double dStdev ) {

		return Normal01CDF( ( dX - dMean ) / dStdev ); }

	/*!
	 * \brief
	 * CDF of a standard normal distribution.
	 * 
	 * \param dX
	 * Sample point.
	 * 
	 * \returns
	 * NCDF(dX), for NCDF a normal CDF with mean 0, variance 1.
	 */
	static double Normal01CDF( double dX ) {

		return CStatisticsImpl::Normal01CDF( dX ); }

	/*!
	 * \brief
	 * Return a random sample from a chi-squared distribution with the given degrees of freedom.
	 * 
	 * \param iDF
	 * Degrees of freedom of the chi-squared distribution.
	 * 
	 * \returns
	 * Random sample from a chi-squared distribution with the requested degrees of freedom.
	 */
	static double SampleChi2( size_t iDF ) {

		return ( 2 * SampleGamma( 1, iDF / 2 ) ); }

	/*!
	 * \brief
	 * Return a random sample from a gamma distribution with the given location and shape parameters.
	 * 
	 * \param dLocation
	 * Location parameter of gamma function to sample.
	 * 
	 * \param dShape
	 * Shape parameter of gamma function to sample.
	 * 
	 * \returns
	 * Random sample from a gamma function with the given shape and location parameters.
	 */
	static double SampleGamma( double dLocation, double dShape ) {

		return ( SampleGammaStandard( dShape ) / dLocation ); }

	/*!
	 * \brief
	 * Calculate a beta probability density at the given point with the given minimum, maximu, alpha, and
	 * beta parameters.
	 * 
	 * \param dX
	 * Point at which to sample the beta probability distribution.
	 * 
	 * \param dMinimum
	 * Minimum parameter of the beta distribution.
	 * 
	 * \param dMaximum
	 * Maximum parameter of the beta distribution.
	 * 
	 * \param dAlpha
	 * Alpha parameter of the beta distribution.
	 * 
	 * \param dBeta
	 * Beta parameter of the beta distribution.
	 * 
	 * \returns
	 * Beta probability density at the given point for a distribution with the requested parameters.
	 * 
	 * \remarks
	 * Implementation courtesy of Press WH, Teukolsky SA, Vetterling WT, Flannery BP.  Numerical Recipes in C,
	 * 1992, Cambridge University Press.
	 */
	static double BetaPDF( double dX, double dMinimum, double dMaximum, double dAlpha, double dBeta ) {
		double	dFunc, dLocation, dScale;

		dLocation = dMinimum;
		dScale = dMaximum - dMinimum;
		dX = ( dX - dLocation ) / dScale;
		dFunc = exp( SampleGammaLogStandard( dAlpha ) + SampleGammaLogStandard( dBeta ) -
			SampleGammaLogStandard( dAlpha + dBeta ) );
		return ( pow( dX, dAlpha - 1 ) * pow( 1 - dX, dBeta - 1 ) / dFunc / dScale ); }

	/*!
	 * \brief
	 * Calculate a normal probability density at the given point for the given mean and variance.
	 * 
	 * \param dX
	 * Point at which to calculate the normal distribution.
	 * 
	 * \param dMu
	 * Mean of the normal distribution.
	 * 
	 * \param dSigma
	 * Variance of the normal distribution.
	 * 
	 * \returns
	 * exp(-(dX - dMu)^2 / (2 * dSigma^2)) / (dSigma * sqrt(2*PI))
	 */
	static double NormalPDF( double dX, double dMu, double dSigma ) {
		static const double	c_dS2P	= sqrt( 2 * 3.1415926535898 );
		double	d;

		d = dX - dMu;
		return ( exp( -( d * d ) / ( 2 * dSigma * dSigma ) ) / ( dSigma * c_dS2P ) ); }
};

}

#endif // STATISTICS_H
