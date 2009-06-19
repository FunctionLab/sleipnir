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

#include <limits>

#include "statisticsi.h"

#include "fullmatrix.h"

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
	 * sum((x - dMean)^2) / n
	 * 
	 * \see
	 * Average
	 */
	template<class tType>
	static double Variance( const tType Begin, const tType End, double dMean ) {
		double	dRet;
		size_t	iN;

		Sums( Begin, End, NULL, &dRet, &iN );
		return ( iN ? max( ( dRet / iN ) - ( dMean * dMean ), 0.0 ) : 0 ); }

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
	 * \returns
	 * sum((x - dMean)^2) / (n - 1)
	 * 
	 * \see
	 * Average
	 */
	template<class tType>
	static double Variance( const tType Begin, const tType End ) {
		double	dSum, dRet;
		size_t	iN;

		Sums( Begin, End, &dSum, &dRet, &iN );
		if( !iN )
			return 0;
		dSum /= iN;
		return max( ( dRet / iN ) - ( dSum * dSum ), 0.0 ); }

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

		return Variance( vecValues.begin( ), vecValues.end( ) ); }

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
		double	dRet;
		size_t	iN;

		Sums( Begin, End, &dRet, NULL, &iN );
		return ( iN ? ( dRet / iN ) : dRet ); }

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
		std::sort( pBegin, pEnd );
		while( iSize && CMeta::IsNaN( pBegin[ iSize - 1 ] ) )
			--iSize;
		if( !iSize )
			return CMeta::GetNaN( );
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

	/*!
	 * \brief
	 * Calculates Cohen's D, a modified z-score effect size measure based on unpooled variance.
	 * 
	 * \param vecOne
	 * First array to be compared.
	 * 
	 * \param vecTwo
	 * Second array to be compared.
	 * 
	 * \param pdAverage
	 * If non-null, output average value of elements in both arrays.
	 * 
	 * \returns
	 * Cohen's D effect size measurement of the difference between input arrays.
	 * 
	 * Calculate's Cohen's D, a z-score-like effect size measurement based on unpooled variance.
	 * For two input arrays A1 and A2 with means m(Ai) and standard deviations s(Ai), the standard
	 * z-score is (m(A1)-m(A2))/s(A1 u A2).  Cohen's D uses the unpooled variance 2(m(A1)-m(A2))/(s(A1)+s(A2)).
	 * 
	 * \see
	 * ZScore
	 */
	template<class tType>
	static double CohensD( const std::vector<tType>& vecOne, const std::vector<tType>& vecTwo,
		double* pdAverage = NULL ) {

		return CohensD( vecOne.begin( ), vecOne.end( ), vecTwo.begin( ), vecTwo.end( ), pdAverage ); }

	/*!
	 * \brief
	 * Calculates Cohen's D, a modified z-score effect size measure based on unpooled variance.
	 * 
	 * \param BeginOne
	 * First element of the first array to be compared.
	 * 
	 * \param EndOne
	 * End of last element of the first array to be compared.
	 * 
	 * \param BeginTwo
	 * First element of the second array to be compared.
	 * 
	 * \param EndTwo
	 * End of last element of the second array to be compared.
	 * 
	 * \param pdAverage
	 * If non-null, output average value of elements in both arrays.
	 * 
	 * \returns
	 * Cohen's D effect size measurement of the difference between input arrays.
	 * 
	 * Calculate's Cohen's D, a z-score-like effect size measurement based on unpooled variance.
	 * For two input arrays A1 and A2 with means m(Ai) and standard deviations s(Ai), the standard
	 * z-score is (m(A1)-m(A2))/s(A1 u A2).  Cohen's D uses the unpooled variance 2(m(A1)-m(A2))/(s(A1)+s(A2)).
	 * 
	 * \see
	 * ZScore
	 */
	template<class tType>
	static double CohensD( const tType BeginOne, const tType EndOne, const tType BeginTwo, const tType EndTwo,
		double* pdAverage = NULL ) {
		double	dAveOne, dAveTwo, dVarOne, dVarTwo, dStd;
		size_t	iOne, iTwo;

		Sums( BeginOne, EndOne, &dAveOne, &dVarOne, &iOne );
		if( iOne ) {
			dAveOne /= iOne;
			dVarOne = ( dVarOne / iOne ) - ( dAveOne * dAveOne ); }
		if( pdAverage )
			*pdAverage = dAveOne;
		Sums( BeginTwo, EndTwo, &dAveTwo, &dVarTwo, &iTwo );
		if( iTwo ) {
			dAveTwo /= iTwo;
			dVarTwo = ( dVarTwo / iTwo ) - ( dAveTwo * dAveTwo ); }

		dStd = sqrt( ( dVarOne + dVarTwo ) / 2 );
		return ( dStd ? ( ( dAveOne - dAveTwo ) / dStd ) :
			( ( dAveOne == dAveTwo ) ? 0 : DBL_MAX ) ); }

	/*!
	 * \brief
	 * Calculates the z-score effect size measure.
	 * 
	 * \param vecOne
	 * First array to be compared.
	 * 
	 * \param vecTwo
	 * Second array to be compared.
	 * 
	 * \returns
	 * Z-score effect size measurement of the difference between input arrays.
	 * 
	 * For two input arrays A1 and A2 with means m(Ai) and standard deviations s(Ai), the standard
	 * z-score is (m(A1)-m(A2))/s(A1 u A2).
	 * 
	 * \see
	 * CohensD
	 */
	template<class tType>
	static double ZScore( const std::vector<tType>& vecOne, const std::vector<tType>& vecTwo ) {

		return ZScore( vecOne.begin( ), vecOne.end( ), vecTwo.begin( ), vecTwo.end( ) ); }

	/*!
	 * \brief
	 * Calculates the z-score effect size measure.
	 * 
	 * \param BeginOne
	 * First element of the first array to be compared.
	 * 
	 * \param EndOne
	 * End of last element of the first array to be compared.
	 * 
	 * \param BeginTwo
	 * First element of the second array to be compared.
	 * 
	 * \param EndTwo
	 * End of last element of the second array to be compared.
	 * 
	 * \param pdAverage
	 * If non-null, output average value of elements in both arrays.
	 * 
	 * \returns
	 * Z-score effect size measurement of the difference between input arrays.
	 * 
	 * For two input arrays A1 and A2 with means m(Ai) and standard deviations s(Ai), the standard
	 * z-score is (m(A1)-m(A2))/s(A1 u A2).
	 * 
	 * \see
	 * CohensD
	 */
	template<class tType>
	static double ZScore( const tType BeginOne, const tType EndOne, const tType BeginTwo, const tType EndTwo,
		double* pdAverage = NULL ) {
		double	dAveOne, dAveTwo, dVarOne, dVarTwo, dAve, dStd;
		size_t	iOne, iTwo;

		Sums( BeginOne, EndOne, &dAveOne, &dVarOne, &iOne );
		Sums( BeginTwo, EndTwo, &dAveTwo, &dVarTwo, &iTwo );
		dAve = ( iOne || iTwo ) ? ( ( dAveOne + dAveTwo ) / ( iOne + iTwo ) ) : 0;
		if( iOne )
			dAveOne /= iOne;
		if( pdAverage )
			*pdAverage = dAveOne;
		dStd = ( iOne || iTwo ) ? sqrt( ( ( dVarOne + dVarTwo ) / ( iOne + iTwo ) ) - ( dAve * dAve ) ) : 0;

		return ( dStd ? ( ( dAveOne - dAve ) / dStd ) :
			( ( dAveOne == dAve ) ? 0 : DBL_MAX ) ); }

	/*!
	 * \brief
	 * Returns the z-test p-value of the given z-score and element count.
	 * 
	 * \param dZScore
	 * Z-score to be converted into a p-value.
	 * 
	 * \param iN
	 * Number of elements used to generate the given z-score.
	 * 
	 * \returns
	 * P-value equivalent to the given z-score and element count.
	 * 
	 * \remarks
	 * Equivalent to 1 - normcdf( 0, sqrt( iN ), |dZScore| )
	 * 
	 * \see
	 * ZScore
	 */
	static double ZTest( double dZScore, size_t iN ) {

		return ( 1 - Normal01CDF( fabs( dZScore ) * sqrt( (double)iN ) ) ); }

	/*!
	 * \brief
	 * Returns the root-mean-square error distance between two input arrays.
	 * 
	 * \param BeginOne
	 * First element of the first array to be compared.
	 * 
	 * \param EndOne
	 * End of last element of the first array to be compared.
	 * 
	 * \param BeginTwo
	 * First element of the second array to be compared.
	 * 
	 * \param EndTwo
	 * End of last element of the second array to be compared.
	 * 
	 * \returns
	 * Root-mean-square error between the two input arrays.
	 * 
	 * Calculates the root-mean-square error (RMSE) between two arrays, equal to
	 * sqrt( sum( (A1[i] - A2[i])^2 ) / min(|A1|, |A2|) ).
	 * 
	 * \remarks
	 * If the two input arrays are of different lengths m and n, only min(m, n) elements of either
	 * array are compared.
	 */
	template<class tType>
	static double RootMeanSquareError( tType BeginOne, tType EndOne, tType BeginTwo, tType EndTwo ) {
		tType	CurOne, CurTwo;
		double	d, dRet;
		size_t	iN;

		for( dRet = 0,CurOne = BeginOne,CurTwo = BeginTwo; ( CurOne != EndOne ) && ( CurTwo != EndTwo );
			++CurOne,++CurTwo ) {
			d = *CurOne - *CurTwo;
			dRet += d * d; }
		iN = min( EndOne - BeginOne, EndTwo - BeginTwo );

		return ( iN ? sqrt( dRet / iN ) : 0 ); }

	/*!
	 * \brief
	 * Returns the Jensen-Shannon divergence between two discrete probability distributions
	 * represented as arrays.
	 * 
	 * \param BeginOne
	 * First element of the first array to be compared.
	 * 
	 * \param EndOne
	 * End of last element of the first array to be compared.
	 * 
	 * \param BeginTwo
	 * First element of the second array to be compared.
	 * 
	 * \param EndTwo
	 * End of last element of the second array to be compared.
	 * 
	 * \returns
	 * Jensen-Shannon divergence between the given probability distributions in bits.
	 * 
	 * The Jensen-Shannon or JS-divergence between two discrete probability distributions A1 and A2
	 * is a symmetrized Kullback-Leibler divergence defined as ( KL(A1, A2) + KL(A2, A1 ) / 2.
	 * 
	 * \remarks
	 * The elements of each array should sum to one, and the two arrays should be of the same length.
	 * 
	 * \see
	 * KullbackLeiblerDivergence
	 */
	template<class tType>
	static double JensenShannonDivergence( tType BeginOne, tType EndOne, tType BeginTwo, tType EndTwo ) {

		return ( ( KullbackLeiblerDivergence( BeginOne, EndOne, BeginTwo, EndTwo ) +
			KullbackLeiblerDivergence( BeginTwo, EndTwo, BeginOne, EndOne ) ) / 2 ); }

	/*!
	 * \brief
	 * Returns the Kullback-Leibler divergence between two discrete probability distributions
	 * represented as arrays.
	 * 
	 * \param BeginOne
	 * First element of the first array to be compared.
	 * 
	 * \param EndOne
	 * End of last element of the first array to be compared.
	 * 
	 * \param BeginTwo
	 * First element of the second array to be compared.
	 * 
	 * \param EndTwo
	 * End of last element of the second array to be compared.
	 * 
	 * \returns
	 * Kullback-Leibler divergence between the given probability distributions in bits.
	 * 
	 * The Kullback-Leibler or KL-divergence between two discrete probability distributions A1 and A2
	 * is defined as sum( A1[i] * log( A1[i]/A2[i] ) ).
	 * 
	 * \remarks
	 * The elements of each array should sum to one, and the two arrays should be of the same length.
	 * 
	 * \see
	 * JensenShannonDivergence
	 */
	template<class tType>
	static double KullbackLeiblerDivergence( tType BeginOne, tType EndOne, tType BeginTwo, tType EndTwo ) {
		double	dRet;
		tType	CurOne, CurTwo;

		if( ( EndOne - BeginOne ) != ( EndTwo - BeginTwo ) )
			return CMeta::GetNaN( );

		for( dRet = 0,CurOne = BeginOne,CurTwo = BeginTwo; ( CurOne != EndOne ) && ( CurTwo != EndTwo );
			++CurOne,++CurTwo )
			dRet += *CurOne * log( *CurOne / *CurTwo );

		return ( dRet / log( 2.0 ) ); }

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
		static const double	c_dEpsilon	= 1e-10;
		double	dT, dF;

		if( iN < 2 )
			return 1;
		if( ( 1 - dR ) < c_dEpsilon )
			return 0;
		dF = iN - 2;
		dT = dR * sqrt( dF / ( 1 - ( dR * dR ) ) );
		return IncompleteBeta( dF / 2, 0.5, dF / ( dF + ( dT * dT ) ) ); }

	/*!
	 * \brief
	 * Returns the p-value corresponding to a D-score obtained from a Kolmogorov-Smirnov test.
	 * 
	 * \param dD
	 * D-value obtained from a Kolmogorov-Smirnov test.
	 * 
	 * \param iM
	 * Number of elements in first tested array.
	 * 
	 * \param iN
	 * Number of elements in second tested array.
	 * 
	 * \returns
	 * P-value equivalent to the given D-score and element count.
	 * 
	 * \see
	 * CMeasureKolmogorovSmirnov
	 */
	static double PValueKolmogorovSmirnov( double dD, size_t iM, size_t iN ) {
		static const float	c_dEpsilon1	= 0.001f;
		static const float	c_dEpsilon2	= 1e-8f;
		static const float	c_dEpsilon3	= 0.475f;
		double	d, dRet, dCur, dPrev;
		size_t	i, iIterations;

		if( !dD )
			return 1;

		d = sqrt( (double)( iM * iN ) / ( iM + iN ) );
		iIterations = max( (size_t)250, (size_t)( 1 / dD ) );
		dD = -2 * pow( dD * d, 2 );
// This is from NR, but it disagrees with R's results
//		dD = -2 * pow( dD * ( d + 0.12 + ( 0.11 / d ) ), 2 );
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
		if( dRet < 1 )
			dRet = min( 2 * dRet, 1.0 );

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

	/*!
	 * \brief
	 * Return the p-value of an f-test between the two given array statistics to determine equality of variance.
	 * 
	 * \param dVarianceOne
	 * Variance of the first sample.
	 * 
	 * \param iNOne
	 * Number of elements in the first sample.
	 * 
	 * \param dVarianceTwo
	 * Variance of the second sample.
	 * 
	 * \param iNTwo
	 * Number of elements in the second sample.
	 * 
	 * \returns
	 * P-value of F = dVarianceOne / dVarianceTwo.
	 */
	static double FTest( double dVarianceOne, size_t iNOne, double dVarianceTwo, size_t iNTwo ) {
		double	dRet, dF;
		size_t	iDF1, iDF2;

		if( dVarianceOne < dVarianceTwo ) {
			std::swap( dVarianceOne, dVarianceTwo );
			std::swap( iNOne, iNTwo ); }
		dF = dVarianceOne / dVarianceTwo;
		iDF1 = iNOne - 1;
		iDF2 = iNTwo - 1;

		dRet = 2 * IncompleteBeta( 0.5 * iDF2, 0.5 * iDF1, iDF2 / ( iDF2 + ( iDF1 * dF ) ) );
		if( dRet > 1 )
			dRet = 2 - dRet;

		return dRet; }

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
	 * Calculate the Skellam probability distribution given a point and the two distribution parameters.
	 * 
	 * \param iX
	 * Value at which PDF should be evaluated.
	 * 
	 * \param dMu1
	 * First distribution parameter (first expected value).
	 * 
	 * \param dMu2
	 * Second distribution parameter (second expected value).
	 * 
	 * \returns
	 * Value of a Skellam distribution with the requested parameters at the given point.
	 */
	static double SkellamPDF( size_t iX, double dMu1, double dMu2 ) {

		return ( exp( -( dMu1 + dMu2 ) ) * pow( dMu1 / dMu2, 0.5 * iX ) *
			ModifiedBesselI( iX, 2 * sqrt( dMu1 * dMu2 ) ) ); }

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

		return ( 1 - IncompleteBeta( 0.5 * iDF, 0.5, iDF / ( iDF + ( dT * dT ) ) ) ); }

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
	 * NCDF((dX - dMean) / dStdev), for NCDF a normal CDF with mean 0, standard deviation 1.
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
	 * Inverse CDF of a standard normal distribution.
	 * 
	 * \param dX
	 * Sample point.
	 * 
	 * \returns
	 * Value of dY, for dX = NCDF(dY) and NCDF a normal CDF with mean 0, variance 1.
	 */
	static double InverseNormal01CDF( double dX );

	/*!
	 * \brief
	 * CDF of a multivariate normal distribution with the given parameters at the given point.
	 * 
	 * \param vecdX
	 * Sample point.
	 * 
	 * \param vecdMu
	 * Mean of normal.
	 * 
	 * \param MatSigmaCholesky
	 * Cholesky decomposition of covariance matrix of normal.
	 * 
	 * \param iN
	 * Number of elements used to calculate sample point.
	 * 
	 * \param dMaxError
	 * Performance parameter; maximum error tolerance of return value.
	 * 
	 * \param dMaxCI
	 * Performance parameter; confidence interval of return value.
	 * 
	 * \param iMaxIterations
	 * Performance parameter; maximum iteratios to calculate return value.
	 * 
	 * \returns
	 * NCDF((vecdX - vecdMu) / MatSigma), for NCDF a normal CDF with mean 0, standard deviation 1.
	 * 
	 * \remarks
	 * Implementation courtesy of Numerical Recipes.
	 * 
	 * \see
	 * CholeskyDecomposition
	 */
	static double MultivariateNormalCDF( const std::vector<float>& vecdX, const std::vector<float>& vecdMu,
		const CDataMatrix& MatSigmaCholesky, size_t iN = 1, float dMaxError = 0.01, float dMaxCI = 0.99,
		size_t iMaxIterations = 300 ) {
		std::vector<double>	vecdDiff, vecdY, vecdF;
		size_t				i, j, iIterations;
		double				d, dRet, dVar, dError, dAlpha, dQ, dN;

		if( vecdX.empty( ) || ( vecdX.size( ) != vecdMu.size( ) ) ||
			( vecdX.size( ) != MatSigmaCholesky.GetRows( ) ) ||
			( vecdX.size( ) != MatSigmaCholesky.GetColumns( ) ) )
			return CMeta::GetNaN( );

		dN = sqrt( (double)iN );
		vecdDiff.resize( vecdX.size( ) );
		for( i = 0; i < vecdDiff.size( ); ++i )
			vecdDiff[ i ] = vecdX[ i ] - vecdMu[ i ];
		dAlpha = InverseNormal01CDF( dMaxCI );
		vecdF.resize( vecdX.size( ) );
		vecdF[ 0 ] = Normal01CDF( dN * vecdDiff[ 0 ] / MatSigmaCholesky.Get( 0, 0 ) );
		vecdY.resize( vecdX.size( ) );

		dRet = dVar = 0;
		dError = 2 * dMaxError;
		for( iIterations = 1; ( iIterations <= iMaxIterations ) && ( dError > dMaxError ); ++iIterations ) {
			for( i = 1; i < vecdY.size( ); ++i ) {
				vecdY[ i - 1 ] = InverseNormal01CDF( vecdF[ i - 1 ] * rand( ) / RAND_MAX );
				dQ = 0;
				for( j = 0; j < i; ++j )
					dQ += MatSigmaCholesky.Get( j, i ) * vecdY[ j ] / dN;
				vecdF[ i ] = Normal01CDF( dN * ( vecdDiff[ i ] - dQ ) / MatSigmaCholesky.Get( i, i ) ) *
					vecdF[ i - 1 ]; }
			d = ( vecdF.back( ) - dRet ) / iIterations;
			dRet += d;
			dVar = ( ( iIterations - 2 ) * dVar / iIterations ) + ( d * d );
			dError = dAlpha * sqrt( dVar ); }

		return dRet; }

	/*!
	 * \brief
	 * Calculates the value of a multivariate normal probability density function at the requested point.
	 * 
	 * \param vecdX
	 * Point at which multivariate normal distribution is evaluated.
	 * 
	 * \param vecdMu
	 * Mean of multivariate normal distribution.
	 * 
	 * \param MatSigma
	 * Covariance matrix of multivariate normal distribution.
	 * 
	 * \returns
	 * Value of multivariate normal distribution at the requested point.
	 * 
	 * \remarks
	 * vecdX and vecdMu must be of the same nonzero length, and MatSigma must be a square matrix with
	 * dimensions equal to that length.
	 * 
	 * \see
	 * MultivariateNormalCDF
	 */
	static double MultivariateNormalPDF( const std::vector<float>& vecdX, const std::vector<float>& vecdMu,
		const CDataMatrix& MatSigma ) {
		CDataMatrix		MatLU, MatInv;
		vector<size_t>	veciIndices;
		bool			fEven;
		double			dDet;

		if( !MatSigma.GetRows( ) || ( MatSigma.GetRows( ) != MatSigma.GetColumns( ) ) )
			return CMeta::GetNaN( );

		MatLU.Initialize( MatSigma.GetRows( ), MatSigma.GetColumns( ) );
		MatLU.Open( MatSigma );
		MatrixLUDecompose( MatLU, veciIndices, fEven );
		MatrixLUInvert( MatLU, veciIndices, MatInv );
		dDet = MatrixLUDeterminant( MatLU, fEven );

		return MultivariateNormalPDF( vecdX, vecdMu, sqrt( dDet ), MatInv ); }

	/*!
	 * \brief
	 * Calculates the value of a multivariate normal probability density function at the requested point.
	 * 
	 * \param vecdX
	 * Point at which multivariate normal distribution is evaluated.
	 * 
	 * \param vecdMu
	 * Mean of multivariate normal distribution.
	 * 
	 * \param dSigmaDetSqrt
	 * Square root of the determinant of the requested distribution's covariance matrix.
	 * 
	 * \param MatSigmaInv
	 * Inverse of the requested distribution's covariance matrix.
	 * 
	 * \returns
	 * Value of multivariate normal distribution at the requested point.
	 * 
	 * \remarks
	 * vecdX and vecdMu must be of the same nonzero length, and MatSigmaInv must be a square matrix with
	 * dimensions equal to that length.
	 * 
	 * \see
	 * MultivariateNormalCDF
	 */
	static double MultivariateNormalPDF( const std::vector<float>& vecdX, const std::vector<float>& vecdMu,
		double dSigmaDetSqrt, const CDataMatrix& MatSigmaInv ) {
		size_t			i;
		vector<float>	vecdXmMu, vecdXmMutS;
		double			d;

		if( !MatSigmaInv.GetRows( ) || ( MatSigmaInv.GetRows( ) != MatSigmaInv.GetColumns( ) ) ||
			vecdX.empty( ) || ( vecdX.size( ) != vecdMu.size( ) ) )
			return CMeta::GetNaN( );

		vecdXmMu.resize( vecdX.size( ) );
		for( i = 0; i < vecdXmMu.size( ); ++i )
			if( CMeta::IsNaN( vecdXmMu[ i ] = vecdX[ i ] - vecdMu[ i ] ) )
				vecdXmMu[ i ] = 0;
		MatrixMultiply( vecdXmMu, MatSigmaInv, vecdXmMutS );
		d = MatrixMultiply( vecdXmMutS, vecdXmMu );

		return ( ( 1 / pow( 2 * 3.1415926535898, vecdX.size( ) / 2.0 ) / dSigmaDetSqrt ) * exp( -d / 2 ) ); }

	/*!
	 * \brief
	 * Return a random sample from a chi-squared distribution with the given degrees of freedom.
	 * 
	 * \param iDF
	 * Degrees of freedom of the chi-squared distribution.
	 * 
	 * \returns
	 * Random sample from a chi-squared distribution with the requested degrees of freedom.
	 * 
	 * \see
	 * Chi2CDF
	 */
	static double SampleChi2( size_t iDF ) {

		return ( 2 * SampleGamma( 1, iDF / 2 ) ); }

	/*!
	 * \brief
	 * Return the value of a chi-squared cumulative density function at the requested point.
	 * 
	 * \param dC2
	 * Point at which chi-squared CDF should be evaluated.
	 * 
	 * \param iDF
	 * Degrees of freedom of the chi-squared distribution.
	 * 
	 * \returns
	 * Value of the chi-squared CDF at the given point with the requested degrees of freedom.
	 * 
	 * \see
	 * SampleChi2
	 */
	static double Chi2CDF( double dC2, size_t iDF ) {

		return CStatisticsImpl::Chi2CDF( sqrt( dC2 ), 0, 1, iDF ); }

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
	 * Calculate a normal probability density at the given point for the given mean and standard deviation.
	 * 
	 * \param dX
	 * Point at which to calculate the normal distribution.
	 * 
	 * \param dMu
	 * Mean of the normal distribution.
	 * 
	 * \param dSigma
	 * Standard deviation of the normal distribution.
	 * 
	 * \returns
	 * exp(-(dX - dMu)^2 / (2 * dSigma^2)) / (dSigma * sqrt(2*PI))
	 */
	static double NormalPDF( double dX, double dMu, double dSigma ) {
		static const double	c_dS2P	= sqrt( 2 * 3.1415926535898 );
		double	d;

		d = dX - dMu;
		return ( exp( -( d * d ) / ( 2 * dSigma * dSigma ) ) / ( dSigma * c_dS2P ) ); }

	/*!
	 * \brief
	 * Returns the value of an exponential probability density function at the requested point.
	 * 
	 * \param dX
	 * Point at which the exponential PDF is evaluated.
	 * 
	 * \param dLambda
	 * Rate parameter of the exponential PDF to be evaluated.
	 * 
	 * \returns
	 * Value of an exponeitial PDF at the requested point.
	 */
	static double ExponentialPDF( double dX, double dLambda ) {

		return ( dLambda * exp( -dLambda * dX ) ); }

	/*!
	 * \brief
	 * Calculates the Cholesky decomposition of the given matrix, overwriting it in the process.
	 * 
	 * \param Matrix
	 * Matrix to be decomposed.
	 * 
	 * \returns
	 * True if the decomposition was successful, false otherwise.
	 * 
	 * \remarks
	 * Matrix must be square.
	 */
	static bool CholeskyDecomposition( CDataMatrix& Matrix ) {
		std::vector<float>	vecdP;
		size_t				i, j, k;
		float				dSum;

		if( Matrix.GetRows( ) != Matrix.GetColumns( ) )
			return false;

		vecdP.resize( Matrix.GetRows( ) );
		for( i = 0; i < vecdP.size( ); ++i )
			for( j = i; j < vecdP.size( ); ++j ) {
				dSum = Matrix.Get( i, j );
				for( k = 0; k < i; ++k )
					dSum -= Matrix.Get( i, k ) * Matrix.Get( j, k );
				if( i == j )
					vecdP[ i ] = ( dSum <= 0 ) ? std::numeric_limits<float>::epsilon( ) : sqrt( dSum );
				else
					Matrix.Set( j, i, dSum / vecdP[ i ] ); }

		for( i = 0; i < vecdP.size( ); ++i )
			for( j = i; j < vecdP.size( ); ++j )
				Matrix.Set( i, j, ( i == j ) ? vecdP[ i ] : Matrix.Get( j, i ) );
		return true; }

	// Matrix operations
	static bool MatrixLUDecompose( CDataMatrix& Mat, std::vector<size_t>& veciIndices, bool& fEven );
	static bool MatrixLUInvert( CDataMatrix& MatLU, const std::vector<size_t>& veciIndices,
		CDataMatrix& MatInv );

	/*!
	 * \brief
	 * Calculate the determinant of a matrix given its LU decomposition.
	 * 
	 * \param MatLU
	 * LU decomposition of matrix of interest.
	 * 
	 * \param fEven
	 * True if LU decomposition of even rank; false otherwise.
	 * 
	 * \returns
	 * Determinant of original matrix.
	 * 
	 * \remarks
	 * MatLU must be square.  Implementation courtesy of Press WH, Teukolsky SA, Vetterling WT, Flannery BP.
	 * Numerical Recipes in C, 1992, Cambridge University Press.
	 * 
	 * \see
	 * MatrixLUDecompose
	 */
	template<class tType>
	static double MatrixLUDeterminant( const CFullMatrix<tType>& MatLU, bool fEven ) {
		double	dRet;
		size_t	i;

		if( MatLU.GetRows( ) != MatLU.GetColumns( ) )
			return CMeta::GetNaN( );

		dRet = fEven ? 1 : -1;
		for( i = 0; i < MatLU.GetRows( ); ++i )
			dRet *= MatLU.Get( i, i );

		return dRet; }
};

}

#endif // STATISTICS_H
