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
#ifndef MEASURE_H
#define MEASURE_H

#include <float.h>
#include <math.h>

namespace Sleipnir {

/*!
 * \brief
 * Encapsulates any similarity (or occasionally distance) measure operating over two vectors.
 *
 * IMeasures are generally similarity measures consuming two floating point vectors and producing a single
 * continuous result.  There are a few distance measures implemented, but since everything in Sleipnir
 * assumes measures to be similarity-based (i.e. higher is closer), several IMeasures exist to reverse
 * distance measures (e.g. CMeasureNegate and CMeasureInvert).
 *
 * Most measures require input vectors of the same length, but some don't; likewise, most measures will
 * ignore weights for the vector elements if they're provided, but some won't.  Most measures will deal
 * appropriately with missing values (NaNs) in their inputs.  Many measures support different forms of
 * centering, although some will ignore centering and return the raw measure value.  For more information,
 * see EMap.
 *
 * \remarks
 * Most measure objects don't actually use any memory, but a few do, depending on their implementation,
 * and care should be taken to manage memory for these appropriately.  The Clone method exists to allow
 * a measure to copy itself without breaking the interface, and most measures that require memory
 * management take a constructor argument indicating whether they are responsible for freeing memory
 * references.
 */
class IMeasure {
public:
	/*!
	 * \brief
	 * Indicates how the result of a measure should be centered.
	 */
	enum EMap {
		/*!
		 * \brief
		 * Perform no centering and return the original measure result.
		 */
		EMapNone	= 0,
		/*!
		 * \brief
		 * Center the measure result by scaling to the range [0, 1].
		 */
		EMapCenter	= EMapNone + 1,
		/*!
		 * \brief
		 * Return the absolute value of the original measure result.
		 */
		EMapAbs		= EMapCenter + 1
	};

	virtual ~IMeasure( ) { };

	/*!
	 * \brief
	 * Return the human-readable unique identifier of the measure type.
	 *
	 * \returns
	 * The string identifier of the measure type.
	 */
	virtual const char* GetName( ) const = 0;
	/*!
	 * \brief
	 * Return true if the measure requires rank-based integer inputs.
	 *
	 * \returns
	 * True if the measure requires rank-based input vectors.
	 *
	 * \remarks
	 * The input vectors for rank measures should contain only floating point values with no fractional
	 * part; behavior is undefined if they don't.
	 */
	virtual bool IsRank( ) const = 0;
	/*!
	 * \brief
	 * Create a copy of the current measure object.
	 *
	 * \returns
	 * Copy of the current measure object.
	 *
	 * \remarks
	 * Caller is, of course, responsible for destroying the created object.
	 */
	virtual IMeasure* Clone( ) const = 0;
	/*!
	 * \brief
	 * Calculate the measure between two given vectors with optional element weights.
	 *
	 * \param adX
	 * First array of values.
	 *
	 * \param iN
	 * Length of first array.
	 *
	 * \param adY
	 * Second array of values.
	 *
	 * \param iM
	 * Length of second array.
	 *
	 * \param eMap
	 * Way in which returned value should be centered (implementation-specific).
	 *
	 * \param adWX
	 * If non-null, weights of elements in the first array.
	 *
	 * \param adWY
	 * If non-null, weights of elements in the second array.
	 *
	 * \returns
	 * Measure calculated between the two input vectors and, optionally, weights.
	 *
	 * \remarks
	 * Pretty much every implementation will puke if given bad input; bounds checking etc. is minimal.
	 */
	virtual double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const = 0;
};

}

#include "measurei.h"

namespace Sleipnir {

/*!
 * \brief
 * Inverts an underlying measure using a sigmoid function.
 *
 * Inverts and transforms the result of an underlying measure to the range [0, 1] by calculating:
 * \code
 * 2 * (1 - (1 / (1 + exp(-result * dMultiplier))))
 * \endcode
 *
 * \see
 * CMeasureNegate | CMeasureInvert
 */
class CMeasureSigmoid : CMeasureSigmoidImpl, public IMeasure {
public:
	/*!
	 * \brief
	 * Construct a new sigmoid measure wrapping the given underlying measure with the specified multiplier.
	 *
	 * \param pMeasure
	 * Measure whose result should be sigmoid transformed.
	 *
	 * \param fMemory
	 * If true, the new sigmoid measure is responsible for releasing the underlying measure's memory.
	 *
	 * \param dMultiplier
	 * Multiplier used when calculating the sigmoid function.
	 *
	 * \remarks
	 * If fMemory is true, the sigmoid measure will delete pMeasure when it is destroyed; otherwise, it will
	 * only hold a reference.
	 */
	CMeasureSigmoid( const IMeasure* pMeasure, bool fMemory, float dMultiplier ) :
		CMeasureSigmoidImpl( pMeasure, fMemory, dMultiplier ) { }

	const char* GetName( ) const {

		return m_pMeasure->GetName( ); }

	bool IsRank( ) const {

		return m_pMeasure->IsRank( ); }

	IMeasure* Clone( ) const {

		return new CMeasureSigmoid( m_pMeasure->Clone( ), true, m_dMult ); }

	double Measure( const float* adX, size_t iM, const float* adY, size_t iN, EMap eMap, const float* adWX,
		const float* adWY ) const {
		double	dRet;

		dRet = m_pMeasure->Measure( adX, iM, adY, iN, eMap, adWX, adWY );
		return ( 2 * ( 1 - ( 1 / ( 1 + exp( -dRet * m_dMult ) ) ) ) ); }
};

/*!
 * \brief
 * Inverts an underlying measure by negating its result.
 *
 * \see
 * CMeasureSigmoid | CMeasureInvert
 */
class CMeasureNegate : CMeasureImpl, public IMeasure {
public:
	/*!
	 * \brief
	 * Construct a new negation measure wrapping the given underlying measure.
	 *
	 * \param pMeasure
	 * Measure whose result should be negated.
	 *
	 * \param fMemory
	 * If true, the new negation measure is responsible for releasing the underlying measure's memory.
	 *
	 * \remarks
	 * If fMemory is true, the negation measure will delete pMeasure when it is destroyed; otherwise, it will
	 * only hold a reference.
	 */
	CMeasureNegate( const IMeasure* pMeasure, bool fMemory ) : CMeasureImpl( pMeasure, fMemory ) { }

	const char* GetName( ) const {

		return m_pMeasure->GetName( ); }

	bool IsRank( ) const {

		return m_pMeasure->IsRank( ); }

	IMeasure* Clone( ) const {

		return new CMeasureNegate( m_pMeasure->Clone( ), true ); }

	double Measure( const float* adX, size_t iM, const float* adY, size_t iN, EMap eMap, const float* adWX,
		const float* adWY ) const {

		return -m_pMeasure->Measure( adX, iM, adY, iN, eMap, adWX, adWY ); }
};

/*!
 * \brief
 * Inverts an underlying measure by inverting its result (dividing one by the value).
 *
 * \see
 * CMeasureSigmoid | CMeasureNegate
 */
class CMeasureInvert : CMeasureImpl, public IMeasure {
public:
	/*!
	 * \brief
	 * Construct a new inversion measure wrapping the given underlying measure.
	 *
	 * \param pMeasure
	 * Measure whose result should be inverted.
	 *
	 * \param fMemory
	 * If true, the new inversion measure is responsible for releasing the underlying measure's memory.
	 *
	 * \remarks
	 * If fMemory is true, the inversion measure will delete pMeasure when it is destroyed; otherwise, it will
	 * only hold a reference.
	 */
	CMeasureInvert( const IMeasure* pMeasure, bool fMemory ) : CMeasureImpl( pMeasure, fMemory ) { }

	const char* GetName( ) const {

		return m_pMeasure->GetName( ); }

	bool IsRank( ) const {

		return m_pMeasure->IsRank( ); }

	IMeasure* Clone( ) const {

		return new CMeasureInvert( m_pMeasure->Clone( ), true ); }

	double Measure( const float* adX, size_t iM, const float* adY, size_t iN, EMap eMap, const float* adWX,
		const float* adWY ) const {
		double	d;

		d = m_pMeasure->Measure( adX, iM, adY, iN, eMap, adWX, adWY );
		return ( d ? ( 1 / d ) : DBL_MAX ); }
};

/*!
 * \brief
 * Autocorrelates an underlying measure by rotating input vectors and returning the minimum result.
 *
 * An autocorrelation measure takes an underlying measure U and two input vectors [X0, ..., XN] and
 * [Y0, ..., YM] and returns the maximum of:
 * \code
 * U([X0, ..., XN], [Y0, ..., YM])
 * U([X0, ..., XN], [Y1, ..., YM, Y0])
 * U([X0, ..., XN], [Y2, ..., YM, Y0, Y1])
 * ...
 * U([X0, ..., XN], [YM, Y0, ..., YM-1])
 * \endcode
 * This is useful for dealing with periodic signals in data (e.g. cell cycle microarrays).
 */
class CMeasureAutocorrelate : CMeasureImpl, public IMeasure {
public:
	/*!
	 * \brief
	 * Construct a new autocorrelation measure wrapping the given underlying measure.
	 *
	 * \param pMeasure
	 * Measure whose result should be autocorrelated.
	 *
	 * \param fMemory
	 * If true, the new autocorrelation measure is responsible for releasing the underlying measure's memory.
	 *
	 * \remarks
	 * If fMemory is true, the autocorrelation measure will delete pMeasure when it is destroyed; otherwise,
	 * it will only hold a reference.
	 */
	CMeasureAutocorrelate( const IMeasure* pMeasure, bool fMemory ) : CMeasureImpl( pMeasure, fMemory ) { }

	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return m_pMeasure->GetName( ); }

	bool IsRank( ) const {

		return m_pMeasure->IsRank( ); }

	IMeasure* Clone( ) const {

		return new CMeasureAutocorrelate( m_pMeasure->Clone( ), true ); }
};

/*!
 * \brief
 * Calculates the Euclidean distance between the two vectors.
 *
 * Calculates Euclidean distance between two vectors; if weights are given, each pairwise product is also
 * multiplied by the appropriate elements' weights.  Centering is ignored.
 */
class CMeasureEuclidean : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "euclidean"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureEuclidean( ); }
};


class CMeasureEuclideanScaled : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "euclid_scaled"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureEuclidean( ); }
};


/*!
 * \brief
 * Calculates the Pearson correlation between the two vectors.
 *
 * Calculates Pearson correlation between two vectors; if weights are given, the means and each pairwise
 * product are also multiplied by the appropriate elements' weights.  Centering is performed as per EMap.
 *
 * \see
 * CMeasureQuickPearson | CMeasurePearNorm
 */
class CMeasurePearson : public IMeasure {
public:
	static double Pearson( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap,
		const float* adWX = NULL, const float* adWY = NULL, size_t* piCount = NULL );

	const char* GetName( ) const {

		return "pearson"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasurePearson( ); }

	double Measure( const float* adX, size_t iM, const float* adY, size_t iN, EMap eMap, const float* adWX,
		const float* adWY ) const {

		return CMeasurePearson::Pearson( adX, iM, adY, iN, eMap, adWX, adWY ); }
};

/*!
 * \brief
 * Calculates the Pearson correlation between the two vectors.
 *
 * Calculates a minimalistic Pearson correlation between two vectors for maximum efficiency: vectors must
 * be unweighted, of the same length, and contain no missing values.  Centering is performed as per EMap.
 *
 * \see
 * CMeasurePearson
 */
class CMeasureQuickPearson : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "quickpear"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureQuickPearson( ); }
};

/*!
 * \brief
 * Calculates the Kolmogorov-Smirnov p-value of difference between two vectors (centering and weights
 * ignored).
 *
 * \remarks
 * Low p-value indicates low similarity, since the vectors are probably different; the KS-test thus operates
 * as a reasonable similarity measure.
 */
class CMeasureKolmogorovSmirnov : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "kolm-smir"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureKolmogorovSmirnov( ); }
};

/*!
 * \brief
 * Calculates the Kendall's Tau correlation between two vectors (centering as per EMap, weights used in
 * pairwise products).
 */
class CMeasureKendallsTau : CMeasureKendallsTauImpl, public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "kendalls"; }

	bool IsRank( ) const {

		return true; }

	IMeasure* Clone( ) const {

		return new CMeasureKendallsTau( ); }
};

/*!
 * \brief
 * Calculates Spearman's rank correlation between two vectors (centering as per EMap, weights ignored).
 */
class CMeasureSpearman : CMeasureSpearmanImpl, public IMeasure {
public:
	/*!
	 * \brief
	 * Construct a new Spearman correlation measure with the indicated ranking behavior.
	 *
	 * \param fTransformed
	 * If true, all inputs are assumed to be pre-rank transformed; otherwise, rank-transforming is performed
	 * for each Measure call.
	 *
	 * \remarks
	 * If you're going to measure all pairwise combinations of some set (e.g. turn a PCL into a CDat), it is
	 * much more efficient to pre-transform all input vectors and construct a Spearman measure with
	 * fTransformed set to true.  This avoids repeatedly re-ranking each vector.
	 */
	CMeasureSpearman( bool fTransformed ) : CMeasureSpearmanImpl( fTransformed ) { }

	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "spearman"; }

	bool IsRank( ) const {

		return true; }

	IMeasure* Clone( ) const {

		return new CMeasureSpearman( m_fTransformed ); }
};

/*!
 * \brief
 * Calculates the Fisher's z-transformed Pearson correlation between the two vectors.
 *
 * Calculates Pearson correlation between two vectors and transforms the result using Fisher's z-transform.
 * This is done using the formula:
 * \code
 * log( (1 + dP) / (1 - dP) ) / 2
 * \endcode
 * If constructed with a known average and standard deviation, the resulting correlations will also be
 * z-scored (i.e. the resulting normal distribution of scores will be shifted to have mean zero and
 * standard deviation one).
 *
 * \see
 * CMeasurePearson
 */
class CMeasurePearNorm : CMeasurePearNormImpl, public IMeasure {
public:
	/*!
	 * \brief
	 * Construct a measure which will calculate Fisher's z-transformed Pearson correlations with no z-scoring.
	 */
	CMeasurePearNorm( ) : CMeasurePearNormImpl( HUGE_VAL, HUGE_VAL ) { }

	/*!
	 * \brief
	 * Construct a measure which will calculate z-scored Fisher's z-transformed Pearson correlations.
	 *
	 * \param dAverage
	 * Average used in z-scoring.
	 *
	 * \param dStdev
	 * Standard deviation used in z-scoring.
	 *
	 * After z-transformation, the z-score is calculated as:
	 * \code
	 * (dP - dAverage) / dStdev
	 * \endcode
	 */
	CMeasurePearNorm( double dAverage, double dStdev ) : CMeasurePearNormImpl( dAverage, dStdev ) { }

	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "pearnorm"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasurePearNorm( m_dAverage, m_dStdDev ); }
};

/*!
 * \brief
 * Calculate the hypergeometric p-value of overlap between two boolean vectors (centering and weights
 * ignored).
 *
 * Returns the hypergeometric p-value of overlap given two boolean vectors, counting the total number of
 * nonzero entries, the number of nonzero entries in each vector, and the number of positions in which
 * both vectors are nonzero.
 */
class CMeasureHypergeometric : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "hypergeom"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureHypergeometric( ); }
};

/*!
 * \brief
 * Calculates the inner product of two vectors (centering ignored, weights used in pairwise products).
 */
class CMeasureInnerProduct : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "innerprod"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureInnerProduct( ); }
};

/*!
 * \brief
 * Calculates the number of positions in which both vectors have nonzero elements (centering and weights
 * ignored).
 */
class CMeasureBinaryInnerProduct : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "bininnerprod"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureBinaryInnerProduct( ); }
};

/*!
 * \brief
 * Calculates the mutual information in bits between two integer vectors (centering and weights ignored).
 *
 * \remarks
 * Results will be pretty meaningless if the vectors contain non-integral values or are particularly
 * short.
 */
class CMeasureMutualInformation : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "mutinfo"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureMutualInformation( ); }
};

/*!
 * \brief
 * Calculates the difference in relative absolute sums of two vectors (centering and weights ignored).
 *
 * For two vectors [X0, ..., XN] and [Y0, ..., YN], calculates:
 * \code
 * 1 - sum(|Xi - Yi|)/sum(|Xi| + |Yi|)
 * \endcode
 * This is useful for detecting whether two vectors vary in approximately the same way (as per Pearson
 * correlation) and at approximately the same values (as per Euclidean distance).
 */
class CMeasureRelativeAUC : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "relauc"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureRelativeAUC( ); }
};

/*!
 * \brief
 * Calculates the p-value of Pearson correlation between two vectors (centering ignored, weights used for
 * correlation).
 *
 * \see
 * CMeasurePearson
 */
class CMeasurePearsonSignificance : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "pearsig"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasurePearsonSignificance( ); }
};

/*!
 * \brief
 * Calculates the continuous of distance correlation coefficient between two vectors, unweighted only!
 */
class CMeasureDistanceCorrelation : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "dcor"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureDistanceCorrelation( ); }
};


/*!
 * \brief
 * Calculates the continuous of distance correlation coefficient between two vectors, unweighted only!
 * Original version of distance correlation is unsigned. This signed version takes the sign from pearson correlation.
 */
class CMeasureSignedDistanceCorrelation : public IMeasure {
public:
	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapNone,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "sdcor"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureSignedDistanceCorrelation( ); }
};

/*!
 * \brief
 * Calculates the continuous of Dice coefficient between two vectors, dot(x, y) /
 * ( dot(x, y) + a*||x-y|| + (1-a)*||y-x|| ).
 */
class CMeasureDice : public IMeasure {
public:
	CMeasureDice( float dAlpha = 0.5 ) : m_dAlpha(dAlpha) { }

	double Measure( const float* adX, size_t iN, const float* adY, size_t iM, EMap eMap = EMapCenter,
		const float* adWX = NULL, const float* adWY = NULL ) const;

	const char* GetName( ) const {

		return "dice"; }

	bool IsRank( ) const {

		return false; }

	IMeasure* Clone( ) const {

		return new CMeasureDice( ); }

private:
	float	m_dAlpha;
};

}

#endif // MEASURE_H
