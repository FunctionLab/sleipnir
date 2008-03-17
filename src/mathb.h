#ifndef MATH_H
#define MATH_H

#include "mathbi.h"

namespace Sleipnir {

/*!
 * \brief
 * Utility class containing static basic math functions.
 */
class CMath : CMathImpl {
public:
	static double LogFact( size_t iN );

	/*!
	 * \brief
	 * Calculates the sigmoid function with the given parameters.
	 * 
	 * \param dHeight
	 * Height (multiplier) of sigmoid.
	 * 
	 * \param dShift
	 * Horizontal shift (difference) of sigmoid.
	 * 
	 * \param dSlope
	 * Slope (sharpness) of sigmoid.
	 * 
	 * \param dVertical
	 * Vertical shift (constant addend) of sigmoid.
	 * 
	 * \param dX
	 * Point at which sigmoid should be evaluated.
	 * 
	 * \returns
	 * dVertical + (dHeight / (1 + exp(-dSlope * (dX - dShift))))
	 */
	static double Sigmoid( double dHeight, double dShift, double dSlope, double dVertical, double dX ) {

		return ( dVertical + ( dHeight / ( 1 + exp( -dSlope * ( dX - dShift ) ) ) ) ); }

	/*!
	 * \brief
	 * Return greatest common denominator of A and B.
	 * 
	 * \param iA
	 * First integer.
	 * 
	 * \param iB
	 * Second integer.
	 * 
	 * \returns
	 * Greatest common denominator of A and B.
	 * 
	 * \remarks
	 * Mmm, Euclid's algorithm.
	 */
	static size_t GCD( size_t iA, size_t iB ) {
		size_t	i;

		while( iB ) {
			i = iA;
			iA = iB;
			iB = i % iB; }

		return iA; }

	/*!
	 * \brief
	 * Return given value rounded to the nearest unsigned integer.
	 * 
	 * \param d
	 * Floating point value to round.
	 * 
	 * \returns
	 * Given value rounded to the nearest unsigned integer.
	 * 
	 * \remarks
	 * Calculated as the floor of the given value plus 0.5.  It would be peachy if Microsoft would get on
	 * top of C99 and implement round correctly...
	 */
	static size_t Round( double d ) {

		return (size_t)floor( d + 0.5 ); }
};

}

#endif // MATH_H
