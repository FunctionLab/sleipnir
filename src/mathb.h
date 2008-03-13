#ifndef MATH_H
#define MATH_H

#include "mathbi.h"

namespace Sleipnir {

class CMath : CMathImpl {
public:
	static size_t GCD( size_t, size_t );
	static double Min( double, double );
	static size_t Round( double );
	static double LogFact( size_t );

	static double Sigmoid( double dHeight, double dShift, double dSlope, double dVertical, double dX ) {

		return ( dVertical + ( dHeight / ( 1 + exp( -dSlope * ( dX - dShift ) ) ) ) ); }
};

}

#endif // MATH_H
