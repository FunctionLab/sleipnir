#ifndef MATH_H
#define MATH_H

#include "mathbi.h"

namespace libBioUtils {

class CMath : CMathImpl {
public:
	static size_t GCD( size_t, size_t );
	static double Min( double, double );
	static size_t Round( double );
	static double LogFact( size_t );
};

}

#endif // MATH_H
