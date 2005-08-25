#ifndef MATH_H
#define MATH_H

#include "mathbi.h"

namespace libBioUtils {

class CMath : CMathImpl {
public:
	static unsigned int GCD( unsigned int, unsigned int );
	static double Min( double, double );
	static size_t Round( double );
	static double LogFact( size_t );
};

}

#endif // MATH_H
