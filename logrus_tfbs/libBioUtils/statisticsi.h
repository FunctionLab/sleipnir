#ifndef STATISTICSI_H
#define STATISTICSI_H

#include "mathb.h"

namespace libBioUtils {

class CStatisticsImpl : protected CMath {
protected:
	static double Chi2CDF( double, double, double, double );
	static double IncompleteGamma( double, double );
	static double Normal01CDF( double );
	static double GammaLog( double );
	static double EpsilonDouble( );
	static double IncompleteBeta( double, double, double );
	static double IncompleteBetaCF( double, double, double );
};

}

#endif // STATISTICSI_H
