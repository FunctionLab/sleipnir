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
#ifndef STATISTICSI_H
#define STATISTICSI_H

#include "mathb.h"

namespace Sleipnir {

class CStatisticsImpl : protected CMath {
protected:
	static double Chi2CDF( double, double, double, double );
	static double IncompleteGamma( double, double );
	static double Normal01CDF( double );
	static double GammaLog( double );
	static double EpsilonDouble( );
	static double IncompleteBeta( double, double, double );
	static double IncompleteBetaCF( double, double, double );
	static double ModifiedBesselI( size_t, double );

	template<class tType>
	static void Sums( tType Begin, tType End, double* pdSum, double* pdSumSq ) {
		tType	Cur;

		if( pdSum )
			*pdSum = 0;
		if( pdSumSq )
			*pdSumSq = 0;
		for( Cur = Begin; Cur != End; ++Cur ) {
			if( pdSum )
				*pdSum += *Cur;
			if( pdSumSq )
				*pdSumSq += *Cur * *Cur; } }
};

}

#endif // STATISTICSI_H
