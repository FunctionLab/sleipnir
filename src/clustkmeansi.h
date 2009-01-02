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
#ifndef CLUSTKMEANSI_H
#define CLUSTKMEANSI_H

#include <vector>

#include "fullmatrix.h"
#include "halfmatrix.h"
#include "meta.h"

namespace Sleipnir {

class CClustKMeansImpl {
protected:
	static void Randomize( CDataMatrix&, size_t, const CDataMatrix& );

	static float GetClean( size_t iOne, size_t iTwo, float dMin, float dMax, const CDistanceMatrix& Mat ) {
		float	dRet;

		return ( ( iOne == iTwo ) ? dMax : ( CMeta::IsNaN( dRet = Mat.Get( iOne, iTwo ) ) ? dMin : dRet ) ); }
};

}

#endif // CLUSTQTC_H
