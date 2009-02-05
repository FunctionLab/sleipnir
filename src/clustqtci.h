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
#ifndef CLUSTQTCI_H
#define CLUSTQTCI_H

#include <vector>

#include "fullmatrix.h"
#include "halfmatrix.h"

namespace Sleipnir {

class CPCL;
class IMeasure;

class CClustQTCImpl {
protected:
	static void InitializeDistances( const CDataMatrix&, const IMeasure*, CDistanceMatrix&,
		const CDataMatrix* );
	static double GetJackDistance( const float*, const float*, size_t, float*, float*, const IMeasure*,
		const float*, const float*, float*, float* );
	static uint16_t QualityThresholdAll( size_t, float, size_t, const CDistanceMatrix&,
		std::vector<uint16_t>& );
	static void QualityThresholdLargest( size_t, float, const CDistanceMatrix&, const std::vector<bool>&,
		std::vector<uint16_t>& );
	static void QualityThresholdGene( size_t, size_t, float, const CDistanceMatrix&, std::vector<bool>&,
		std::vector<float>&, std::vector<uint16_t>& );
};

}

#endif // CLUSTQTCI_H
