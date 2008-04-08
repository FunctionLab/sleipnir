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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef CLUSTQTC_H
#define CLUSTQTC_H

#include "clustqtci.h"

namespace Sleipnir {

/*!
 * \brief
 * Utility class containing static quality threshold clustering methods.
 */
class CClustQTC : CClustQTCImpl {
public:
	static uint16_t Cluster( const CDataMatrix& MatData, const IMeasure* pMeasure, float dDiameter,
		size_t iSize, std::vector<uint16_t>& vecsClusters, const CDataMatrix* pMatWeights = NULL );
	static void Cluster( const CDataMatrix& MatData, const IMeasure* pMeasure, float dMinDiameter,
		float dMaxDiameter, float dDeltaDiameter, size_t iSize, CDistanceMatrix& MatResults,
		const CDataMatrix* pMatWeights = NULL );
};

}

#endif // CLUSTQTC_H
