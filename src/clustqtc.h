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
