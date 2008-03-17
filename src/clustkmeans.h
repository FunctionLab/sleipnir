#ifndef CLUSTKMEANS_H
#define CLUSTKMEANS_H

#include "clustkmeansi.h"

namespace Sleipnir {

class IMeasure;

/*!
 * \brief
 * Utility class containing static k-means clustering methods.
 */
class CClustKMeans : protected CClustKMeansImpl {
public:
	static bool Cluster( const CDataMatrix& MatData, const IMeasure* pMeasure, size_t iK,
		std::vector<uint16_t>& vecsClusters, const CDataMatrix* pMatWeights = NULL );
};

}

#endif // CLUSTQTC_H
