#ifndef CLUSTPIVOT_H
#define CLUSTPIVOT_H

#include <vector>

#include "halfmatrix.h"
#include "typesi.h"

namespace Sleipnir {

class CDat;

/*!
 * \brief
 * Utility class containing static pivot clustering methods.
 */
class CClustPivot {
public:
	static uint16_t Cluster( const CDistanceMatrix& MatSimilarities, float dCutoff,
		std::vector<uint16_t>& vecsClusters );
};

}

#endif // CLUSTPIVOT_H
