#ifndef CLUSTPIVOT_H
#define CLUSTPIVOT_H

#include <vector>

#include "halfmatrix.h"
#include "typesi.h"

namespace Sleipnir {

class CDat;

class CClustPivot {
public:
	static uint16_t Cluster( const CDistanceMatrix&, float, std::vector<uint16_t>& );
};

}

#endif // CLUSTPIVOT_H
