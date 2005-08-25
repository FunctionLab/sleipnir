#ifndef CLUSTPIVOT_H
#define CLUSTPIVOT_H

#include <vector>

#include "halfmatrix.h"

namespace libBioUtils {

class CDat;

class CClustPivot {
public:
	static size_t Cluster( const CDistanceMatrix&, float, std::vector<size_t>& );
};

}

#endif // CLUSTPIVOT_H
