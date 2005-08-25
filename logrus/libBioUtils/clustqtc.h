#ifndef CLUSTQTC_H
#define CLUSTQTC_H

#include "clustqtci.h"

namespace libBioUtils {

class CClustQTC : CClustQTCImpl {
public:
	static size_t Cluster( const CDataMatrix&, const IMeasure*, float, size_t, bool,
		std::vector<size_t>& );
};

}

#endif // CLUSTQTC_H
