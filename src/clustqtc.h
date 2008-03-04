#ifndef CLUSTQTC_H
#define CLUSTQTC_H

#include "clustqtci.h"

namespace libBioUtils {

class CClustQTC : CClustQTCImpl {
public:
	static uint16_t Cluster( const CDataMatrix&, const IMeasure*, float, size_t, std::vector<uint16_t>&,
		const CDataMatrix* = NULL );
	static void Cluster( const CDataMatrix&, const IMeasure*, float, float, float, size_t, CDistanceMatrix&,
		const CDataMatrix* = NULL );
};

}

#endif // CLUSTQTC_H
