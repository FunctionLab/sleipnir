#ifndef CLUSTKMEANS_H
#define CLUSTKMEANS_H

#include "clustkmeansi.h"

namespace libBioUtils {

class IMeasure;

class CClustKMeans : protected CClustKMeansImpl {
public:
	static bool Cluster( const CDataMatrix&, const IMeasure*, size_t, std::vector<uint16_t>&,
		const CDataMatrix* = NULL );
};

}

#endif // CLUSTQTC_H
