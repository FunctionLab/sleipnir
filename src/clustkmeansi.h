#ifndef CLUSTKMEANSI_H
#define CLUSTKMEANSI_H

#include <vector>

#include "fullmatrix.h"
#include "typesi.h"

namespace libBioUtils {

class CClustKMeansImpl {
protected:
	static void Randomize( CDataMatrix&, size_t, const CDataMatrix& );
};

}

#endif // CLUSTQTC_H
