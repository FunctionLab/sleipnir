#ifndef CLUSTQTCI_H
#define CLUSTQTCI_H

#include <vector>

#include "fullmatrix.h"
#include "halfmatrix.h"

namespace libBioUtils {

class CPCL;
class IMeasure;

class CClustQTCImpl {
protected:
	static void InitializeDistances( const CDataMatrix&, const IMeasure*, bool,
		CDistanceMatrix&, const CDataMatrix* );
	static double GetJackDistance( const float*, const float*, size_t, bool, float*, float*,
		const IMeasure*, const float*, const float*, float*, float* );
	static uint16_t QualityThresholdAll( const CDataMatrix&, float, size_t,
		const CDistanceMatrix&, std::vector<uint16_t>& );
	static void QualityThresholdLargest( const CDataMatrix&, float, const CDistanceMatrix&,
		const std::vector<bool>&, std::vector<uint16_t>& );
	static void QualityThresholdGene( size_t, const CDataMatrix&, float,
		const CDistanceMatrix&, std::vector<bool>&, std::vector<float>&,
		std::vector<uint16_t>& );
};

}

#endif // CLUSTQTCI_H
