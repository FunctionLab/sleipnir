#ifndef DATAPAIRI_H
#define DATAPAIRI_H

#include "dat.h"

namespace libBioUtils {

class CDataPairImpl : public CDat {
protected:
	static const char	c_szQuantExt[];

	void Reset( bool );

	std::vector<float>	m_vecdQuant;
	bool				m_fContinuous;
};

}

#endif // DATAPAIRI_H
