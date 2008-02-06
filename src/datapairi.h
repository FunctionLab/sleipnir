#ifndef DATAPAIRI_H
#define DATAPAIRI_H

#include "dat.h"

namespace libBioUtils {

class CPairImpl {
protected:
	static const char	c_szQuantExt[];

	static bool Open( const char*, std::ifstream& );
	static bool Open( const char*, std::vector<float>& );
};

class CDataPairImpl : protected CPairImpl, public CDat {
protected:
	void Reset( bool );

	bool				m_fContinuous;
	std::vector<float>	m_vecdQuant;
};

class CPCLPairImpl : protected CPairImpl, public CPCL {
protected:
	std::vector<std::vector<float> >	m_vecvecdQuants;
};

}

#endif // DATAPAIRI_H
