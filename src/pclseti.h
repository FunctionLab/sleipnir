#ifndef PCLSETI_H
#define PCLSETI_H

#include <map>

#include "fullmatrix.h"

namespace Sleipnir {

class CPCL;

class CPCLSetImpl {
protected:
	typedef std::map<string,size_t>	TMapStrI;

	CPCLSetImpl( );
	~CPCLSetImpl( );

	void Reset( );

	CPCL*					m_aPCLs;
	size_t					m_iPCLs;
	TMapStrI				m_mapGenes;
	std::vector<string>		m_vecstrGenes;
	CFullMatrix<size_t>		m_Genes;
};

}

#endif // PCLSETI_H
