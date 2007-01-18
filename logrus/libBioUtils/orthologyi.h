#ifndef ORTHOLOGYI_H
#define ORTHOLOGYI_H

#include <map>

#include "file.h"

namespace libBioUtils {

class CGene;
class CGenome;

class COrthologyImpl : protected CFile {
protected:
	static const char	c_cOrgSep	= '|';

	typedef std::map<CGene*,std::vector<size_t> >	TMapGeneVecI;
	typedef std::map<std::string,CGenome*>			TMapStrGenome;

	~COrthologyImpl( );

	void Reset( );

	TMapStrGenome						m_mapGenomes;
	TMapGeneVecI						m_mapOrthology;
	std::vector<std::vector<CGene*> >	m_vecvecpGenes;
};

}

#endif // ORTHOLOGYI_H
