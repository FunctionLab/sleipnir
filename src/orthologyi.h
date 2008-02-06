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

	typedef std::map<CGene*,size_t>	TMapGeneI;

	~COrthologyImpl( );

	void Reset( );

	std::vector<std::string>			m_vecstrOrganisms;
	std::vector<CGenome*>				m_vecpGenomes;
	TMapGeneI							m_mapGenes;
	std::vector<std::vector<CGene*> >	m_vecvecpGenes;
};

}

#endif // ORTHOLOGYI_H
