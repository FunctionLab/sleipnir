#ifndef FASTAI_H
#define FASTAI_H

#include <map>
#include <string>
#include <vector>

#include "file.h"

namespace Sleipnir {

class CFASTAImpl : public CFile {
protected:
	static const char	c_acComment[];
	static const char	c_acHeader[];

	typedef std::map<std::string, size_t>		TMapStrI;
	typedef std::map<std::string, std::string>	TMapStrStr;

	CFASTAImpl( ) {

		m_szBuffer = new char[ c_iBufferSize ]; }

	virtual ~CFASTAImpl( ) {

		delete[] m_szBuffer; }

	mutable std::ifstream		m_ifsm;
	TMapStrI					m_mapstriGenes;
	std::vector<std::string>	m_vecstrGenes;
	std::vector<TMapStrStr>		m_vecmapstrstrHeaders;
	std::vector<TMapStrI>		m_vecmapstriSequences;
	char*						m_szBuffer;
};

}

#endif // FASTAI_H
