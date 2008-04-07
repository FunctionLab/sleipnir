#ifndef PCLI_H
#define PCLI_H

#include <set>
#include <vector>

#include "file.h"
#include "fullmatrix.h"

namespace Sleipnir {

class CPCLImpl : protected CFile {
protected:
	static const size_t	c_iSkip			= 2;
	static const char	c_szEWEIGHT[];
	static const char	c_szGENE[];
	static const char	c_szGID[];
	static const char	c_szGWEIGHT[];
	static const char	c_szNAME[];
	static const char	c_szOne[];

	typedef std::vector<std::string>	TVecStr;
	typedef std::set<size_t>			TSetI;

	CPCLImpl( bool fHeader ) : m_fHeader(fHeader) { }
	~CPCLImpl( );

	bool OpenExperiments( std::istream&, size_t, char*, size_t );
	bool OpenGene( std::istream&, std::vector<float>&, char*, size_t );
	void Reset( );

	CDataMatrix				m_Data;
	TVecStr					m_vecstrGenes;
	TVecStr					m_vecstrExperiments;
	TVecStr					m_vecstrFeatures;
	std::vector<TVecStr>	m_vecvecstrFeatures;
	TSetI					m_setiGenes;
	bool					m_fHeader;
};

}

#endif // PCLI_H