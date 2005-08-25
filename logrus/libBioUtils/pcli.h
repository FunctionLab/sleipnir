#ifndef PCLI_H
#define PCLI_H

#include <set>
#include <vector>

#include "file.h"
#include "fullmatrix.h"

namespace libBioUtils {

class CPCLImpl : CFile {
protected:
	static const size_t	c_iSkip			= 2;
	static const char	c_szEWEIGHT[];
	static const char	c_szGENE[];
	static const char	c_szGID[];

	typedef std::vector<std::string>	TVecStr;
	typedef std::set<size_t>			TSetI;

	~CPCLImpl( );

	void Reset( );
	bool OpenExperiments( std::istream&, size_t );
	bool OpenGene( std::istream&, std::vector<float>& );

	CDataMatrix				m_Data;
	TVecStr					m_vecstrGenes;
	TVecStr					m_vecstrExperiments;
	TVecStr					m_vecstrFeatures;
	std::vector<TVecStr>	m_vecvecstrFeatures;
	TSetI					m_setiGenes;
};

}

#endif // PCLI_H
