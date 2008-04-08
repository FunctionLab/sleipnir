/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
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
