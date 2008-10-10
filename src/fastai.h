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
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef FASTAI_H
#define FASTAI_H

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "file.h"

namespace Sleipnir {

struct SFASTABase {
	std::string	m_strType;
	bool		m_fIntronFirst;
};


struct SFASTASequence : SFASTABase {
	std::vector<std::string>	m_vecstrSequences;
};

struct SFASTAWiggle : SFASTABase {
	std::vector<float>	m_vecdValues;
};

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

	bool Get( size_t, std::vector<SFASTASequence>*, std::vector<SFASTAWiggle>* ) const;
	bool Get( size_t, std::vector<SFASTASequence>&, size_t, const std::string&, SFASTASequence& ) const;
	bool Get( size_t, std::vector<SFASTAWiggle>&, size_t, SFASTAWiggle& ) const;

	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	mutable std::ifstream		m_ifsm;
	TMapStrI					m_mapstriGenes;
	std::vector<std::string>	m_vecstrGenes;
	std::vector<TMapStrStr>		m_vecmapstrstrHeaders;
	std::vector<TMapStrI>		m_vecmapstriSequences;
	char*						m_szBuffer;
	std::set<std::string>		m_setstrTypes;
};

}

#endif // FASTAI_H
