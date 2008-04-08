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
#ifndef GENOMEI_H
#define GENOMEI_H

#include <map>
#include <vector>

namespace Sleipnir {

class CGene;
class CGenome;
class IOntology;

class CGeneImpl {
protected:
	CGeneImpl( const std::string& );
	~CGeneImpl( );

	CGeneImpl& operator=( const CGeneImpl& );
	void IncrementOntologies( const IOntology* );

	std::string				m_strName;
	size_t					m_iOntologies;
	const IOntology**		m_apOntologies;
	std::vector<size_t>**	m_apveciAnnotations;
	size_t					m_iSynonyms;
	std::string*			m_astrSynonyms;
	bool					m_fRNA;
	bool					m_fDubious;
	std::string				m_strGloss;
};

class CGenomeImpl {
protected:
	typedef std::map<std::string,size_t>	TMapStrI;

	static const char	c_szDubious[];
	static const char	c_szORF[];
	static const char*	c_aszRNA[];

	~CGenomeImpl( );

	std::vector<CGene*>	m_vecpGenes;
	TMapStrI			m_mapGenes;
};

class CGenesImpl {
protected:
	static const char	c_cComment	= '#';

	typedef std::map<std::string,size_t>	TMapStrI;

	CGenesImpl( CGenome& );

	CGenome&			m_Genome;
	std::vector<CGene*>	m_vecpGenes;
	TMapStrI			m_mapGenes;
};

}

#endif // GENOMEI_H
