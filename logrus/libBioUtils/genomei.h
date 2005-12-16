#ifndef GENOMEI_H
#define GENOMEI_H

#include <map>
#include <vector>

namespace libBioUtils {

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
	typedef std::map<std::string,size_t>	TMapStrI;

	CGenesImpl( CGenome& );

	CGenome&			m_Genome;
	std::vector<CGene*>	m_vecpGenes;
	TMapStrI			m_mapGenes;
};

}

#endif // GENOMEI_H
