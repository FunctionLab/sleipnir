#ifndef ANNOTATIONI_H
#define ANNOTATIONI_H

#include <map>
#include <iostream>
#include <set>
#include <stack>
#include <vector>

#include "file.h"

namespace libBioUtils {

class CGene;
class CGenes;
class CGenome;
class IOntology;

class COntologyImpl {
protected:
	typedef std::map<std::string,size_t>	TMapStrI;
	typedef std::set<const CGene*>			TSetPGenes;

	struct SNode {
		SNode( );

		void Reset( );

		std::string		m_strID;
		std::string		m_strGloss;
		size_t			m_iParents;
		size_t*			m_aiParents;
		size_t			m_iChildren;
		size_t*			m_aiChildren;
		size_t			m_iGenes;
		const CGene**	m_apGenes;
		size_t			m_iCacheGenes;
		const CGene**	m_apCacheGenes;
	};

	struct SParser {
		static const size_t	c_iBuffer	= 4096;

		SParser( std::istream&, CGenome& );

		bool GetLine( );
		bool IsStart( const char* ) const;

		istream&	m_istm;
		CGenome&	m_Genome;
		char		m_szLine[ c_iBuffer ];
		std::string	m_strGloss;
	};

	COntologyImpl( const std::string& );
	~COntologyImpl( );

	size_t GetNode( const std::string& ) const;
	bool IsAnnotated( size_t, const CGene&, bool ) const;
	size_t GetNodes( ) const;
	size_t GetParents( size_t ) const;
	size_t GetParent( size_t, size_t ) const;
	size_t GetChildren( size_t ) const;
	size_t GetChild( size_t, size_t ) const;
	size_t GetGenes( size_t, bool ) const;
	const CGene& GetGene( size_t, size_t ) const;
	const std::string& GetID( size_t ) const;
	const std::string& GetID( ) const;
	const std::string& GetGloss( size_t ) const;
	void GetGeneNames( std::vector<std::string>& ) const;
	void Reset( );
	void CollectGenes( size_t ) const;
	void CollectGenes( size_t, TSetPGenes& );
	void TermFinder( const CGenes&, std::vector<TPrID>&, bool, bool, bool ) const;

	const IOntology*	m_pOntology;
	std::string			m_strID;
	size_t				m_iNodes;
	TMapStrI			m_mapNodes;
	SNode*				m_aNodes;
};

class COntologyKEGGImpl : protected COntologyImpl {
protected:
	static const char	c_szKEGG[];
	static const char	c_szEntry[];
	static const char	c_szName[];
	static const char	c_szDefinition[];
	static const char	c_szClass[];
	static const char	c_szPath[];
	static const char	c_szDBLinks[];
	static const char	c_szGenes[];
	static const char	c_szYeast[];
	static const char	c_szEnd[];
	static const size_t	c_iKEGG		= 10000;

	struct SParserKEGG : SParser {
		SParserKEGG( std::istream&, CGenome& );

		void Reset( );

		bool								m_fYeast;
		std::vector<CGene*>					m_vecpGenes;
		std::vector<std::string>			m_vecstrIDs;
		std::map<std::string,std::string>	m_mapGlosses;
	};

	COntologyKEGGImpl( );

	bool Open( SParserKEGG& );
	bool OpenEntry( SParserKEGG& );
	bool OpenName( SParserKEGG& );
	bool OpenDefinition( SParserKEGG& );
	bool OpenClass( SParserKEGG& );
	bool OpenDBLinks( SParserKEGG& );
	bool OpenGenes( SParserKEGG& );
	bool OpenOrganism( SParserKEGG& );
	char* OpenGene( SParserKEGG&, char* );
	bool OpenEnd( SParserKEGG& );
	bool OpenGloss( SParserKEGG& );
};

class COntologyGOImpl : protected COntologyImpl {
protected:
	static const char	c_szAltID[];
	static const char	c_szBP[];
	static const char	c_szBioProc[];
	static const char	c_szBroadSyn[];
	static const char	c_szCC[];
	static const char	c_szCelComp[];
	static const char	c_szComment[];
	static const char	c_szDef[];
	static const char	c_szExactSyn[];
	static const char	c_szGO[];
	static const char	c_szGOC[];
	static const char	c_szID[];
	static const char	c_szIsA[];
	static const char	c_szIsObsolete[];
	static const char	c_szMF[];
	static const char	c_szName[];
	static const char	c_szNarrowSyn[];
	static const char	c_szMolFunc[];
	static const char	c_szNamespace[];
	static const char	c_szNOT[];
	static const char	c_szPartOf[];
	static const char	c_szRelationship[];
	static const char	c_szRelatedSyn[];
	static const char	c_szSGD[];
	static const char	c_szSubset[];
	static const char	c_szSynonym[];
	static const char	c_szTerm[];
	static const char	c_szXRefAnalog[];
	static const char	c_szXRefUnknown[];

	struct SParserGO : SParser {
		typedef std::set<const CGene*>	TSetPGene;

		SParserGO( std::istream&, CGenome& );

		void Reset( );

		const char*					m_szTarget;
		bool						m_fObsolete;
		std::string					m_strNamespace;
		std::vector<std::string>	m_vecstrIDs;
		std::vector<size_t>			m_veciParents;
		std::vector<SNode>			m_vecNodes;
		std::vector<TSetPGene>		m_vecsetpGenes;
	};

	COntologyGOImpl( );

	bool OpenOntology( SParserGO& );
	bool OpenHeader( SParserGO& );
	bool OpenBlock( SParserGO& );
	bool OpenTerm( SParserGO& );
	bool OpenID( SParserGO& );
	bool OpenName( SParserGO& );
	bool OpenNamespace( SParserGO& );
	bool OpenRelationship( SParserGO& );
	bool OpenParent( SParserGO& );
	bool OpenAltID( SParserGO& );
	bool OpenObsolete( SParserGO& );
	bool OpenGenes( SParserGO& );
	bool OpenGene( SParserGO& );
};

class COntologyMIPSImpl : protected COntologyImpl {
protected:
	static const char	c_szMIPS[];

	struct SParserMIPS : SParser {
		SParserMIPS( std::istream&, CGenome& );

		std::vector<size_t>						m_veciParents;
		std::vector<std::string>				m_vecstrIDs;
		std::vector<std::string>				m_vecstrGlosses;
		std::stack<size_t>						m_stakiHier;
		std::vector<std::vector<const CGene*> >	m_vecpGenes;
	};

	COntologyMIPSImpl( );

	bool OpenOntology( SParserMIPS& );
	bool OpenCategory( SParserMIPS& );
	size_t OpenID( SParserMIPS& );
	bool OpenGenes( SParserMIPS& );
	bool OpenGene( SParserMIPS& );
};

class CSlimImpl : protected CFile {
protected:
	void Reset( const IOntology* );

	std::vector<std::string>				m_vecstrSlims;
	std::vector<std::vector<size_t> >		m_vecveciTerms;
	std::vector<std::vector<const CGene*> >	m_vecvecpGenes;
	const IOntology*						m_pOntology;
};

}

#endif // ANNOTATIONI_H
