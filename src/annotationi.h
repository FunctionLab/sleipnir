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
#ifndef ANNOTATIONI_H
#define ANNOTATIONI_H

#include <map>
#include <iostream>
#include <set>
#include <stack>
#include <string>
#include <vector>

#include "file.h"

namespace Sleipnir {

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
		static const size_t	c_iBuffer	= 65536;

		SParser( std::istream&, CGenome& );

		bool GetLine( );
		bool IsStart( const char* ) const;

		std::istream&	m_istm;
		CGenome&		m_Genome;
		char			m_szLine[ c_iBuffer ];
		std::string		m_strGloss;
		size_t			m_iLine;
	};

	COntologyImpl( const std::string& strID ) : m_strID(strID), m_iNodes(0), m_aNodes(NULL) { }

	~COntologyImpl( ) {

		Reset( ); }

	size_t GetNode( const std::string& ) const;
	bool IsAnnotated( size_t, const CGene&, bool ) const;
	const CGene& GetGene( size_t, size_t ) const;
	void GetGeneNames( std::vector<std::string>& ) const;
	void Reset( );
	void CollectGenes( size_t, TSetPGenes& );
	void TermFinder( const CGenes&, std::vector<STermFound>&, bool, bool, bool, float, const CGenes* ) const;

	size_t GetNodes( ) const {

		return m_iNodes; }

	size_t GetParents( size_t iNode ) const {

		return m_aNodes[ iNode ].m_iParents; }

	size_t GetParent( size_t iNode, size_t iParent ) const {

		return m_aNodes[ iNode ].m_aiParents[ iParent ]; }

	size_t GetChildren( size_t iNode ) const {

		return m_aNodes[ iNode ].m_iChildren; }

	size_t GetChild( size_t iNode, size_t iChild ) const {

		return m_aNodes[ iNode ].m_aiChildren[ iChild ]; }

	size_t GetGenes( size_t iNode, bool fKids ) const {
		size_t	iRet;

		iRet = m_aNodes[ iNode ].m_iGenes;
		if( fKids ) {
			CollectGenes( iNode );
			iRet += m_aNodes[ iNode ].m_iCacheGenes; }

		return iRet; }

	const std::string& GetID( ) const {

		return m_strID; }

	const std::string& GetID( size_t iNode ) const {

		return m_aNodes[ iNode ].m_strID; }

	const std::string& GetGloss( size_t iNode ) const {

		return m_aNodes[ iNode ].m_strGloss; }

	void CollectGenes( size_t iNode ) const {
		TSetPGenes	setpGenes;

		if( m_aNodes[ iNode ].m_iCacheGenes == -1 )
			((COntologyImpl*)this)->CollectGenes( iNode, setpGenes ); }

	bool GetChildren( size_t iNode, std::set<size_t>& setiChildren ) const {
		size_t	i, iChild;

		if( setiChildren.find( iNode ) != setiChildren.end( ) )
			return true;

		for( i = 0; i < GetChildren( iNode ); ++i ) {
			if( !GetChildren( iChild = GetChild( iNode, i ), setiChildren ) )
				return false;
			setiChildren.insert( iChild ); }

		return true; }

	bool GetParents( size_t iNode, std::set<size_t>& setiParents ) const {
		size_t	i, iParent;

		if( setiParents.find( iNode ) != setiParents.end( ) )
			return true;

		for( i = 0; i < GetParents( iNode ); ++i ) {
			if( !GetParents( iParent = GetParent( iNode, i ), setiParents ) )
				return false;
			setiParents.insert( iParent ); }

		return true; }

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
	static const char	c_szReference[];
	static const char	c_szDisease[];
	static const char	c_szPathway[];
	static const char	c_szModule[];
	static const char	c_szBR[];
	static const char	c_szDBLinks[];
	static const char	c_szGenes[];
	static const char	c_szEnd[];
	static const size_t	c_iKEGG		= 10000;

	struct SParserKEGG : SParser {
		SParserKEGG( std::istream&, CGenome&, const std::string&, bool fSynonyms );

		void Reset( );

		const std::string&					m_strOrganism;
		bool								m_fOrganism;
		bool								m_fPathing;
		bool								m_fSynonyms;
		std::vector<CGene*>					m_vecpGenes;
		std::vector<std::string>			m_vecstrIDs;
		std::map<std::string,std::string>	m_mapGlosses;
	};

	COntologyKEGGImpl( );

	bool Open( SParserKEGG& );
	bool OpenEntry( SParserKEGG& );
	bool OpenReferences( SParserKEGG& );
	bool OpenReference( SParserKEGG& );
	bool OpenName( SParserKEGG& );
	bool OpenDisease( SParserKEGG& );
	bool OpenPathway( SParserKEGG& );
	bool OpenModule( SParserKEGG& );
	bool OpenDefinition( SParserKEGG& );
	bool OpenClass( SParserKEGG& );
	bool OpenDBLinks( SParserKEGG& );
	bool OpenGenes( SParserKEGG& );
	bool OpenOrganism( SParserKEGG& );
	char* OpenGene( SParserKEGG&, char* );
	bool OpenEnd( SParserKEGG& );
	bool OpenGloss( SParserKEGG& );
};

class COntologyOBOImpl : protected COntologyImpl {
protected:
	static const char	c_szAltID[];
	static const char	c_szOBO[];
	static const char	c_szHUMAN[];
	static const char	c_szID[];
	static const char	c_szIsA[];
	static const char	c_szIsObsolete[];
	static const char	c_szName[];
	static const char	c_szNamespace[];
	static const char	c_szNOT[];
	static const char	c_szPartOf[];
	static const char	c_szRelationship[];
	static const char	c_szSGD[];
	static const char	c_szTerm[];
	
	struct SParserOBO : SParser {
		typedef std::set<const CGene*>	TSetPGene;

		SParserOBO( std::istream&, CGenome&, bool = false, bool = false );

		void Reset( );

		const char*					m_szTarget;
		std::vector<std::vector<std::string> >	m_vecvecstrParents;
		bool						m_fObsolete;
		bool						m_fDBIDs;
		bool						m_fSynonyms;
		std::string					m_strNamespace;
		std::vector<std::string>	m_vecstrIDs;
		std::vector<SNode>			m_vecNodes;
		std::vector<TSetPGene>		m_vecsetpGenes;
	};

	COntologyOBOImpl( );
	
	bool OpenOntology( SParserOBO& );
	bool OpenHeader( SParserOBO& );
	bool OpenBlock( SParserOBO& );
	bool OpenTerm( SParserOBO& );
	bool OpenID( SParserOBO& );
	bool OpenName( SParserOBO& );
	bool OpenNamespace( SParserOBO& );
	bool OpenRelationship( SParserOBO& );
	bool OpenParent( SParserOBO& );
	bool OpenAltID( SParserOBO& );
	bool OpenObsolete( SParserOBO& );
	bool OpenGenes( SParserOBO& );
	bool OpenGene( SParserOBO& );
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
