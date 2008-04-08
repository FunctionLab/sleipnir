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
#ifndef BAYESNETI_H
#define BAYESNETI_H

#ifndef NO_SMILE
#include <smile.h>
#include <syscoord.h>
#endif // NO_SMILE

#include "dataset.h"
//#include "trie.h"

namespace pnl {

class CBNet;

}

namespace Sleipnir {

class CBayesNetPNL;
class CBayesNetSmile;
class CPCLPair;

class CBayesNetImpl {
protected:
	static const char	c_cMissing	= '_';
	static const char	c_cBase		= 'A';
	static const char	c_szFR[];
	static const char	c_szZero[];

	typedef std::vector<std::vector<float> >	TVecVecD;
	typedef std::map<std::string, size_t>		TMapData;
//	typedef CTrie<size_t>						TTrieData;

	CBayesNetImpl( bool );

	static void EncodeData( const IDataset*, TMapData& );
	static std::string EncodeDatum( const IDataset*, size_t, size_t );
	static std::string EncodeDatum( const CPCLPair&, size_t, const std::vector<size_t>& );
	static void DecodeDatum( const std::string&, std::vector<size_t>& );
	static bool IsAnswer( const std::string& );

	bool	m_fGroup;
};

#ifndef NO_SMILE

class CBayesNetSmileImpl : protected CBayesNetImpl {
protected:
	friend class CBayesNetFN;

	static const size_t	c_iMinimum	= 25;
	static const char	c_szGaussian[];

	static bool IsGaussian( const DSL_network& );
	static bool IsNaive( const DSL_network& );
	static float ELRDot( const TVecVecD&, const TVecVecD& );
	static float ELRAvoidZero( float );
	static void ELRComputeDirection( float, const TVecVecD&, TVecVecD& );
	static bool GetCPT( DSL_node*, CDataMatrix& );

	CBayesNetSmileImpl( bool );

	bool ConvertGraph( CBayesNetPNL& ) const;
	bool ConvertCPTs( CBayesNetPNL& ) const;
	void LearnExpected( DSL_node*, DSL_Dmatrix*, size_t = 1 );
	bool Evaluate( const IDataset*, CDat*, TVecVecD*, bool ) const;
	bool FillCPTs( const IDataset*, size_t, size_t, bool, bool );
	bool FillCPTs( const std::vector<bool>&, const std::string&, bool, bool, bool = false );
	bool FillCPTs( const std::vector<bool>&, const std::vector<unsigned char>&, bool, bool, bool = false );
	bool LearnGrouped( const IDataset*, size_t, bool );
	bool LearnUngrouped( const IDataset*, size_t, bool );
	bool LearnNaive( const IDataset*, bool );
	bool LearnELR( const IDataset*, size_t, bool );
	size_t ELRCountParameters( ) const;
	void ELRCopyParameters( TVecVecD& );
	void ELRComputeGradient( const std::vector<bool>&, const TMapData&, bool, TVecVecD& );
	void ELRUpdateGradient( float, TVecVecD& );
	void ELRNormalizeDirection( TVecVecD& ) const;
	float ELRLineSearch( const std::vector<bool>&, const TMapData&, const TVecVecD&, const TVecVecD&, TVecVecD&,
		float&, float&, bool );
	float ELREvalFunction( const std::vector<bool>&, const TMapData&, float, const TVecVecD&, const TVecVecD&,
		TVecVecD&, bool );
	void ELRBracket( const std::vector<bool>&, const TMapData&, const TVecVecD&, const TVecVecD&, TVecVecD&,
		float&, float&, float&, float&, float&, float&, bool );
	float ELRConditionalLikelihood( const std::vector<bool>&, const TMapData&, bool );
	float ELRBrent( const std::vector<bool>&, const TMapData&, const TVecVecD&, const TVecVecD&, TVecVecD&,
		float&, float&, float, float, float, float, bool );

	bool IsNaive( ) const {

		return ( m_fSmileNet ? CBayesNetSmileImpl::IsNaive( m_SmileNet ) : false ); }

	bool IsContinuous( ) const {

		return ( m_fSmileNet ? IsGaussian( m_SmileNet ) : false ); }

	bool					m_fSmileNet;
	DSL_network				m_SmileNet;
	const CBayesNetSmile*	m_pDefaults;
};

#endif // NO_SMILE

class CBayesNetPNLImpl : protected CBayesNetImpl {
protected:
	friend class CBayesNetSmileImpl;

	static const char	c_szBN[];

	CBayesNetPNLImpl( bool );
	~CBayesNetPNLImpl( );

	bool Evaluate( const IDataset*, CDat*, TVecVecD*, bool ) const;
	bool IsContinuous( ) const;

	pnl::CBNet*	m_pPNLNet;
};

}

#include "bayesnetfni.h"

#endif // BAYESNETI_H
