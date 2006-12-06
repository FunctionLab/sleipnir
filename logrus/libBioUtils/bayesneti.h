#ifndef BAYESNETI_H
#define BAYESNETI_H

#include <smile.h>
#include <syscoord.h>

#include "dataset.h"
//#include "trie.h"

namespace pnl {

class CBNet;

}

namespace libBioUtils {

class CBayesNetPNL;
class CBayesNetSmile;
class CPCLPair;

class CBayesNetImpl {
protected:
	static const char	c_cMissing	= '_';
	static const char	c_cBase		= 'A';
	static const char	c_szFR[];
	static const char	c_szZero[];

	typedef std::vector<std::vector<float> >	TVecVecF;
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

class CBayesNetSmileImpl : protected CBayesNetImpl {
protected:
	friend class CBayesNetFN;

	static const size_t	c_iMinimum	= 10;
	static const char	c_szGaussian[];

	static bool IsGaussian( const DSL_network& );
	static bool IsNaive( const DSL_network& );
	static float ELRDot( const TVecVecF&, const TVecVecF& );
	static float ELRAvoidZero( float );
	static void ELRComputeDirection( float, const TVecVecF&, TVecVecF& );
	static bool GetCPT( DSL_node*, CDataMatrix& );

	CBayesNetSmileImpl( bool );

	bool ConvertGraph( CBayesNetPNL& ) const;
	bool ConvertCPTs( CBayesNetPNL& ) const;
	void LearnExpected( DSL_node*, DSL_Dmatrix*, size_t = 1 );
	bool Evaluate( const IDataset*, CDat*, TVecVecF*, bool ) const;
	bool IsContinuous( ) const;
	bool IsNaive( ) const;
	bool FillCPTs( const IDataset*, size_t, size_t, bool, bool );
	bool FillCPTs( const std::vector<bool>&, const std::string&, bool, bool, bool = false );
	bool FillCPTs( const std::vector<bool>&, const std::vector<unsigned char>&, bool, bool );
	bool LearnGrouped( const IDataset*, size_t, bool );
	bool LearnUngrouped( const IDataset*, size_t, bool );
	bool LearnNaive( const IDataset*, bool );
	bool LearnELR( const IDataset*, size_t, bool );
	size_t ELRCountParameters( ) const;
	void ELRCopyParameters( TVecVecF& );
	void ELRComputeGradient( const std::vector<bool>&, const TMapData&, bool, TVecVecF& );
	void ELRUpdateGradient( float, TVecVecF& );
	void ELRNormalizeDirection( TVecVecF& ) const;
	float ELRLineSearch( const std::vector<bool>&, const TMapData&, const TVecVecF&, const TVecVecF&, TVecVecF&,
		float&, float&, bool );
	float ELREvalFunction( const std::vector<bool>&, const TMapData&, float, const TVecVecF&, const TVecVecF&,
		TVecVecF&, bool );
	void ELRBracket( const std::vector<bool>&, const TMapData&, const TVecVecF&, const TVecVecF&, TVecVecF&,
		float&, float&, float&, float&, float&, float&, bool );
	float ELRConditionalLikelihood( const std::vector<bool>&, const TMapData&, bool );
	float ELRBrent( const std::vector<bool>&, const TMapData&, const TVecVecF&, const TVecVecF&, TVecVecF&,
		float&, float&, float, float, float, float, bool );

	bool					m_fSmileNet;
	DSL_network				m_SmileNet;
	const CBayesNetSmile*	m_pDefaults;
};

class CBayesNetPNLImpl : protected CBayesNetImpl {
protected:
	friend class CBayesNetSmileImpl;

	static const char	c_szBN[];

	CBayesNetPNLImpl( bool );
	~CBayesNetPNLImpl( );

	bool Evaluate( const IDataset*, CDat*, TVecVecF*, bool ) const;
	bool IsContinuous( ) const;

	pnl::CBNet*	m_pPNLNet;
};

}

#include "bayesnetfni.h"

#endif // BAYESNETI_H
