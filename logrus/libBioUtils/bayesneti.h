#ifndef BAYESNETI_H
#define BAYESNETI_H

#pragma warning (disable: 4244 4267)
#include <pnl_dll.hpp>
#pragma warning (default: 4244 4267)
#include <smile.h>

#include "dataset.h"
#include "trie.h"

namespace libBioUtils {

class CBayesNetPNL;

class CBayesNetImpl {
protected:
	static const char	c_cMissing	= '_';
	static const char	c_cBase		= 'A';

	typedef std::map<std::string, size_t>	TMapData;
	typedef CTrie<size_t>					TTrieData;

	CBayesNetImpl( bool );

	static void EncodeData( const IDataset*, TMapData& );
	static std::string EncodeDatum( const IDataset*, size_t, size_t );
	static void DecodeDatum( const std::string&, std::vector<size_t>& );
	static bool IsAnswer( const std::string& );

	bool	m_fGroup;
};

class CBayesNetSmileImpl : protected CBayesNetImpl {
protected:
	typedef std::vector<std::vector<float> >	TVecVecF;

	static const char	c_szGaussian[];

	static bool IsGaussian( const DSL_network& );
	static float ELRDot( const TVecVecF&, const TVecVecF& );
	static float ELRAvoidZero( float );
	static void ELRComputeDirection( float, const TVecVecF&, TVecVecF& );

	CBayesNetSmileImpl( bool );

	bool ConvertGraph( CBayesNetPNL& ) const;
	bool ConvertCPTs( CBayesNetPNL& ) const;
	void LearnExpected( DSL_node*, DSL_Dmatrix*, size_t = 1 );
	bool Evaluate( const IDataset*, CDat*, vector<vector<float> >*, bool ) const;
	bool IsContinuous( ) const;
	bool IsNaive( ) const;
	bool FillCPTs( const IDataset*, size_t, size_t, bool, bool );
	bool FillCPTs( const IDataset*, const std::string&, bool, bool );
	bool FillCPTs( const IDataset*, const std::vector<unsigned char>&, bool, bool );
	bool LearnGrouped( const IDataset*, size_t, bool );
	bool LearnUngrouped( const IDataset*, size_t, bool );
	bool LearnNaive( const IDataset* );
	bool LearnELR( const IDataset*, size_t, bool );
	size_t ELRCountParameters( ) const;
	void ELRCopyParameters( TVecVecF& );
	void ELRComputeGradient( const IDataset*, const TMapData&, bool, TVecVecF& );
	void ELRUpdateGradient( float, TVecVecF& );
	void ELRNormalizeDirection( TVecVecF& ) const;
	float ELRLineSearch( const IDataset*, const TMapData&, const TVecVecF&, const TVecVecF&, TVecVecF&,
		float&, float&, bool );
	float ELREvalFunction( const IDataset*, const TMapData&, float, const TVecVecF&, const TVecVecF&,
		TVecVecF&, bool );
	void ELRBracket( const IDataset*, const TMapData&, const TVecVecF&, const TVecVecF&, TVecVecF&,
		float&, float&, float&, float&, float&, float&, bool );
	float ELRConditionalLikelihood( const IDataset*, const TMapData&, bool );
	float ELRBrent( const IDataset*, const TMapData&, const TVecVecF&, const TVecVecF&, TVecVecF&,
		float&, float&, float, float, float, float, bool );

	bool		m_fSmileNet;
	DSL_network	m_SmileNet;
};

class CBayesNetPNLImpl : protected CBayesNetImpl {
protected:
	friend class CBayesNetSmileImpl;

	static const char	c_szBN[];

	CBayesNetPNLImpl( bool );
	~CBayesNetPNLImpl( );

	bool Evaluate( const IDataset*, CDat*, vector<vector<float> >*, bool ) const;
	bool IsContinuous( ) const;

	pnl::CBNet*	m_pPNLNet;
};

}

#endif // BAYESNETI_H
