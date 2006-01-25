#ifndef BAYESNETI_H
#define BAYESNETI_H

#pragma warning (disable: 4244 4267)
#include <pnl_dll.hpp>
#pragma warning (default: 4244 4267)
#include <smile.h>

#include "dataset.h"

namespace libBioUtils {

class CBayesNetPNL;

class CBayesNetImpl {
protected:
	CBayesNetImpl( bool );

	static std::string EncodeDatum( const IDataset*, size_t, size_t );
	static void DecodeDatum( const std::string&, std::vector<size_t>& );

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
	bool LearnGrouped( const IDataset*, size_t, bool );
	bool LearnUngrouped( const IDataset*, size_t, bool );
	bool LearnNaive( const IDataset* );
	bool LearnELR( const IDataset*, size_t, bool );
	size_t ELRCountParameters( ) const;
	void ELRCopyParameters( TVecVecF& );
	void ELRComputeGradient( const IDataset*, bool, TVecVecF& );
	void ELRUpdateGradient( float, TVecVecF& );
	void ELRNormalizeDirection( TVecVecF& ) const;
	float ELRLineSearch( const IDataset*, const TVecVecF&, const TVecVecF&, TVecVecF&, float&, float&, bool );
	float ELREvalFunction( const IDataset*, float, const TVecVecF&, const TVecVecF&, TVecVecF&, bool );
	void ELRBracket( const IDataset*, const TVecVecF&, const TVecVecF&, TVecVecF&, float&, float&, float&, float&,
		float&, float&, bool );
	float ELRConditionalLikelihood( const IDataset*, bool );
	float ELRBrent( const IDataset*, const TVecVecF&, const TVecVecF&, TVecVecF&,
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
