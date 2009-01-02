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
#ifndef COALESCESTRUCTSI_H
#define COALESCESTRUCTSI_H

#ifdef _MSC_VER
#include <hash_map>
#else // _MSC_VER
#include <ext/hash_map>

#define hash_value	hash<const char*>( )
#define stdext		__gnu_cxx
#endif // _MSC_VER

#include "coalescebasei.h"
#include "fastai.h"
#include "fullmatrix.h"

namespace Sleipnir {

class CCoalesceMotifLibrary;
class CFASTA;
class CHierarchy;
class CPCL;

struct SCoalesceModifiers {

	void Initialize( const CPCL& );

	bool Add( const CFASTA* pWiggle ) {

		if( !pWiggle )
			return false;
		m_vecpWiggles.push_back( pWiggle );
		return true; }

	std::vector<const CFASTA*>			m_vecpWiggles;
	std::vector<std::vector<size_t> >	m_vecveciPCL2Wiggles;
};

struct SCoalesceModifierCache {
	SCoalesceModifierCache( const SCoalesceModifiers& Modifiers ) : m_Modifiers(Modifiers) { }

	void Get( size_t );
	void InitializeWeight( size_t, size_t );
	void AddWeight( size_t, size_t, size_t );
	void SetType( const std::string& );

	float GetWeight( size_t iK ) const {
		size_t	i, iWiggles;

		for( iWiggles = i = 0; i < m_veciWiggleTypes.size( ); ++i )
			if( m_veciWiggleTypes[ i ] != -1 )
				iWiggles++;

		return ( iWiggles ? ( m_dWeight / iWiggles ) : iK ); }

	const SCoalesceModifiers&				m_Modifiers;
	std::vector<std::vector<SFASTAWiggle> >	m_vecvecsWiggles;
	std::vector<size_t>						m_veciWiggleTypes;
	float									m_dWeight;
};

struct SMotifMatch {
	SMotifMatch( ) { }

	SMotifMatch( uint32_t iMotif, const std::string& strType,
		CCoalesceSequencerBase::ESubsequence eSubsequence, float dZ, float dAverage ) : m_iMotif(iMotif),
		m_strType(strType), m_eSubsequence(eSubsequence), m_dZ(dZ), m_dAverage(dAverage) { }

	bool Open( std::istream&, CCoalesceMotifLibrary& );
	uint32_t Open( const CHierarchy&, const std::vector<SMotifMatch>&, CCoalesceMotifLibrary&, size_t& );
	uint32_t Open( const SMotifMatch&, const SMotifMatch&, CCoalesceMotifLibrary& );
	std::string Save( const CCoalesceMotifLibrary*, bool = false, float = 0, float = 0, float = 0,
		bool = false ) const;

	bool operator==( const SMotifMatch& sMotif ) const {

		return ( ( m_iMotif == sMotif.m_iMotif ) && ( m_strType == sMotif.m_strType ) &&
			( m_eSubsequence == sMotif.m_eSubsequence ) ); }

	bool operator!=( const SMotifMatch& sMotif ) const {

		return !( *this == sMotif ); }

	bool operator<( const SMotifMatch& sMotif ) const {

		if( m_iMotif == sMotif.m_iMotif ) {
			if( m_strType == sMotif.m_strType )
				return ( m_eSubsequence < sMotif.m_eSubsequence );
			return ( m_strType < sMotif.m_strType ); }

		return ( m_iMotif < sMotif.m_iMotif ); }

	size_t GetHash( ) const {
		size_t	iMotif, iType, iSubsequence;

		iMotif = m_iMotif * ( (size_t)-1 / 20000 );
		iType = stdext::hash_value( m_strType.c_str( ) ) * ( (size_t)-1 / 5 );
		iSubsequence = m_eSubsequence * ( (size_t)-1 / CCoalesceSequencerBase::ESubsequenceEnd );

		return ( iMotif ^ iType ^ iSubsequence ); }

	uint32_t								m_iMotif;
	std::string								m_strType;
	CCoalesceSequencerBase::ESubsequence	m_eSubsequence;
	float									m_dZ;
	float									m_dAverage;
};

struct SCoalesceDataset {
	template<class tType>
	SCoalesceDataset( const tType& Conditions ) {

		m_veciConditions.resize( Conditions.size( ) );
		std::copy( Conditions.begin( ), Conditions.end( ), m_veciConditions.begin( ) ); }

	SCoalesceDataset( size_t iCondition ) {

		m_veciConditions.resize( 1 );
		m_veciConditions[ 0 ] = iCondition; }

	bool CalculateCovariance( const CPCL& );

	bool IsCondition( size_t iCondition ) const {

		return ( std::find( m_veciConditions.begin( ), m_veciConditions.end( ), iCondition ) !=
			m_veciConditions.end( ) ); }

	size_t GetConditions( ) const {

		return m_veciConditions.size( ); }

	size_t GetCondition( size_t iCondition ) const {

		return m_veciConditions[ iCondition ]; }

	std::vector<size_t>	m_veciConditions;
	CDataMatrix			m_MatSigmaChol;
	CDataMatrix			m_MatSigmaInv;
	double				m_dSigmaDetSqrt;
	std::vector<float>	m_vecdStdevs;
};

}

#endif // COALESCESTRUCTSI_H
