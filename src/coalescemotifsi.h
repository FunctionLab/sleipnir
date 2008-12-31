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
#ifndef COALESCEMOTIFSI_H
#define COALESCEMOTIFSI_H

#include "fullmatrix.h"

namespace Sleipnir {

class CPST;

class CCoalesceMotifLibraryImpl {
protected:
	static const char	c_szBases[];
	static const char	c_szComplements[];
	static const size_t	c_iShift		= 2; // ceil( log2( strlen( c_szBases ) ) )
	static const char	c_cSeparator	= '|';

	enum EType {
		ETypeKMer,
		ETypeRC,
		ETypePST
	};

	static size_t CountKMers( size_t iK ) {

		return ( 1 << ( iK << 1 ) ); }

	static size_t CountRCs( size_t iK ) {
		size_t	iKMers;

		iKMers = CountKMers( iK ) >> 1;
		return ( iKMers - ( ( iK % 2 ) ? 0 : ( CountKMers( iK >> 1 ) >> 1 ) ) ); }

	static std::string GetComplement( const std::string& strKMer ) {
		std::string	strRet;
		size_t		i;

		strRet.resize( strKMer.length( ) );
		for( i = 0; i < strRet.length( ); ++i )
			strRet[ i ] = GetComplement( strKMer[ i ] );

		return strRet; }

	static char GetComplement( char cBase ) {
		const char*	pc;

		return ( ( pc = strchr( c_szBases, cBase ) ) ? c_szComplements[ pc - c_szBases ] : cBase ); }

	static uint32_t KMer2ID( const std::string& strKMer, bool fRC = false ) {
		size_t			i, iIndex;
		uint32_t		iRet;
		const char*		pc;
		unsigned char	c;

		for( i = iRet = 0; i < strKMer.length( ); ++i ) {
			iIndex = fRC ? ( strKMer.length( ) - i - 1 ) : i;
			if( !( pc = strchr( c_szBases, strKMer[ iIndex ] ) ) )
				return -1;
			c = (unsigned char)( pc - c_szBases );
			if( fRC ) {
				if( !( pc = strchr( c_szBases, c_szComplements[ c ] ) ) )
					return -1;
				c = (unsigned char)( pc - c_szBases ); }
			iRet = ( iRet << c_iShift ) | c; }

		return iRet; }

	static std::string ID2KMer( uint32_t iID, size_t iK ) {
		static const size_t	c_iMask	= ( 1 << c_iShift ) - 1;
		std::string	strRet;
		size_t		i;

		strRet.resize( iK );
		for( i = 0; i < iK; ++i ) {
			strRet[ iK - i - 1 ] = c_szBases[ iID & c_iMask ];
			iID >>= c_iShift; }

		return strRet; }

	static bool IsIgnorableKMer( const std::string& strKMer ) {
		size_t	i;

		for( i = 0; i < strKMer.size( ); ++i )
			if( !strchr( c_szBases, strKMer[ i ] ) )
				return true;

		return false; }

	static bool GetPWM( const std::string& strKMer, CFullMatrix<size_t>& MatPWM ) {
		size_t	i, j;

		if( ( MatPWM.GetColumns( ) != strlen( c_szBases ) ) || ( MatPWM.GetRows( ) != strKMer.length( ) ) ) {
			MatPWM.Initialize( strlen( c_szBases ), strKMer.length( ) );
			MatPWM.Clear( ); }
		for( i = 0; i < strKMer.length( ); ++i ) {
			for( j = 0; j < MatPWM.GetRows( ); ++j )
				if( strKMer[ i ] == c_szBases[ j ] )
					break;
			if( j >= MatPWM.GetRows( ) )
				return false;
			MatPWM.Get( j, i )++; }

		return true; }

	static std::string GetReverseComplement( const std::string& strKMer ) {
		std::string	strReverse;

		strReverse = strKMer;
		std::reverse( strReverse.begin( ), strReverse.end( ) );
		return GetComplement( strReverse ); }

	static float GetInformation( const CFullMatrix<size_t>& MatPWM ) {
		CDataMatrix			MatProbs;
		size_t				iPos, iFrom, iTo;
		std::vector<size_t>	veciTotals;
		float				d, dIC, dFrom, dTo;

		if( MatPWM.GetColumns( ) < 2 )
			return 0;

		veciTotals.resize( MatPWM.GetColumns( ) );
		for( iPos = 0; iPos < MatPWM.GetColumns( ); ++iPos )
			for( iFrom = 0; iFrom < MatPWM.GetRows( ); ++iFrom )
				veciTotals[ iPos ] += MatPWM.Get( iFrom, iPos );

		MatProbs.Initialize( MatPWM.GetRows( ), MatPWM.GetRows( ) );
		MatProbs.Clear( );
		for( iPos = 1; iPos < MatPWM.GetColumns( ); ++iPos )
			for( iFrom = 0; iFrom < MatPWM.GetRows( ); ++iFrom ) {
				if( !( dFrom = ( (float)MatPWM.Get( iFrom, iPos - 1 ) / veciTotals[ iPos - 1 ] ) ) )
					continue;
				for( iTo = 0; iTo < MatPWM.GetRows( ); ++iTo )
					if( dTo = ( (float)MatPWM.Get( iTo, iPos ) / veciTotals[ iPos ] ) )
						MatProbs.Get( iFrom, iTo ) += dFrom * dTo; }
		for( iFrom = 0; iFrom < MatProbs.GetRows( ); ++iFrom ) {
			for( d = 0,iTo = 0; iTo < MatProbs.GetColumns( ); ++iTo )
				d += MatProbs.Get( iFrom, iTo );
			if( !d )
				continue;
			for( iTo = 0; iTo < MatProbs.GetColumns( ); ++iTo )
				MatProbs.Get( iFrom, iTo ) /= d; }
		dIC = 0;
		for( iFrom = 0; iFrom < MatProbs.GetRows( ); ++iFrom )
			for( iTo = 0; iTo < MatProbs.GetColumns( ); ++iTo )
				if( d = MatProbs.Get( iFrom, iTo ) / ( MatPWM.GetColumns( ) - 1 ) )
					dIC += d * log( d );

		return ( dIC / ( -log( 2.0f ) * MatProbs.GetColumns( ) ) ); }

	CCoalesceMotifLibraryImpl( size_t iK ) : m_iK(iK), m_dPenaltyGap(1), m_dPenaltyMismatch(2) {
		uint32_t	i, j, iRC;

// TODO: if I was smart, I could do this with a direct encoding...
		m_veciKMer2RC.resize( GetKMers( ) );
		m_veciRC2KMer.resize( GetRCs( ) );
		std::fill( m_veciKMer2RC.begin( ), m_veciKMer2RC.end( ), -1 );
		for( i = j = 0; i < m_veciKMer2RC.size( ); ++i )
			if( m_veciKMer2RC[ i ] == -1 ) {
				iRC = KMer2ID( ID2KMer( i, m_iK ), true );
				if( iRC != i ) {
					m_veciKMer2RC[ i ] = m_veciKMer2RC[ iRC ] = j;
					m_veciRC2KMer[ j++ ] = i; } } }

	virtual ~CCoalesceMotifLibraryImpl( );

	std::string GetMotif( uint32_t ) const;
	CPST* CreatePST( uint32_t& );
	uint32_t MergeKMers( const std::string&, const std::string&, float );
	uint32_t MergeKMerRC( uint32_t, uint32_t, float );
	uint32_t MergeKMerPST( const std::string&, const CPST&, float );
	uint32_t MergeRCs( uint32_t, uint32_t, float );
	uint32_t MergeRCPST( uint32_t, const CPST&, float );
	uint32_t MergePSTs( const CPST&, const CPST&, float );
	float AlignKMers( const std::string&, const std::string&, float ) const;
	float AlignKMerRC( const std::string&, uint32_t, float ) const;
	float AlignKMerPST( const std::string&, const CPST&, float ) const;
	float AlignRCs( uint32_t, uint32_t, float ) const;
	float AlignRCPST( uint32_t, const CPST&, float ) const;
	float AlignPSTs( const CPST&, const CPST&, float ) const;
	uint32_t RemoveRCs( const CPST&, float, float );

	EType GetType( uint32_t iMotif ) const {

		if( iMotif < GetKMers( ) )
			return ETypeKMer;
		if( iMotif < GetBasePSTs( ) )
			return ETypeRC;

		return ETypePST; }

	uint32_t GetKMers( ) const {

		return (uint32_t)CountKMers( m_iK ); }

	uint32_t GetRCs( ) const {

			return (uint32_t)CountRCs( m_iK ); }

	uint32_t GetBaseRCs( ) const {

		return GetKMers( ); }

	uint32_t GetBasePSTs( ) const {

		return ( GetBaseRCs( ) + GetRCs( ) ); }

	const CPST* GetPST( uint32_t iMotif ) const {

		return m_vecpPSTs[ iMotif - GetBasePSTs( ) ]; }

	uint32_t GetPSTs( ) const {

		return (uint32_t)m_vecpPSTs.size( ); }

	float Align( const std::string& strFixed, const std::string& strMobile, float dCutoff, int& iRet ) const {
		const std::string&	strShort	= ( strFixed.length( ) < strMobile.length( ) ) ? strFixed : strMobile;
		const std::string&	strLong		= ( strFixed.length( ) < strMobile.length( ) ) ? strMobile : strFixed;
		size_t				iBegin, iEnd, iOffset, i, iLength, iMin;
		float				dRet, dCur;

		dRet = FLT_MAX;
		dCur = dCutoff / m_dPenaltyGap;
		iLength = strShort.length( ) + strLong.length( );
		iBegin = ( dCur < iLength ) ? (size_t)ceil( ( iLength - dCur ) / 2 ) : 0;
		iEnd = iLength + 1 - iBegin;
		for( iMin = 0,iOffset = iBegin; iOffset < iEnd; ++iOffset ) {
			i = ( iOffset <= strShort.length( ) ) ? iOffset : min( strShort.length( ), iLength - iOffset );
			dCur = m_dPenaltyGap * ( iLength - ( 2 * i ) );
			for( i = ( iOffset <= strShort.length( ) ) ? 0 : ( iOffset - strShort.length( ) );
				i < min( strShort.length( ), iOffset ); ++i )
				if( strShort[ i ] != strLong[ strLong.length( ) - iOffset + i ] )
					dCur += m_dPenaltyMismatch;
			if( dCur < dRet ) {
				dRet = dCur;
				iMin = iOffset; } }

		iRet = (int)strLong.length( ) - (int)iMin;
		if( strFixed.length( ) < strMobile.length( ) )
			iRet = -iRet;
		return dRet; }

	std::string GetRCOne( uint32_t iMotif ) const {

		return ID2KMer( (uint32_t)m_veciRC2KMer[ iMotif - GetBaseRCs( ) ], m_iK ); }

	float										m_dPenaltyGap;
	float										m_dPenaltyMismatch;
	size_t										m_iK;
	std::vector<uint32_t>						m_veciKMer2RC;
	std::vector<uint32_t>						m_veciRC2KMer;
	std::vector<CPST*>							m_vecpPSTs;
	std::set<std::pair<uint32_t, uint32_t> >	m_setpriiMerged;
};

}

#endif // COALESCEMOTIFSI_H
