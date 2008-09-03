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
#ifndef PSTI_H
#define PSTI_H

#include <sstream>
#include <string>
#include <vector>

#include "typesi.h"

namespace Sleipnir {

class CPSTImpl {
protected:
	struct SNode {
		SNode( ) : m_iTotal(0) { }

		SNode( unsigned char cCharacter ) : m_cCharacter(cCharacter), m_iCount(1), m_iTotal(0) { }

		unsigned char		m_cCharacter;
		uint16_t			m_iCount;
		uint16_t			m_iTotal;
		std::vector<SNode>	m_vecsChildren;
	};

	static std::string GetMotif( const SNode& sNode ) {
		std::ostringstream	ossm;
		size_t				i;

		switch( sNode.m_vecsChildren.size( ) ) {
			case 0:
				break;

			case 1:
				ossm << sNode.m_vecsChildren[ 0 ].m_cCharacter << GetMotif( sNode.m_vecsChildren[ 0 ] );
				break;

			default:
				for( i = 0; i < sNode.m_vecsChildren.size( ); ++i ) {
					const SNode&	sChild	= sNode.m_vecsChildren[ i ];

					ossm << ( i ? "|" : "" ) << '(' << sChild.m_iCount << ':' << sChild.m_cCharacter <<
						GetMotif( sChild ) << ')'; } }

		return ossm.str( ); }

	static size_t Add( unsigned char cCharacter, SNode& sNode ) {
		size_t	i;

		for( i = 0; i < sNode.m_vecsChildren.size( ); ++i )
			if( sNode.m_vecsChildren[ i ].m_cCharacter == cCharacter )
				break;
		if( i < sNode.m_vecsChildren.size( ) )
			sNode.m_vecsChildren[ i ].m_iCount++;
		else
			sNode.m_vecsChildren.push_back( SNode( cCharacter ) );
		sNode.m_iTotal++;

		return i; }

	static size_t Add( const std::string& str, SNode& sNode ) {

		return ( str.empty( ) ? 0 :
			( 1 + Add( str.substr( 1 ), sNode.m_vecsChildren[ Add( str[ 0 ], sNode ) ] ) ) ); }

	static size_t Add( const SNode& sIn, SNode& sOut ) {
		size_t	i, j, iCur, iRet;

		sOut.m_iCount += sIn.m_iCount;
		sOut.m_iTotal += sIn.m_iTotal;
		for( iRet = i = 0; i < sIn.m_vecsChildren.size( ); ++i ) {
			const SNode&	sChildIn	= sIn.m_vecsChildren[ i ];

			for( j = 0; j < sOut.m_vecsChildren.size( ); ++j )
				if( sChildIn.m_cCharacter == sOut.m_vecsChildren[ j ].m_cCharacter )
					break;
			if( j >= sOut.m_vecsChildren.size( ) ) {
				sOut.m_vecsChildren.push_back( SNode( sChildIn.m_cCharacter ) );
				sOut.m_vecsChildren.back( ).m_iCount = 0; }
			if( ( iCur = ( 1 + Add( sChildIn, sOut.m_vecsChildren[ j ] ) ) ) > iRet )
				iRet = iCur; }

		return iRet; }

	template<class tType>
	static size_t Add( const tType& Addition, size_t iOffset, SNode& sNode ) {
		size_t	i, iRet, iCur;

		if( !iOffset )
			return Add( Addition, sNode );

		iOffset--;
		iRet = 0;
		for( i = 0; i < sNode.m_vecsChildren.size( ); ++i ) {
			iCur = ( 1 + Add( Addition, iOffset, sNode.m_vecsChildren[ i ] ) );
			if( iCur > iRet )
				iRet = iCur; }

		return iRet; }

	CPSTImpl( size_t iArity ) : m_iDepth(0), m_iArity(iArity) { }

	long double GetMatch( const std::string& strTarget, const SNode& sNode, size_t iOffset,
		size_t& iMatched ) const {
		size_t		i, iCur, iMax;
		long double	dRet, dCur;

		if( strTarget.empty( ) || sNode.m_vecsChildren.empty( ) )
			return 1;

		if( iOffset ) {
			iOffset--;
			dRet = 0;
			for( iMax = i = 0; i < sNode.m_vecsChildren.size( ); ++i ) {
				iCur = 0;
				dCur = GetMatch( strTarget, sNode.m_vecsChildren[ i ], iOffset, iCur );
				if( dCur > dRet ) {
					dRet = dCur;
					iMax = iCur; } }
			iMatched = iMax;
			return dRet; }

		for( i = 0; i < sNode.m_vecsChildren.size( ); ++i ) {
			const SNode&	sChild	= sNode.m_vecsChildren[ i ];

			if( strTarget[ 0 ] == sChild.m_cCharacter ) {
				iMatched++;
				return ( sChild.m_iCount * GetMatch( strTarget.substr( 1 ), sChild, 0, iMatched ) /
					sNode.m_iTotal ); } }

		return 1; }

	template<class tType>
	float Align( const tType& Data, size_t iLength, float dPenaltyGap, float dPenaltyMismatch, float dCutoff,
		int& iOffset ) const {
		size_t	iBegin, iEnd, iCur, iMin;
		float	dRet, dCur;

		dRet = FLT_MAX;
		dCur = dCutoff / dPenaltyGap;
		iBegin = (size_t)ceil( iLength - dCur );
		iEnd = iLength + 1 - iBegin;
		for( iMin = 0,iCur = iBegin; iCur < iEnd; ++iCur ) {
			dCur = dPenaltyGap * ( ( iCur < iLength ) ? ( iLength - iCur ) : 0 );
			dCur += dPenaltyMismatch * CPSTImpl::Align( Data, iLength, iCur, m_sRoot ); 
			if( dCur < dRet ) {
				dRet = dCur;
				iMin = iCur; } }

		iOffset = (int)iLength - (int)iMin;
		return dRet; }

	size_t Align( const std::string& strSequence, size_t, size_t iOffset, const SNode& sNode ) const {
		size_t	i, iCur, iRet;

		if( strSequence.empty( ) || !iOffset )
			return 0;

		iRet = strSequence.length( );
		if( iOffset > iRet ) {
			iOffset--;
			for( i = 0; i < sNode.m_vecsChildren.size( ); ++i )
				if( ( iCur = Align( strSequence, 0, iOffset, sNode.m_vecsChildren[ i ] ) ) < iRet )
					iRet = iCur;
			return iRet; }

		for( i = 0; i < sNode.m_vecsChildren.size( ); ++i ) {
			const SNode&	sChild	= sNode.m_vecsChildren[ i ];

			iCur = ( ( sChild.m_cCharacter == strSequence[ strSequence.length( ) - iOffset ] ) ? 0 : 1 ) +
				Align( strSequence.substr( strSequence.length( ) - iOffset + 1 ), 0, iOffset - 1, sChild );
			if( iCur < iRet )
				iRet = iCur; }

		return iRet; }

	size_t Align( const SNode& sThem, size_t iDepth, size_t iOffset, const SNode& sUs ) const {
		size_t	i, j, iCur, iRet;

		if( sThem.m_vecsChildren.empty( ) )
			return 0;

		iRet = iDepth;
		if( iOffset > iRet ) {
			iOffset--;
			for( i = 0; i < sUs.m_vecsChildren.size( ); ++i )
				if( ( iCur = Align( sThem, iDepth, iOffset, sUs.m_vecsChildren[ i ] ) ) < iRet )
					iRet = iCur;
			return iRet; }
		if( iOffset < iRet ) {
			iDepth--;
			for( i = 0; i < sThem.m_vecsChildren.size( ); ++i )
				if( ( iCur = Align( sThem.m_vecsChildren[ i ], iDepth, iOffset, sUs ) ) < iRet )
					iRet = iCur;
			return iRet; }

		iDepth--;
		iOffset--;
		for( i = 0; i < sUs.m_vecsChildren.size( ); ++i ) {
			const SNode&	sChildUs	= sUs.m_vecsChildren[ i ];

			for( j = 0; j < sThem.m_vecsChildren.size( ); ++j ) {
				const SNode&	sChildThem	= sThem.m_vecsChildren[ j ];

				iCur = ( ( sChildUs.m_cCharacter == sChildThem.m_cCharacter ) ? 0 : 1 ) +
					Align( sChildThem, iDepth, iOffset, sChildUs );
				if( iCur < iRet )
					iRet = iCur; } }

		return iRet; }

	SNode	m_sRoot;
	size_t	m_iDepth;
	size_t	m_iArity;
};

}

#endif // PSTI_H
