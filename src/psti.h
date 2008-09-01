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

	static size_t Add( const std::string& strOne, const std::string& strTwo, size_t iOffset, SNode& sNode ) {
		size_t	i, j;

		if( !iOffset ) {
			i = Add( strOne, sNode );
			j = Add( strTwo, sNode );
			return max( i, j ); }

		i = Add( strOne[ 0 ], sNode );
		return ( 1 + Add( strOne.substr( 1 ), strTwo, iOffset - 1, sNode.m_vecsChildren[ i ] ) ); }

	CPSTImpl( size_t iArity ) : m_iDepth(0), m_iArity(iArity) { }

	long double GetMatch( const std::string& strTarget, const SNode& sNode ) const {
		size_t	i;

		if( strTarget.empty( ) || sNode.m_vecsChildren.empty( ) )
			return 1;

		for( i = 0; i < sNode.m_vecsChildren.size( ); ++i ) {
			const SNode&	sChild	= sNode.m_vecsChildren[ i ];

			if( strTarget[ 0 ] == sChild.m_cCharacter )
				return ( sChild.m_iCount * GetMatch( strTarget.substr( 1 ), sChild ) / sNode.m_iTotal ); }

		return 0; }

	SNode	m_sRoot;
	size_t	m_iDepth;
	size_t	m_iArity;
};

}

#endif // PSTI_H
