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
#ifndef HMMI_H
#define HMMI_H

#include "fullmatrix.h"

namespace Sleipnir {

class CHMMImpl {
protected:
	size_t GetStates( ) const {
		size_t	i, iRet;

		iRet = 1;
		for( i = 0; i < m_iDegree; ++i )
			iRet *= GetSymbols( );

		return iRet; }

	size_t GetSymbols( ) const {

		return ( m_strAlphabet.length( ) + 1 ); }

	size_t Encode( const std::string& strData, size_t iCount ) const {
		size_t	i, iRet;

		iRet = 0;
		for( i = 0; ( i < iCount ) && ( i < strData.size( ) ); ++i )
			iRet = ( iRet * GetSymbols( ) ) + Encode( strData[ i ] ) + 1;

		return iRet; }

	std::string Decode( size_t iState ) const {
		std::string	strRet;
		size_t		i, iCur;

		for( i = 0; i < m_iDegree; ++i ) {
			iCur = iState % GetSymbols( );
			strRet = ( iCur ? m_strAlphabet[ iCur - 1 ] : '_' ) + strRet;
			iState = ( iState - iCur ) / GetSymbols( ); }

		return strRet; }

	size_t Encode( char cDatum ) const {
		size_t	i;

		for( i = 0; i < m_strAlphabet.length( ); ++i )
			if( cDatum == m_strAlphabet[ i ] )
				return i;

// I can't think of a better way to handle all exception cases...
		return ( rand( ) % m_strAlphabet.length( ) ); }

	std::string			m_strAlphabet;
	size_t				m_iDegree;
	CFullMatrix<size_t>	m_MatTransitions;
};

}

#endif // HMMI_H
