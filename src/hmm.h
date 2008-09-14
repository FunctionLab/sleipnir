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
#ifndef HMM_H
#define HMM_H

#include "hmmi.h"

namespace Sleipnir {

class CHMM : protected CHMMImpl {
public:
	void Open( size_t iDegree, const std::string& strAlphabet ) {

		m_iDegree = iDegree;
		m_strAlphabet = strAlphabet;
		m_MatTransitions.Initialize( GetStates( ), GetSymbols( ) - 1 );
		m_MatTransitions.Clear( ); }

	void Save( std::ostream& ostm ) const {
		size_t	i, j;

		ostm << m_iDegree << endl;
		ostm << m_strAlphabet << endl;
		for( i = 0; i < m_MatTransitions.GetRows( ); ++i ) {
			ostm << Decode( i );
			for( j = 0; j < m_MatTransitions.GetColumns( ); ++j )
				ostm << '\t' << m_MatTransitions.Get( i, j );
			ostm << endl; } }

	bool Add( const std::string& strData ) {
		size_t	i, iState;

		iState = 0;
		for( i = 0; i < strData.length( ); ++i ) {
			m_MatTransitions.Get( iState, Encode( strData[ i ] ) )++;
			iState = ( i < m_iDegree ) ? Encode( strData, i + 1 ) :
				Encode( strData.substr( i - m_iDegree + 1, m_iDegree ), m_iDegree );
			if( iState == -1 )
				return false; }

		return true; }

	std::string Get( size_t iLength ) const {
		std::string	strRet;
		size_t		i, iState, iTotal, iCur;

		for( iState = 0; strRet.length( ) < iLength; ) {
			for( iTotal = i = 0; i < m_MatTransitions.GetColumns( ); ++i )
				iTotal += m_MatTransitions.Get( iState, i );
			iCur = (size_t)( ( (float)rand( ) / ( RAND_MAX + 1 ) ) * iTotal );
			for( i = 0; ( i + 1 ) < m_MatTransitions.GetColumns( ); ++i ) {
				if( iCur < m_MatTransitions.Get( iState, i ) )
					break;
				iCur -= m_MatTransitions.Get( iState, i ); }
			strRet += m_strAlphabet[ i ];
			iState = ( ( iState * GetSymbols( ) ) + i + 1 ) % m_MatTransitions.GetRows( ); }

		return strRet; }

	void SetUniform( ) {
		size_t	i, j;

		for( i = 0; i < m_MatTransitions.GetRows( ); ++i )
			for( j = 0; j < m_MatTransitions.GetColumns( ); ++j )
				m_MatTransitions.Set( i, j, 1 ); }
};

}

#endif // HMM_H
