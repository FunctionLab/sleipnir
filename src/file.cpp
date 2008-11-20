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
#include "stdafx.h"
#include "file.h"

namespace Sleipnir {

bool CFileImpl::IsNewline( char c ) {

	return ( ( c == '\n' ) || ( c == '\r' ) ); }

/*!
 * \brief
 * Return the next tab-delimited token from the given input stream.
 * 
 * \param istm
 * Input stream from which the token is read.
 * 
 * \returns
 * String containing all characters up to (but excluding) the next tab or newline.
 */
string CFile::OpenToken( std::istream& istm ) {
	string	strRet;
	char	c;

	while( isspace( c = istm.get( ) ) && ( c != '\t' ) && ( c != EOF ) );
	for( ; c && ( c != -1 ) && ( c != '\t' ) && !IsNewline( c ); c = istm.get( ) )
		strRet += c;
	if( IsNewline( c ) )
		istm.unget( );

	return strRet; }

/*!
 * \brief
 * Return the next tab-delimited token from the given string.
 * 
 * \param szInput
 * String from which the token is read.
 * 
 * \param pcEnd
 * If non-null, outputs a pointer to the end of the token in the given string.
 * 
 * \returns
 * String containing all characters up to (but excluding) the next tab or newline.
 */
string CFile::OpenToken( const char* szInput, const char** pcEnd ) {
	string	strRet;
	char	c;

	do
		c = *(szInput++);
	while( c && ( c >= 0 ) && ( c < 256 ) && ( c != '\t' ) && isspace( c ) );
	for( ; c && ( c != -1 ) && ( c != '\t' ) && !IsNewline( c ); c = *(szInput++) )
		strRet += c;
	if( pcEnd )
		*pcEnd = szInput - ( ( !c || IsNewline( c ) ) ? 1 : 0 );

	return strRet; }

}
