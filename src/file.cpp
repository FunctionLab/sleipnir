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
string CFile::OpenToken( istream& istm ) {
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

	while( isspace( c = *(szInput++) ) && ( c != '\t' ) && c );
	for( ; c && ( c != -1 ) && ( c != '\t' ) && !IsNewline( c ); c = *(szInput++) )
		strRet += c;
	if( pcEnd )
		*pcEnd = szInput - ( ( !c || IsNewline( c ) ) ? 1 : 0 );

	return strRet; }

}
