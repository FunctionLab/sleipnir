#include "stdafx.h"
#include "file.h"

namespace libBioUtils {

bool CFileImpl::IsNewline( char c ) {

	return ( ( c == '\n' ) || ( c == '\r' ) ); }

string CFile::OpenToken( istream& istmInput ) {
	string	strRet;
	char	c;

	while( isspace( c = istmInput.get( ) ) && ( c != '\t' ) && ( c != EOF ) );
	for( ; c && ( c != -1 ) && ( c != '\t' ) && !IsNewline( c ); c = istmInput.get( ) )
		strRet += c;
	if( IsNewline( c ) )
		istmInput.unget( );

	return strRet; }

}
