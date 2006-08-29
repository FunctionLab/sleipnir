#include "stdafx.h"
#include "meta.h"

namespace libBioUtils {

const char	CMeta::c_szWS[]		= " \t\r\n";

void CMeta::Startup( int iVerbosity, size_t iRand ) {
	Category&			CatBioUtils	= Category::getInstance( c_szBioUtils );
	OstreamAppender*	pAppOstm	= new OstreamAppender( "cerr", &cerr );

	srand( ( iRand == -1 ) ?
#ifdef _MSC_VER
		GetTickCount( )
#else
		time( NULL )
#endif // _MSC_VER
		: iRand );
	pAppOstm->setLayout( new BasicLayout( ) );
	CatBioUtils.setAdditivity( false );
	CatBioUtils.setAppender( pAppOstm );
	CatBioUtils.setPriority( iVerbosity * Priority::ALERT ); }

void CMeta::Shutdown( ) {

	Category::shutdown( ); }

string CMeta::Filename( const string& str, char cRepl ) {
	size_t	i;
	string	strRet;
	char	c;

	for( i = 0; i < str.length( ); ++i )
		strRet += isalnum( c = str[ i ] ) ? c : cRepl;

	return strRet; }

void CMeta::Tokenize( const char* szString, vector<string>& vecstrTokens,
	const char* szDelim, bool fTrim ) {
	const char*	pc;
	string		strCur;
	bool		fPush;

	if( !( pc = szString ) )
		return;

	fPush = false;
	while( true ) {
		strCur.clear( );
		if( fTrim )
			for( ; *pc && strchr( szDelim, *pc ); ++pc );
		if( !*pc ) {
			if( !fTrim && fPush )
				vecstrTokens.push_back( strCur );
			return; }
		for( ; *pc && !strchr( szDelim, *pc ); ++pc )
			strCur += *pc;
		if( fPush = !!*pc )
			pc++;
		vecstrTokens.push_back( strCur ); } }

string CMeta::Basename( const char* szPath ) {
	const char*	pchOne;
	const char*	pchTwo;

	if( pchOne = strrchr( szPath, '\\' ) )
		pchOne++;
	if( pchTwo = strrchr( szPath, '/' ) )
		pchTwo++;

	return ( pchOne ? ( pchTwo ? max( pchOne, pchTwo ) : pchOne ) :
		( pchTwo ? pchTwo : szPath ) ); }

string CMeta::Trim( const char* szIn ) {
	size_t	iBeg, iEnd, iLen;

	if( !szIn || !( iLen = strlen( szIn ) ) )
		return "";

	for( iBeg = 0; szIn[ iBeg ]; ++iBeg )
		if( !isspace( szIn[ iBeg ] ) )
			break;
	for( iEnd = 0; szIn[ iLen - iEnd - 1 ]; ++iEnd )
		if( !isspace( szIn[ iLen - iEnd - 1 ] ) )
			break;

	return string( szIn + iBeg, iLen - iBeg - iEnd ); }

string CMeta::Deextension( const string& strName ) {
	size_t	i;

	return ( ( ( i = strName.rfind( c_cPeriod ) ) == string::npos ) ? strName :
		strName.substr( 0, i ) ); }

}
