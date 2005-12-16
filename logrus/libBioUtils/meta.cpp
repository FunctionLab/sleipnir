#include "stdafx.h"
#include "meta.h"

namespace libBioUtils {

const char	CMeta::c_szWS[]		= " \t\r\n";
const float	CMetaImpl::c_dNaN	= (float)HUGE_VAL;

void CMeta::Startup( int iVerbosity, size_t iRand ) {
	Category&			CatBioUtils	= Category::getInstance( c_szBioUtils );
	OstreamAppender*	pAppOstm	= new OstreamAppender( "cerr", &cerr );

	srand( ( iRand == -1 ) ?
#ifdef WIN32
		GetTickCount( )
#else
		time( NULL )
#endif // WIN32
		: iRand );
	pAppOstm->setLayout( new BasicLayout( ) );
	CatBioUtils.setAdditivity( false );
	CatBioUtils.setAppender( pAppOstm );
	CatBioUtils.setPriority( iVerbosity * Priority::ALERT ); }

void CMeta::Shutdown( ) {

	Category::shutdown( ); }

float CMeta::GetNaN( ) {

	return c_dNaN; }

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

}
