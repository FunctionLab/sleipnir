#include "stdafx.h"
#include "meta.h"

namespace libBioUtils {

const char	CMeta::c_szWS[]		= " \t\r\n";

void CMeta::Startup( int iVerbosity, size_t iRand ) {
// MEFIT OFF
	Category&			CatBioUtils	= Category::getInstance( c_szBioUtils );
	OstreamAppender*	pAppOstm	= new OstreamAppender( "cerr", &cerr );
// MEFIT ON

	srand( ( iRand == -1 ) ?
#ifdef _MSC_VER
		GetTickCount( )
#else
		time( NULL )
#endif // _MSC_VER
		: iRand );
// MEFIT OFF
	pAppOstm->setLayout( new BasicLayout( ) );
	CatBioUtils.setAdditivity( false );
	CatBioUtils.setAppender( pAppOstm );
	CatBioUtils.setPriority( iVerbosity * Priority::ALERT );
// MEFIT ON
}

void CMeta::Shutdown( ) {

// MEFIT OFF
	Category::shutdown( );
// MEFIT ON
}

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

bool CMeta::MapRead( unsigned char*& pbData, HANDLE& hndlMap, size_t& iData, const char* szFile ) {

	Unmap( pbData, hndlMap, iData );
#ifdef _MSC_VER
	HANDLE	hndlFile;

	if( !( hndlFile = CreateFile( szFile, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
		FILE_ATTRIBUTE_READONLY, NULL ) ) )
		return false;
	if( !( hndlMap = CreateFileMapping( hndlFile, NULL, PAGE_READONLY, 0,
		(DWORD)( iData = GetFileSize( hndlFile, NULL ) ), szFile ) ) ) {
		CloseHandle( hndlFile );
		return false; }
	CloseHandle( hndlFile );

	if( !( pbData = (unsigned char*)MapViewOfFile( hndlMap, FILE_MAP_READ, 0, 0, 0 ) ) ) {
		CloseHandle( hndlMap );
		return false; }
#else // _MSC_VER
	int			iFile;
	struct stat	sStat;

	if( !( iFile = open( szFile, O_RDONLY ) ) )
		return false;
	fstat( iFile, &sStat );
	iData = sStat.st_size;

	if( ( pbData = (unsigned char*)mmap( NULL, iData, PROT_READ, MAP_SHARED, iFile, 0 ) ) == MAP_FAILED ) {
		g_CatBioUtils.error( "CMeta::MapRead( %s ) %s", szFile, strerror( errno ) );
		pbData = NULL;
		close( iFile );
		return false; }
	close( iFile );
#endif // _MSC_VER

	return true; }

bool CMeta::MapWrite( unsigned char*& pbData, HANDLE& hndlMap, size_t iData, const char* szFile ) {

	Unmap( pbData, hndlMap, iData );
#ifdef _MSC_VER
	HANDLE	hndlFile;

	if( !( hndlFile = CreateFile( szFile, GENERIC_READ | GENERIC_WRITE, 0, NULL, CREATE_ALWAYS,
		FILE_ATTRIBUTE_NORMAL, NULL ) ) )
		return false;
	if( !( hndlMap = CreateFileMapping( hndlFile, NULL, PAGE_READWRITE, 0, (DWORD)iData, NULL ) ) ) {
		CloseHandle( hndlFile );
		return false; }
	CloseHandle( hndlFile );

	if( !( pbData = (unsigned char*)MapViewOfFile( hndlMap, FILE_MAP_WRITE, 0, 0, iData ) ) ) {
		CloseHandle( hndlMap );
		return false; }
#else // _MSC_VER
	int			iFile;
	struct stat	sStat;

	if( !( iFile = open( szFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP | S_IROTH ) ) )
		return false;
	lseek( iFile, iData - 1, SEEK_SET );
	write( iFile, &iData, 1 );
	if( ( pbData = (unsigned char*)mmap( NULL, iData, PROT_READ | PROT_WRITE, MAP_SHARED, iFile, 0 ) ) == MAP_FAILED ) {
		g_CatBioUtils.error( "CMeta::MapWrite( %s ) %s", szFile, strerror( errno ) );
		pbData = NULL;
		close( iFile );
		return false; }
	close( iFile );
#endif // _MSC_VER

	return true; }

void CMeta::Unmap( const unsigned char* pbData, HANDLE hndlMap, size_t iData ) {

#ifdef _MSC_VER
	if( pbData )
		UnmapViewOfFile( pbData );
	if( hndlMap )
		CloseHandle( hndlMap );
#else // _MSC_VER
	if( pbData )
		munmap( (void*)pbData, iData );
#endif // _MSC_VER
}

}
