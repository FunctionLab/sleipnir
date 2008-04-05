#include "stdafx.h"
#include "meta.h"

namespace Sleipnir {

const char	CMeta::c_szWS[]		= " \t\r\n";

/*!
 * \brief
 * Perform Sleipnir initialization.
 * 
 * \param iVerbosity
 * If linked with log4cpp, the verbosity level for logging.
 * 
 * \param iRandomSeed
 * Random seed for use with srand; if -1, the current time is used.
 * 
 * A Startup object should be created at the beginning of any process using Sleipnir, usually at the top of
 * main.  This serves primarily to initialize logging functions, but if you're doing a test that you want to
 * behave consistently, you can also seed the random number generator with a constant value.
 * 
 * \remarks
 * If Sleipnir is configured without log4cpp, iVerbosity will be ignored.
 */
CMeta::Startup::Startup( int iVerbosity, size_t iRandomSeed ) {
#ifndef USE_LOG4CPP_STUB
	Category&			CatSleipnir	= Category::getInstance( c_szSleipnir );
	OstreamAppender*	pAppOstm	= new OstreamAppender( "cerr", &cerr );
#endif // USE_LOG4CPP_STUB

	srand( ( iRandomSeed == -1 ) ?
#ifdef _MSC_VER
		GetTickCount( )
#else
		time( NULL )
#endif // _MSC_VER
		: iRandomSeed );
#ifndef USE_LOG4CPP_STUB
	pAppOstm->setLayout( new BasicLayout( ) );
	CatSleipnir.setAdditivity( false );
	CatSleipnir.setAppender( pAppOstm );
	CatSleipnir.setPriority( iVerbosity * Priority::ALERT );
#endif // USE_LOG4CPP_STUB
}

/*!
 * \brief
 * Perform Sleipnir cleanup.
 */
CMeta::Startup::~Startup( ) {

	Category::shutdown( ); }

/*!
 * \brief
 * Replace all non-alphanumeric characters in a string with the given replacement.
 * 
 * \param strString
 * String in which non-alphanumeric characters are replaced.
 * 
 * \param cReplacement
 * Character used to replace non-alphanumeric characters.
 * 
 * \returns
 * String with non-alphanumeric characters replaced.
 * 
 * This method is intended to clean a string to make it appropriate for use as a file name or other
 * alphanumeric identifier; given a string, non-alphanumeric characters are replaced with a configurable
 * character, usually underscore.
 */
string CMeta::Filename( const std::string& strString, char cReplacement ) {
	size_t	i;
	string	strRet;
	char	c;

	for( i = 0; i < strString.length( ); ++i )
		strRet += isalnum( c = strString[ i ] ) ? c : cReplacement;

	return strRet; }

/*!
 * \brief
 * Tokenize a given string based on one or more delimiter characters.
 * 
 * \param szString
 * String to be tokenized.
 * 
 * \param vecstrTokens
 * Output vector of tokens from the given string.
 * 
 * \param szSeparators
 * One or more separator characters used to split tokens.
 * 
 * \param fNoEmpties
 * If true, discard empty strings between delimiters; otherwise, include them in the output.
 */
void CMeta::Tokenize( const char* szString, std::vector<std::string>& vecstrTokens, const char* szSeparators,
	bool fNoEmpties ) {
	const char*	pc;
	string		strCur;
	bool		fPush;

	if( !( pc = szString ) )
		return;

	fPush = false;
	while( true ) {
		strCur.clear( );
		if( fNoEmpties )
			for( ; *pc && strchr( szSeparators, *pc ); ++pc );
		if( !*pc ) {
			if( !fNoEmpties && fPush )
				vecstrTokens.push_back( strCur );
			return; }
		for( ; *pc && !strchr( szSeparators, *pc ); ++pc )
			strCur += *pc;
		if( fPush = !!*pc )
			pc++;
		vecstrTokens.push_back( strCur ); } }

/*!
 * \brief
 * Attempt to return the filename portion of a path in a platform-independent manner.
 * 
 * \param szPath
 * File path from which filename is extracted.
 * 
 * \returns
 * Filename portion of the given path.
 * 
 * \remarks
 * Actually looks for the last / or \ character in the string and returns everything to the right of that.
 * Of course this won't always work, but it tends to do an awfully good job.
 */
string CMeta::Basename( const char* szPath ) {
	const char*	pchOne;
	const char*	pchTwo;

	if( pchOne = strrchr( szPath, '\\' ) )
		pchOne++;
	if( pchTwo = strrchr( szPath, '/' ) )
		pchTwo++;

	return ( pchOne ? ( pchTwo ? max( pchOne, pchTwo ) : pchOne ) :
		( pchTwo ? pchTwo : szPath ) ); }

/*!
 * \brief
 * Trim whitespace from the beginning and end of the given string.
 * 
 * \param szString
 * String from which whitespace is trimmed.
 * 
 * \returns
 * String with whitespace removed.
 */
string CMeta::Trim( const char* szString ) {
	size_t	iBeg, iEnd, iLen;

	if( !szString || !( iLen = strlen( szString ) ) )
		return "";

	for( iBeg = 0; szString[ iBeg ]; ++iBeg )
		if( !isspace( szString[ iBeg ] ) )
			break;
	for( iEnd = 0; szString[ iLen - iEnd - 1 ]; ++iEnd )
		if( !isspace( szString[ iLen - iEnd - 1 ] ) )
			break;

	return string( szString + iBeg, iLen - iBeg - iEnd ); }

/*!
 * \brief
 * Memory map an existing file read-only in a largely platform-independent manner.
 * 
 * \param pbData
 * Output pointer to mapped file data.
 * 
 * \param hndlMap
 * Output handle to file map; ignored on non-Windows platforms.
 * 
 * \param iSize
 * Output size of mapped file.
 * 
 * \param szFile
 * File name to map.
 * 
 * \returns
 * True if file was memory mapped successfully.
 * 
 * \remarks
 * This has been successfully tested on Windows, Linux, and (minimally) Mac OS.  It plays some mildly
 * ugly tricks to provide a standard memory mapping interface on all three systems, but it should work;
 * your mileage may vary.  On success, pbData will be of size iSize.  An opened map should be closed with
 * Unmap.
 * 
 * \see
 * MapWrite
 */
bool CMeta::MapRead( unsigned char*& pbData, HANDLE& hndlMap, size_t& iSize, const char* szFile ) {

	Unmap( pbData, hndlMap, iSize );
#ifdef _MSC_VER
	HANDLE	hndlFile;

	if( !( hndlFile = CreateFile( szFile, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
		FILE_ATTRIBUTE_READONLY, NULL ) ) )
		return false;
	if( !( hndlMap = CreateFileMapping( hndlFile, NULL, PAGE_READONLY, 0,
		(DWORD)( iSize = GetFileSize( hndlFile, NULL ) ), szFile ) ) ) {
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
	iSize = sStat.st_size;

	if( ( pbData = (unsigned char*)mmap( NULL, iSize, PROT_READ, MAP_SHARED, iFile, 0 ) ) == MAP_FAILED ) {
		g_CatSleipnir.error( "CMeta::MapRead( %s ) %s", szFile, strerror( errno ) );
		pbData = NULL;
		close( iFile );
		return false; }
	close( iFile );
#endif // _MSC_VER

	return true; }

/*!
 * \brief
 * Create a new writeable memory mapped file in a largely platform-independent manner.
 * 
 * \param pbData
 * Output pointer to mapped file data.
 * 
 * \param hndlMap
 * Output handle to file map; ignored on non-Windows platforms.
 * 
 * \param iSize
 * Size of desired memory map.
 * 
 * \param szFile
 * File name to map.
 * 
 * \returns
 * True if file was memory mapped successfully.
 * 
 * This function creates a new file, sizes it to the requested number of bytes, and memory maps it
 * writeably.  Using it on an existing file will generally destroy it and overwrite it with new data.
 * 
 * \remarks
 * This has been successfully tested on Windows, Linux, and (minimally) Mac OS.  It plays some mildly
 * ugly tricks to provide a standard memory mapping interface on all three systems, but it should work;
 * your mileage may vary.  On success, pbData will be of size iSize.  An opened map should be closed with
 * Unmap.
 * 
 * \see
 * MapRead
 */
bool CMeta::MapWrite( unsigned char*& pbData, HANDLE& hndlMap, size_t iSize, const char* szFile ) {

	Unmap( pbData, hndlMap, iSize );
#ifdef _MSC_VER
	HANDLE	hndlFile;

	if( !( hndlFile = CreateFile( szFile, GENERIC_READ | GENERIC_WRITE, 0, NULL, CREATE_ALWAYS,
		FILE_ATTRIBUTE_NORMAL, NULL ) ) )
		return false;
	if( !( hndlMap = CreateFileMapping( hndlFile, NULL, PAGE_READWRITE, 0, (DWORD)iSize, NULL ) ) ) {
		CloseHandle( hndlFile );
		return false; }
	CloseHandle( hndlFile );

	if( !( pbData = (unsigned char*)MapViewOfFile( hndlMap, FILE_MAP_WRITE, 0, 0, iSize ) ) ) {
		CloseHandle( hndlMap );
		return false; }
#else // _MSC_VER
	int			iFile;
	struct stat	sStat;

	if( !( iFile = open( szFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP | S_IROTH ) ) )
		return false;
	lseek( iFile, iSize - 1, SEEK_SET );
	write( iFile, &iSize, 1 );
	if( ( pbData = (unsigned char*)mmap( NULL, iSize, PROT_READ | PROT_WRITE, MAP_SHARED, iFile, 0 ) ) == MAP_FAILED ) {
		g_CatSleipnir.error( "CMeta::MapWrite( %s ) %s", szFile, strerror( errno ) );
		pbData = NULL;
		close( iFile );
		return false; }
	close( iFile );
#endif // _MSC_VER

	return true; }

/*!
 * \brief
 * Unmap a memory map in a largely platform-independent manner.
 * 
 * \param pbData
 * Pointer to mapped file data.
 * 
 * \param hndlMap
 * Handle to file map; ignored on non-Windows platforms.
 * 
 * \param iSize
 * Size of memory map.
 * 
 * \remarks
 * Should be used to close any maps opened with MapRead or MapWrite.
 */
void CMeta::Unmap( const unsigned char* pbData, HANDLE hndlMap, size_t iSize ) {

#ifdef _MSC_VER
	if( pbData )
		UnmapViewOfFile( pbData );
	if( hndlMap )
		CloseHandle( hndlMap );
#else // _MSC_VER
	if( pbData )
		munmap( (void*)pbData, iSize );
#endif // _MSC_VER
}

}
