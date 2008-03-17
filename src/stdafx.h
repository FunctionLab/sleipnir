#ifndef STDAFX_H
#define STDAFX_H

#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#ifdef _MSC_VER
#include <io.h>
#include <winsock2.h>

#define socklen_t	int
#else // _MSC_VER
#include <dirent.h>
#include <errno.h>
#define __STDC_LIMIT_MACROS
#include <netinet/in.h>
#include <signal.h>
#include <stdarg.h>
#include <sys/mman.h>
#include <sys/socket.h>
#include <sys/time.h>

#define _lseek				lseek
#define _read				read
#define _write				write
#define closesocket			close
#define SOCKET				int
#define strcpy_s(a,b,c)		strcpy(a,c)
#define strncpy_s(a,b,c,d)	strncpy(a,c,d)

inline int sprintf_s( char* szDest, size_t iSize, const char* szFormat,
	... ) {
	va_list	valArgs;

	va_start( valArgs, szFormat );
	return vsprintf( szDest, szFormat, valArgs ); }

inline size_t max( size_t iOne, size_t iTwo ) {

	return ( ( iOne > iTwo ) ? iOne : iTwo ); }
#endif // _MSC_VER

#pragma warning (disable: 4267)
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <strstream>
#include <vector>
using namespace std;

#ifndef USE_LOG4CPP_STUB
#include <log4cpp/Category.hh>
#include <log4cpp/OstreamAppender.hh>
using namespace log4cpp;
#endif // USE_LOG4CPP_STUB

#if _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else // _MSC_VER
#define ios_base	ios
#define UINT		unsigned int
#endif // _MSC_VER

namespace pnl { }
using namespace pnl;

namespace Sleipnir {

#ifdef USE_LOG4CPP_STUB

struct Category {

	static void shutdown( ) { }

	static void log4cpp( const char* szTag, const char* szFormat, const va_list& valArgs ) {

		fprintf( stderr, "%d ", time( NULL ) );
		fprintf( stderr, szTag );
		fprintf( stderr, " : " );
		vfprintf( stderr, szFormat, valArgs );
		fprintf( stderr, "\n" ); }

	void error( const char* szFormat, ... ) const {
		va_list	valArgs;

		va_start( valArgs, szFormat );
		log4cpp( "ERROR", szFormat, valArgs ); }

	void info( const char* szFormat, ... ) const {
		va_list	valArgs;

		va_start( valArgs, szFormat );
		log4cpp( "INFO", szFormat, valArgs ); }

	void notice( const char* szFormat, ... ) const {
		va_list	valArgs;

		va_start( valArgs, szFormat );
		log4cpp( "NOTICE", szFormat, valArgs ); }

	void warn( const char* szFormat, ... ) const {
		va_list	valArgs;

		va_start( valArgs, szFormat );
		log4cpp( "WARN", szFormat, valArgs ); }

	void debug( const char* szFormat, ... ) const {
		va_list	valArgs;

		va_start( valArgs, szFormat );
		log4cpp( "DEBUG", szFormat, valArgs ); }
};

#endif // USE_LOG4CPP_STUB

extern const char	c_szSleipnir[];
#ifdef USE_LOG4CPP_STUB
extern Category		g_CatSleipnir;
#else // USE_LOG4CPP_STUB
extern Category&	g_CatSleipnir;
#endif // USE_LOG4CPP_STUB

}

#endif // STDAFX_H
