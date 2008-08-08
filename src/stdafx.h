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
#ifndef STDAFX_H
#define STDAFX_H

#include <assert.h>
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
#define _mktemp_s			mktemp
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
#include <sstream>
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

	static void log4cpp( const char* szTag, const char* szFormat, va_list& valArgs ) {

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

#ifdef USE_LOG4CPP_STUB
extern Category		g_CatSleipnir;
#else // USE_LOG4CPP_STUB
extern const char	c_szSleipnir[];
extern Category&	g_CatSleipnir;
#endif // USE_LOG4CPP_STUB

}

#endif // STDAFX_H
