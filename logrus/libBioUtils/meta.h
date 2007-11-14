#ifndef META_H
#define META_H

#include <float.h>
#include <math.h>

#include <vector>

#include "metai.h"

#ifdef _MSC_VER
#include <windows.h>
#else // _MSC_VER
#include <stdarg.h>

#define _finite	isfinite

typedef size_t	HANDLE;

inline int sprintf_s( char* szDest, const char* szFormat, ... ) {
    va_list valArgs;

    va_start( valArgs, szFormat );
    return vsprintf( szDest, szFormat, valArgs ); }
#endif // _MSC_VER

#ifndef ARRAYSIZE
#define ARRAYSIZE(a)	(sizeof(a)/sizeof(*a))
#endif // ARRAYSIZE

namespace libBioUtils {

class CMeta : CMetaImpl {
public:
	static const char	c_szWS[];

	static void Startup( int, size_t = 0 );
	static void Shutdown( );
	static std::string Filename( const std::string&, char = '_' );
	static std::string Basename( const char* );
	static std::string Deextension( const std::string& );
	static void Tokenize( const char*, std::vector<std::string>&, const char* = "\t",
		bool = false );
	static std::string Trim( const char* );
	static bool MapRead( unsigned char*&, HANDLE&, size_t&, const char* );
	static bool MapWrite( unsigned char*&, HANDLE&, size_t, const char* );
	static void Unmap( const unsigned char*, HANDLE, size_t );

	template<class tType>
	static bool IsNaN( tType ) {

		return false; }

	static bool IsNaN( float d ) {

		return !_finite( d ); }

	static float GetNaN( ) {

		return (float)HUGE_VAL; }

	template <class tType>
	static void Permute( std::vector<tType>& vecItems,
		const std::vector<size_t>& veciOrder ) {
		std::vector<size_t>	veciMap, veciCopy;
		size_t				i;
		tType				Temp;

		veciCopy.resize( veciOrder.size( ) );
		veciMap.resize( veciOrder.size( ) );
		for( i = 0; i < veciOrder.size( ); ++i ) {
			veciCopy[ i ] = veciOrder[ i ];
			veciMap[ veciOrder[ i ] ] = i; }

		for( i = 0; i < vecItems.size( ); ++i ) {
			if( veciCopy[ i ] == i )
				continue;
			Temp = vecItems[ i ];
			vecItems[ i ] = vecItems[ veciCopy[ i ] ];
			vecItems[ veciCopy[ i ] ] = Temp;
			veciCopy[ veciMap[ i ] ] = veciCopy[ i ];
			veciMap[ veciCopy[ i ] ] = veciMap[ i ]; } }

	template <class tType>
	static void Permute( tType* const aItems, size_t iItems,
		const std::vector<size_t>& veciOrder ) {
		std::vector<size_t>	veciMap, veciCopy;
		size_t				i;
		tType				Temp;

		veciCopy.resize( veciOrder.size( ) );
		veciMap.resize( veciOrder.size( ) );
		for( i = 0; i < veciOrder.size( ); ++i ) {
			veciCopy[ i ] = veciOrder[ i ];
			veciMap[ veciOrder[ i ] ] = i; }

		for( i = 0; i < iItems; ++i ) {
			if( veciCopy[ i ] == i )
				continue;
			Temp = aItems[ i ];
			aItems[ i ] = aItems[ veciCopy[ i ] ];
			aItems[ veciCopy[ i ] ] = Temp;
			veciCopy[ veciMap[ i ] ] = veciCopy[ i ];
			veciMap[ veciCopy[ i ] ] = veciMap[ i ]; } }

	template <class tType>
	static size_t Quantize( tType Value, const std::vector<tType>& vecQuants ) {
		size_t	i;

		if( IsNaN( Value ) )
			return -1;

		for( i = 0; i < vecQuants.size( ); ++i )
			if( Value <= vecQuants[ i ] )
				break;

		return min( i, vecQuants.size( ) - 1 ); }
};

}

#endif // META_H
