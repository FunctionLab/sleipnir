#ifndef META_H
#define META_H

#include <float.h>
#include <math.h>

#include <vector>

#include "metai.h"

#ifndef _MSC_VER
#define _finite	isfinite
#endif // _MSC_VER

namespace libBioUtils {

class CMeta : CMetaImpl {
public:
	static const char	c_szWS[];

	static void Startup( int, size_t = 0 );
	static void Shutdown( );
	static float GetNaN( );
	static std::string Filename( const std::string&, char = '_' );
	static std::string Basename( const char* );
	static void Tokenize( const char*, std::vector<std::string>&, const char* = "\t",
		bool = false );

	static bool IsNaN( float d ) {

		return !_finite( d ); }

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
};

}

#endif // META_H
