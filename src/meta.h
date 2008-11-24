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
#ifndef META_H
#define META_H

#include <float.h>
#include <math.h>

#include <vector>

#include "metai.h"

#ifdef _MSC_VER
#include <windows.h>
#include <winsock.h>
#else // _MSC_VER
#include <stdarg.h>
#include <sys/time.h>

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

namespace Sleipnir {

/*!
 * \brief
 * Utility class containing static utility functions.
 * 
 * CMeta is critical in that it contains the Startup and Shutdown functions, which should be called at the
 * beginning and end of every process (usually in the main function) using Sleipnir.  These exist primarily
 * to set up and tear down logging, and can also be used to standardize the random seed for a process (useful
 * for testing).  Most other methods in CMeta are generic utilities for string manipulation and a few
 * operating system abstractions (particularly memory mapping).
 */
class CMeta : CMetaImpl {
public:
	/*!
	 * \brief
	 * String constant containing basic whitespace characters: space, tab, newline, return.
	 */
	static const char	c_szWS[];

	static std::string Filename( const std::string& strString, char cReplacement = '_' );
	static std::string Basename( const char* szPath );
	static void Tokenize( const char* szString, std::vector<std::string>& vecstrTokens,
		const char* szSeparators = "\t", bool fNoEmpties = false );
	static std::string Trim( const char* szString );
	static bool MapRead( unsigned char*& pbData, HANDLE& hndlMap, size_t& iSize, const char* szFile );
	static bool MapWrite( unsigned char*& pbData, HANDLE& hndlMap, size_t iSize, const char* szFile );
	static void Unmap( const unsigned char* pbData, HANDLE hndlMap, size_t iSize );

	/*!
	 * \brief
	 * Return true if the given value represents a missing value.
	 * 
	 * \param Value
	 * Value to test.
	 * 
	 * \returns
	 * True if the given value represents a missing value.
	 * 
	 * \remarks
	 * Returns true for either infinite or not-a-number (NaN) values, the latter of which is used as
	 * a standard missing value marker.  Templated to work with either doubles or floats, although the
	 * latter are more standard.  Used since == isn't reliable for NaN.
	 * 
	 * \see
	 * GetNaN
	 */
	template<class tType>
	static bool IsNaN( tType Value ) {

		return !_finite( Value ); }

	/*!
	 * \brief
	 * Return a standard missing value marker.
	 * 
	 * \returns
	 * Standard missing value marker.
	 * 
	 * \remarks
	 * Should be used anywhere a missing value is required, e.g. CPCL or CDat.
	 * 
	 * \see
	 * IsNaN
	 */
	static float GetNaN( ) {

		return (float)HUGE_VAL; }

	/*!
	 * \brief
	 * Given a filename, remove the file type extension (if any).
	 * 
	 * \param strName
	 * Filename to be de-extensioned.
	 * 
	 * \returns
	 * Filename without the trailing extension.
	 * 
	 * \remarks
	 * Actually removes anything after the last . in the given string.
	 */
	static std::string Deextension( const std::string& strName ) {
		size_t	i;

		return ( ( ( i = strName.rfind( c_cPeriod ) ) == std::string::npos ) ? strName :
			strName.substr( 0, i ) ); }

	/*!
	 * \brief
	 * Reorder a given item list based on a target ordering.
	 * 
	 * \param Items
	 * Iterator over items to be reordered.
	 * 
	 * \param veciOrder
	 * Indices at which each item should be placed.
	 * 
	 * Reorders a list of items based on a list of target indices.  For example, suppose the input list is
	 * [A, B, C] and the target order is [1, 0, 2].  Then after permutation, the vector of items will
	 * contain [B, A, C].  The reordering is done without copying more than one element at a time.
	 */
	template <class tIterator>
	static void Permute( tIterator Items, const std::vector<size_t>& veciOrder ) {
		std::vector<size_t>						veciMap, veciCopy;
		size_t									i;
		  typename iterator_traits<tIterator>::value_type	Temp;

		veciCopy.resize( veciOrder.size( ) );
		veciMap.resize( veciOrder.size( ) );
		for( i = 0; i < veciOrder.size( ); ++i ) {
			veciCopy[ i ] = veciOrder[ i ];
			veciMap[ veciOrder[ i ] ] = i; }

		for( i = 0; i < veciOrder.size( ); ++i ) {
			if( veciCopy[ i ] == i )
				continue;
			Temp = Items[ i ];
			Items[ i ] = Items[ veciCopy[ i ] ];
			Items[ veciCopy[ i ] ] = Temp;
			veciCopy[ veciMap[ i ] ] = veciCopy[ i ];
			veciMap[ veciCopy[ i ] ] = veciMap[ i ]; } }

	/*!
	 * \brief
	 * Reorder a given item list based on a target ordering.
	 * 
	 * \param vecItems
	 * Vector of items to be reordered.
	 * 
	 * \param veciOrder
	 * Indices at which each item should be placed.
	 * 
	 * Reorders a list of items based on a list of target indices.  For example, suppose the input list is
	 * [A, B, C] and the target order is [1, 0, 2].  Then after permutation, the vector of items will
	 * contain [B, A, C].  The reordering is done without copying more than one element at a time.
	 * 
	 * \remarks
	 * I'm not smart enough to get this to work seamlessly with both arrays and vectors, so this is my
	 * compromise.  Crazy STL and iterators...
	 */
	template <class tType>
	static void Permute( std::vector<tType>& vecItems, const std::vector<size_t>& veciOrder ) {

		Permute( vecItems.begin( ), veciOrder ); }

	/*!
	 * \brief
	 * Discretize a given continuous value based on a vector of bin edges.
	 * 
	 * \param Value
	 * Continuous value to be discretized.
	 * 
	 * \param vecQuants
	 * Bin edges used to discretize the given value.
	 * 
	 * \returns
	 * Discretized value based on the given bin edges.
	 * 
	 * Given N bin edges, the continuous value will be discretized into the range [0, N-1] depending on
	 * the first bin edge it is less than or equal to.  Note that this means that the last bin edge will be
	 * ignored.  Upper bin edges are inclusive, lower bin edges are exclusive.  This means that for bins
	 * [-0.1, 0.3, 0.6], the given value will be discretized into one of three outputs:
	 * - 0, corresponding to values less than or equal to -0.1.
	 * - 1, corresponding to values greater than -0.1 but less than or equal to 0.3.
	 * - 2, corresponding to values greater than 0.3.
	 * 
	 * \see
	 * CDataPair
	 */
	template <class tType>
	static size_t Quantize( tType Value, const std::vector<tType>& vecQuants ) {
		size_t	i;

		if( IsNaN( Value ) )
			return -1;

		for( i = 0; i < vecQuants.size( ); ++i )
			if( Value <= vecQuants[ i ] )
				break;

		return min( i, vecQuants.size( ) - 1 ); }

	/*!
	 * \brief
	 * Calculates the difference in microseconds between two timevals.
	 * 
	 * \param sBegin
	 * Earlier timeval.
	 * 
	 * \param sEnd
	 * Later timeval.
	 * 
	 * \returns
	 * Difference in microseconds between later and earlier timeval.
	 */
	static size_t GetMicroseconds( const struct timeval& sBegin, const struct timeval& sEnd ) {

		return ( ( sBegin.tv_sec == sEnd.tv_sec ) ? ( sEnd.tv_usec - sBegin.tv_usec ) :
			( ( 1000000 * ( sEnd.tv_sec - sBegin.tv_sec ) ) + sEnd.tv_usec - sBegin.tv_usec ) ); }

	/*!
	 * \brief
	 * Returns true if the given file path ends with the given extension.
	 * 
	 * \param strFile
	 * File path (or name).
	 * 
	 * \param strExtension
	 * Extension to test (including period, if desired).
	 * 
	 * \returns
	 * True if the given file path ends with the given extension.
	 */
	static bool IsExtension( const std::string& strFile, const std::string& strExtension ) {
		size_t	i;

		return ( ( strFile.length( ) >= strExtension.length( ) ) &&
			( ( i = strFile.rfind( strExtension ) ) ==
			( strFile.length( ) - strExtension.length( ) ) ) ); }

	/*!
	 * \brief
	 * Utility constructor that initializes Sleipnir (primarily log4cpp) at construction time and performs
	 * cleanup when destroyed.
	 * 
	 * \param iVerbosity
	 * If linked with log4cpp, the verbosity level for logging.
	 * 
	 * \param iRandomSeed
	 * Random seed for use with srand; if -1, the current time is used.
	 * 
	 * One (and only one) CMeta object should be created in a Sleipnir client's \c main function before
	 * making any library calls.  The object will be automatically destroyed as \c main exits, guaranteeing
	 * proper cleanup of Sleipnir (and log4cpp).
	 * 
	 * \remarks
	 * If Sleipnir is configured without log4cpp, iVerbosity will be ignored.
	 */
	CMeta( int iVerbosity, size_t iRandomSeed = 0 );
	/* 
	 * \brief
	 * Perform Sleipnir cleanup.
	 */
	~CMeta( );
};

}

#endif // META_H
