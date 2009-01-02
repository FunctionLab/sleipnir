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
#ifndef FILEI_H
#define FILEI_H

#include "typesi.h"

namespace Sleipnir {

class CFileImpl {
protected:
	static const size_t c_iBufferSize	= 2097152;

	static void SaveString( std::ostream& ostm, const std::string& str ) {
		uint32_t	iLength;

		iLength = (uint32_t)str.length( );
		ostm.write( (const char*)&iLength, sizeof(iLength) );
		ostm.write( str.c_str( ), iLength ); }

	static void OpenString( std::istream& istm, std::string& str ) {
		uint32_t	iLength;

		istm.read( (char*)&iLength, sizeof(iLength) );
		str.resize( iLength );
		istm.read( &str[ 0 ], iLength ); }

	static bool IsNewline( char c ) {

		return ( ( c == '\n' ) || ( c == '\r' ) ); }
};

}

#endif // FILEI_H
