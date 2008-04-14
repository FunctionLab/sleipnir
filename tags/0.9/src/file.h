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
#ifndef FILE_H
#define FILE_H

#include <iostream>
#include <string>

#include "filei.h"

namespace Sleipnir {

/*!
 * \brief
 * Parent class for types dealing with files, usually for text input.
 * 
 * \remarks
 * There's not often a reason to use CFile directly; it's mainly used internally by Sleipnir classes that
 * have to read tab-delimited text files.
 */
class CFile : protected CFileImpl {
public:
	static std::string OpenToken( std::istream& istm );
	static std::string OpenToken( const char* szInput, const char** pcEnd = NULL );
};

}

#endif // FILE_H
