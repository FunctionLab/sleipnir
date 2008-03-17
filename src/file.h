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
