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

#ifdef _MSC_VER
#define	FOR_EACH_DIRECTORY_FILE(strDir, strFile)											\
	HANDLE			hSearch;																\
	WIN32_FIND_DATA	sEntry;																	\
	bool			fOk;																	\
	for( fOk = ( ( hSearch = FindFirstFile( ( (strDir) + "/*" ).c_str( ), &sEntry ) )		\
		!= INVALID_HANDLE_VALUE ); fOk; fOk = !!FindNextFile( hSearch, &sEntry ) ) {		\
		if( sEntry.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )							\
			continue;																		\
		(strFile) = sEntry.cFileName;
#else // _MSC_VER
#include <dirent.h>
#include <sys/stat.h>

#define	FOR_EACH_DIRECTORY_FILE(strDir, strFile)											\
	DIR*			pDir;																	\
	struct dirent*	psEntry;																\
	struct stat		sStat;																	\
	pDir = opendir( (strDir).c_str( ) );													\
	while( pDir && ( psEntry = readdir( pDir ) ) ) {										\
		stat( ( (strDir) + '/' + psEntry->d_name ).c_str( ), &sStat );						\
		if( S_ISDIR( sStat.st_mode ) )														\
			continue;																		\
		(strFile) = psEntry->d_name;
#endif // _MSC_VER


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

	static size_t GetBufferSize( ) {

		return c_iBufferSize; }

	/*!
	 * \brief
	 * Return the next tab-delimited token from the given string.
	 * 
	 * \param szInput
	 * String from which the token is read.
	 * 
	 * \param pcEnd
	 * If non-null, outputs a pointer to the end of the token in the given string.
	 * 
	 * \returns
	 * String containing all characters up to (but excluding) the next tab or newline.
	 */
	static std::string OpenToken( const char* szInput, const char** ppcEnd = NULL ) {
		const char*	pcStart;
		const char*	pcEnd;
		char		c;

		do
			c = *(szInput++);
		while( ( c > 0 ) && ( c != '\t' ) && isspace( c ) );
		pcStart = szInput - 1;
		for( ; ( c > 0 ) && ( c != '\t' ) && !IsNewline( c ); c = *(szInput++) );
		pcEnd = szInput;
		if( ppcEnd )
			*ppcEnd = pcEnd - ( ( !c || IsNewline( c ) ) ? 1 : 0 );

		return std::string( pcStart, pcEnd - pcStart - 1 ); }
};

}

#endif // FILE_H
