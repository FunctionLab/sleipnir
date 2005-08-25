#ifndef FILE_H
#define FILE_H

#include <iostream>
#include <string>

#include "filei.h"

namespace libBioUtils {

class CFile : protected CFileImpl {
public:
	static std::string OpenToken( std::istream& );
};

}

#endif // FILE_H
