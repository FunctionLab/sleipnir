#ifndef FILEI_H
#define FILEI_H

namespace libBioUtils {

class CFileImpl {
protected:
	static const size_t c_iBufferSize	= 131072;

	static bool IsNewline( char );
};

}

#endif // FILEI_H
