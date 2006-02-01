#ifndef DATAPAIR_H
#define DATAPAIR_H

#include "datapairi.h"

namespace libBioUtils {

class CSlim;

class CDataPair : public CDataPairImpl {
public:
	bool Open( const char*, bool );
	bool Open( const CSlim& );
	size_t Quantify( float ) const;
	bool IsContinuous( ) const;
	unsigned char GetValues( ) const;
};

}

#endif // DATAPAIR_H
