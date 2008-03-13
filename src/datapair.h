#ifndef DATAPAIR_H
#define DATAPAIR_H

#include "datapairi.h"

namespace Sleipnir {

class CSlim;

class CDataPair : public CDataPairImpl {
public:
	bool Open( const char*, bool, bool = false, size_t = 2, bool = false );
	bool Open( const CSlim& );
	bool OpenQuants( const char* );
	void SetQuants( const float*, size_t );
	void SetQuants( const std::vector<float>& );
	size_t Quantize( float ) const;
	bool IsContinuous( ) const;
	unsigned char GetValues( ) const;

	bool Open( const CDat& Dat, const std::vector<CGenes*>& vecpOther,
		const CGenome& Genome, bool fPositives ) {

		return CDat::Open( Dat, vecpOther, Genome, fPositives ); }

	bool Open( const std::vector<std::string>& vecstrGenes, const CDistanceMatrix& Dist ) {

		return CDat::Open( vecstrGenes, Dist ); }
};

class CPCLPair : public CPCLPairImpl {
public:
	bool Open( const char*, size_t );
	size_t Quantize( float, size_t ) const;
};

}

#endif // DATAPAIR_H
