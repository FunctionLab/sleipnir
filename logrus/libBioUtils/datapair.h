#ifndef DATAPAIR_H
#define DATAPAIR_H

#include "datapairi.h"

namespace libBioUtils {

class CSlim;

class CDataPair : public CDataPairImpl {
public:
	bool Open( const char*, bool );
	bool Open( const CSlim& );
	void SetQuants( const float*, size_t );
	void SetQuants( const std::vector<float>& );
	size_t Quantify( float ) const;
	bool IsContinuous( ) const;
	unsigned char GetValues( ) const;

	bool Open( const std::vector<CGenes*>& vecpPositives, const CDat& DatNegatives,
		const CGenome& Genome ) {

		return CDat::Open( vecpPositives, DatNegatives, Genome ); }

	bool Open( const std::vector<std::string>& vecstrGenes, const CDistanceMatrix& Dist ) {

		return CDat::Open( vecstrGenes, Dist ); }
};

class CPCLPair : public CPCLPairImpl {
public:
	bool Open( const char*, size_t );
	size_t Quantify( float, size_t ) const;
};

}

#endif // DATAPAIR_H
