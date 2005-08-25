#ifndef DAT_H
#define DAT_H

#include <iostream>
#include <string>
#include <vector>

#include "dati.h"

namespace libBioUtils {

class CGenes;
class CSlim;

class CDat : protected CDatImpl {
public:
	enum EFilter {
		EFilterInclude	= 0,
		EFilterDiscard	= EFilterInclude + 1,
		EFilterNegate	= EFilterDiscard + 1,
		EFilterExclude	= EFilterNegate + 1
	};

	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetGenes( ) const;
	size_t GetGene( const std::string& ) const;
	std::string GetGene( size_t ) const;
	float Get( size_t, size_t ) const;
	const CDistanceMatrix& Get( ) const;
	bool Set( size_t, size_t, float );
	bool Open( const char* );
	bool Open( std::istream&, bool );
	bool Open( const CSlim& );
	bool Open( const std::vector<std::string>& );
	bool Open( const std::vector<std::string>&, const CDistanceMatrix& );
	bool OpenGenes( std::istream&, bool );
	void Save( std::ostream&, bool ) const;
	const CDistanceMatrix& GetData( ) const;
	void Normalize( );
	void FilterGenes( const CGenes&, EFilter );
};

}

#endif // DAT_H
