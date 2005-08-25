#ifndef PCLSET_H
#define PCLSET_H

#include "pclseti.h"

namespace libBioUtils {

class CPCLSet : CPCLSetImpl {
public:
	bool Open( const std::vector<std::string>&, size_t = 2 );
	size_t GetGene( const std::string& ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	size_t GetPCLs( ) const;
	size_t GetExperiments( size_t ) const;
	float Get( size_t, size_t, size_t ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
};

}

#endif // PCLSET_H
