#ifndef PCL_H
#define PCL_H

#include <iostream>
#include <string>

#include "pcli.h"

namespace libBioUtils {

class CPCL : CPCLImpl {
public:
	static size_t GetSkip( );

	float Get( size_t, size_t ) const;
	const float* Get( size_t ) const;
	const CDataMatrix& Get( ) const;
	const std::string& GetExperiment( size_t ) const;
	size_t GetExperiments( ) const;
	std::string GetFeature( size_t, const char* ) const;
	const std::string& GetFeature( size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGene( const std::string& ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetGenes( ) const;
	void MaskGene( size_t );
	void Open( const CPCL& );
	bool Open( std::istream&, size_t );
	void Open( const std::vector<size_t>&, const std::vector<string>&,
		const std::vector<string>& );
	void Open( const std::vector<string>&, const std::vector<string>& );
	void Save( std::ostream&, const std::vector<size_t>* = NULL ) const;
	void SaveGene( std::ostream&, size_t, size_t = -1 ) const;
	void SaveHeader( std::ostream&, bool = false ) const;
	void Set( size_t, size_t, float );
	void SortGenes( const std::vector<size_t>& );
	void RankTransform( );
};

}

#endif // PCL_H
