#ifndef DAT_H
#define DAT_H

#include <iostream>
#include <string>
#include <vector>

#include "dati.h"

namespace libBioUtils {

class CGenes;

class CDat : protected CDatImpl {
public:
	enum EFilter {
		EFilterInclude	= 0,
		EFilterTerm		= EFilterInclude + 1,
		EFilterExclude	= EFilterTerm + 1
	};

	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetGene( const std::string& ) const;
	bool Open( const char* );
	bool Open( std::istream&, bool );
	bool Open( const CSlim& );
	bool Open( const CSlim&, const CSlim& );
	bool Open( const std::vector<std::string>& );
	bool Open( const std::vector<std::string>&, const CDistanceMatrix& );
	bool OpenGenes( std::istream&, bool );
	void Save( std::ostream&, bool ) const;
	void Normalize( bool = true );
	void Invert( );
	bool FilterGenes( const char*, EFilter );
	void FilterGenes( const CGenes&, EFilter );

	float Get( size_t iX, size_t iY ) const {

		return CDatImpl::Get( iX, iY ); }

	size_t GetGenes( ) const {

		return CDatImpl::GetGenes( ); }

	const CDistanceMatrix& Get( ) const {

		return m_Data; }

	bool Set( size_t iX, size_t iY, float dValue ) {

		return CDatImpl::Set( iX, iY, dValue ); }

	std::string GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	const CDistanceMatrix& GetData( ) const {

		return m_Data; }
};

}

#endif // DAT_H
