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
		EFilterDiscard	= EFilterInclude + 1,
		EFilterNegate	= EFilterDiscard + 1,
		EFilterExclude	= EFilterNegate + 1
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
	void Normalize( );
	void FilterGenes( const CGenes&, EFilter );

	float Get( size_t iX, size_t iY ) const {

		return CDatImpl::Get( iX, iY ); }

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	const CDistanceMatrix& Get( ) const {

		return m_Data; }

	bool Set( size_t iX, size_t iY, float dValue ) {

		return CDatImpl::Set( iX, iY, dValue ); }

	string GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	const CDistanceMatrix& GetData( ) const {

		return m_Data; }
};

}

#endif // DAT_H
