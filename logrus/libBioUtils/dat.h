#ifndef DAT_H
#define DAT_H

#include <iostream>
#include <string>
#include <vector>

#include "dati.h"

namespace libBioUtils {

class CGenes;
class CGenome;

class CDat : protected CDatImpl {
public:
	enum EFilter {
		EFilterInclude	= 0,
		EFilterTerm		= EFilterInclude + 1,
		EFilterExclude	= EFilterTerm + 1
	};

	size_t GetGene( const std::string& ) const;
	bool Open( const char* );
	bool Open( std::istream&, bool, bool = false, float = HUGE_VAL );
	bool Open( const CSlim& );
	bool Open( const CSlim&, const CSlim& );
	bool Open( const std::vector<std::string>& );
	bool Open( const std::vector<std::string>&, const CDistanceMatrix& );
	bool Open( const std::vector<CGenes*>&, const std::vector<CGenes*>&, float, const CGenome& );
	bool Open( const std::vector<CGenes*>&, const CDat&, const CGenome& );
	bool Open( const CPCL&, const IMeasure*, bool );
	bool OpenGenes( std::istream&, bool, bool = false );
	void Save( std::ostream&, bool ) const;
	void SaveDOT( std::ostream&, float = HUGE_VAL, const CGenome* = NULL, bool = false ) const;
	void SaveGDF( std::ostream&, float = HUGE_VAL ) const;
	void SaveNET( std::ostream&, float = HUGE_VAL ) const;
	void Normalize( bool = true );
	void Invert( );
	void Rank( );
	bool FilterGenes( const char*, EFilter );
	void FilterGenes( const CGenes&, EFilter );

	float Get( size_t iX, size_t iY ) const {

		return CDatImpl::Get( iX, iY ); }

	size_t GetGenes( ) const {

		return CDatImpl::GetGenes( ); }

	const CDistanceMatrix& Get( ) const {

		return m_Data; }

	CDistanceMatrix& Get( ) {

		return m_Data; }

	bool Set( size_t iX, size_t iY, float dValue ) {

		return CDatImpl::Set( iX, iY, dValue ); }

	std::string GetGene( size_t iGene ) const {

		return CDatImpl::GetGene( iGene ); }

	const CDistanceMatrix& GetData( ) const {

		return m_Data; }

	const std::vector<std::string>& GetGeneNames( ) const {

		return ( m_pMeasure ? m_pPCL->GetGeneNames( ) : m_vecstrGenes ); }
};

}

#endif // DAT_H
