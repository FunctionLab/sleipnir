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

	enum EFormat {
		EFormatBinary	= 0,
		EFormatText		= EFormatBinary + 1,
		EFormatPCL		= EFormatText + 1,
		EFormatSparse	= EFormatPCL + 1
	};

	bool Open( const char*, bool = false );
	bool Open( std::istream&, EFormat = EFormatBinary, float = HUGE_VAL );
	bool Open( const CSlim& );
	bool Open( const CSlim&, const CSlim& );
	bool Open( const std::vector<std::string>&, bool = true, const char* = NULL );
	bool Open( const std::vector<std::string>&, const CDistanceMatrix& );
	bool Open( const std::vector<CGenes*>&, const std::vector<CGenes*>&, float, const CGenome& );
	bool Open( const CDat&, const std::vector<CGenes*>&, const CGenome&, bool );
	bool Open( const CPCL&, const IMeasure*, bool );

	bool OpenGenes( std::istream&, bool, bool = false );
	void Save( std::ostream&, EFormat = EFormatBinary ) const;
	void Save( const char* ) const;
	void SaveDOT( std::ostream&, float = HUGE_VAL, const CGenome* = NULL, bool = false ) const;
	void SaveGDF( std::ostream&, float = HUGE_VAL ) const;
	void SaveNET( std::ostream&, float = HUGE_VAL ) const;
	void Normalize( bool = true );
	void Invert( );
	void Rank( );
	bool FilterGenes( const char*, EFilter );
	void FilterGenes( const CGenes&, EFilter );

	size_t GetGene( const std::string& strGene ) const {

		return CDatImpl::GetGene( strGene ); }

	float& Get( size_t iX, size_t iY ) const {

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

	void Set( size_t iX, const float* adValues ) {

		m_Data.Set( iX, adValues ); }

	const float* Get( size_t iX ) const {

		return m_Data.Get( iX ); }
};

}

#endif // DAT_H
