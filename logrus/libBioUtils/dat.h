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
		EFilterExclude	= EFilterTerm + 1,
		EFilterPixie	= EFilterExclude + 1,
		EFilterEdge		= EFilterPixie + 1
	};

	enum EFormat {
		EFormatBinary	= 0,
		EFormatText		= EFormatBinary + 1,
		EFormatPCL		= EFormatText + 1,
		EFormatSparse	= EFormatPCL + 1
	};

	static std::string GetColor( float dValue ) {
		float	adColor[ 3 ];
		char	acColor[ 16 ];
		size_t	i;

		for( i = 0; i < ARRAYSIZE(adColor); ++i )
			adColor[ i ] = ( dValue < 0.5 ) ? ( c_adColorMin[ i ] * ( 1 - ( 2 * dValue ) ) ) :
				( c_adColorMax[ i ] * 2 * ( dValue - 0.5f ) );
		sprintf_s( acColor, "%02X%02X%02X", (size_t)( adColor[ 0 ] * 255 ),
			(size_t)( adColor[ 1 ] * 255 ), (size_t)( adColor[ 2 ] * 255 ) );

		return acColor; }

	bool Open( const char*, bool = false, size_t = 2, bool = false );
	bool Open( std::istream&, EFormat = EFormatBinary, float = HUGE_VAL, bool = false, size_t = 2,
		bool = false );
// MEFIT OFF
	bool Open( const CSlim& );
	bool Open( const CSlim&, const CSlim& );
// MEFIT ON
	bool Open( const std::vector<std::string>&, bool = true, const char* = NULL );
	bool Open( const std::vector<std::string>&, const CDistanceMatrix& );
	bool Open( const std::vector<CGenes*>&, const std::vector<CGenes*>&, float, const CGenome& );
	bool Open( const CDat&, const std::vector<CGenes*>&, const CGenome&, bool );
	bool Open( const CPCL&, const IMeasure*, bool );

	bool OpenGenes( std::istream&, bool, bool = false );
	bool OpenGenes( const char*, size_t = 2 );
	void Save( std::ostream&, EFormat = EFormatBinary ) const;
	void Save( const char* ) const;
	void SaveDOT( std::ostream&, float = HUGE_VAL, const CGenome* = NULL, bool = false ) const;
	void SaveGDF( std::ostream&, float = HUGE_VAL ) const;
	void SaveNET( std::ostream&, float = HUGE_VAL ) const;
	void SaveMATISSE( std::ostream&, float = HUGE_VAL, const CGenome* = NULL ) const;
	void Normalize( bool = true );
	void Invert( );
	void Rank( );
	bool FilterGenes( const char*, EFilter, size_t = -1 );
	void FilterGenes( const CGenes&, EFilter, size_t = -1, bool = true );

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

		return CDatImpl::GetGeneNames( ); }

	void Set( size_t iX, const float* adValues ) {

		m_Data.Set( iX, adValues ); }

	const float* Get( size_t iX ) const {

		return m_Data.Get( iX ); }

	void SetGene( size_t iGene, const std::string& strGene ) {

		if( m_pPCL )
			m_pPCL->SetGene( iGene, strGene );
		else
			m_vecstrGenes[ iGene ] = strGene; }
};

}

#endif // DAT_H
