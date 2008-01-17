#ifndef COMPACTMATRIX_H
#define COMPACTMATRIX_H

#include "compactmatrixi.h"

namespace libBioUtils {

class CCompactMatrix : CCompactMatrixImpl {
public:
	bool Open( std::istream& );
	const unsigned char* Open( const unsigned char* );
	void Save( std::ostream& ) const;
	void Randomize( );

	void Set( size_t iX, size_t iY, unsigned char cValue ) {

		HalfIndex( iX, iY );
		CCompactMatrixBase::Set( iX, iY, cValue ); }

	void Initialize( size_t iSize, unsigned char cValues, bool fClear ) {

		m_iSize = (uint32_t)iSize;
		return CCompactMatrixBase::Initialize( cValues, fClear ); }

	size_t GetSize( ) const {

		return m_iSize; }

	unsigned char Get( size_t iX, size_t iY ) const {

		if( iX == iY )
			return 0;

		HalfIndex( iX, iY );
		return CCompactMatrixBase::Get( iX, iY ); }
};

class CCompactFullMatrix : CCompactFullMatrixImpl {
public:
	void Initialize( size_t iRows, size_t iColumns, unsigned char cValues, bool fClear ) {

		m_iRows = (uint32_t)iRows;
		m_iColumns = (uint32_t)iColumns;
		CCompactMatrixBase::Initialize( cValues, fClear ); }

	size_t GetRows( ) const {

		return m_iRows; }

	size_t GetColumns( ) const {

		return m_iColumns; }

	unsigned char Get( size_t iY, size_t iX ) const {

		return CCompactMatrixBase::Get( iY, iX ); }

	void Set( size_t iX, size_t iY, unsigned char cValue ) {

		CCompactMatrixBase::Set( iX, iY, cValue ); }
};

}

#endif // COMPACTMATRIX_H
