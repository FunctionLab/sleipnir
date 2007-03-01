#ifndef HALFMATRIX_H
#define HALFMATRIX_H

#include "halfmatrixi.h"

namespace libBioUtils {

template<class tType>
class CHalfMatrix : protected CHalfMatrixBase {
public:
	static size_t GetSpace( size_t iSize ) {

		return ( ( ( iSize * ( iSize - 1 ) ) / 2 ) * sizeof(tType) ); }

	CHalfMatrix( ) : m_aaData(NULL) { }

	~CHalfMatrix( ) {

		Reset( ); }

	void Reset( ) {
		size_t	i;

		if( m_aaData && m_fMemory ) {
			for( i = 0; ( i + 1 ) < m_iSize; ++i )
				delete[] m_aaData[ i ];
			delete[] m_aaData; }

		m_iSize = 0;
		m_aaData = NULL; }

	const tType* Get( size_t iX ) const {

		return m_aaData[ iX ]; }

	tType& Get( size_t iX, size_t iY ) const {
		static tType	c_Zero	= 0;

		if( iX == iY )
			return c_Zero;

		HalfIndex( iX, iY );
		return m_aaData[ iX ][ iY ]; }

	tType& GetAbs( size_t iIndex ) const {
		size_t	iX, iY;

		UnHalfIndex( iIndex, iX, iY );
		return m_aaData[ iX ][ iY ]; }

	void Set( size_t iX, size_t iY, const tType& Value ) {

		if( iX == iY )
			return;

		HalfIndex( iX, iY );
		m_aaData[ iX ][ iY ] = Value; }

	void Set( size_t iX, const tType* adValues ) {

		memcpy( m_aaData[ iX ], adValues, sizeof(*m_aaData[ iX ]) * ( m_iSize - iX - 1 ) ); }

	void SetAbs( size_t iIndex, const tType& Value ) {
		size_t	iX, iY;

		UnHalfIndex( iIndex, iX, iY );
		m_aaData[ iX ][ iY ] = Value; }

	void Initialize( size_t iSize, const CHalfMatrix& Mat ) {

		Initialize( iSize, Mat.m_aaData ); }

	virtual void Initialize( size_t iSize, tType** aaData = NULL ) {
		size_t	i;

		Reset( );
		m_iSize = iSize;
		if( m_fMemory = !( m_aaData = aaData ) ) {
			m_aaData = iSize ? new tType*[ m_iSize - 1 ] : NULL;
			for( i = 0; ( i + 1 ) < m_iSize; ++i )
				m_aaData[ i ] = new tType[ m_iSize - i - 1 ]; } }

	size_t GetSize( ) const {

		return m_iSize; }

protected:
	bool	m_fMemory;
	tType**	m_aaData;
};

typedef CHalfMatrix<float>	CDistanceMatrix;

class CBinaryMatrix : public CHalfMatrix<unsigned char> {
public:
	bool Get( size_t iX, size_t iY ) const {

		if( iX == iY )
			return false;
		HalfIndex( iX, iY );
		return ( ( m_aaData[ iX ][ iY / 8 ] >> ( iY % 8 ) ) & 1 ); }

	void Set( size_t iX, size_t iY, bool fValue ) {
		unsigned char	c;

		if( iX == iY )
			return;
		HalfIndex( iX, iY );
		c = 1 << (unsigned char)( iY % 8 );
		m_aaData[ iX ][ iY / 8 ] = ( m_aaData[ iX ][ iY / 8 ] & ~c ) | ( fValue ? c : 0 ); }

	void Initialize( size_t, unsigned char** = NULL );
};

class CCompactMatrix : CCompactMatrixImpl {
public:
	bool Open( std::istream& );
	const unsigned char* Open( const unsigned char* );
	void Save( std::ostream& ) const;
	unsigned char Get( size_t, size_t ) const;
	void Set( size_t, size_t, unsigned char );
	void Initialize( size_t, unsigned char, bool = false );
	size_t GetSize( ) const;
};

}

#endif // HALFMATRIX_H
