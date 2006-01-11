#ifndef HALFMATRIX_H
#define HALFMATRIX_H

#include "halfmatrixi.h"

namespace libBioUtils {

template<class tType>
class CHalfMatrix : protected CHalfMatrixBase {
public:
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

		HalfIndex( iX, iY );
		return m_aaData[ iX ][ iY ]; }

	void Set( size_t iX, size_t iY, const tType& Value ) {

		HalfIndex( iX, iY );
		m_aaData[ iX ][ iY ] = Value; }

	void Set( size_t iX, const tType* adValues ) {

		memcpy( m_aaData[ iX ], adValues, sizeof(*m_aaData[ iX ]) * ( m_iSize - iX - 1 ) ); }

	void Initialize( size_t iSize, const CHalfMatrix& Mat ) {

		Initialize( iSize, Mat.m_aaData ); }

	virtual void Initialize( size_t iSize, tType** aaData = NULL ) {
		size_t	i;

		Reset( );
		m_iSize = iSize;
		if( m_fMemory = !( m_aaData = aaData ) ) {
			m_aaData = new tType*[ m_iSize - 1 ];
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
	bool Get( size_t, size_t ) const;
	void Set( size_t, size_t, bool );
	void Initialize( size_t, unsigned char** = NULL );
};

class CCompactMatrix : CCompactMatrixImpl {
public:
	unsigned char Get( size_t, size_t ) const;
	void Set( size_t, size_t, unsigned char );
	void Initialize( size_t, unsigned char, bool = false );
	size_t GetSize( ) const;
};

}

#endif // HALFMATRIX_H
