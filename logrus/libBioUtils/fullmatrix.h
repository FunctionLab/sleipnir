#ifndef FULLMATRIX_H
#define FULLMATRIX_H

#include <sstream>
#include <vector>

#include "pstdint.h"

namespace libBioUtils {

template<class tType>
class CFullMatrix {
public:
	CFullMatrix( ) {

		m_iR = m_iC = 0;
		m_aaData = NULL; }

	~CFullMatrix( ) {

		Reset( ); }

	void Reset( ) {
		size_t	i;

		if( m_aaData && m_fMemory ) {
			for( i = 0; i < m_iR; ++i )
				delete[] m_aaData[ i ];
			delete[] m_aaData; }

		m_iR = m_iC = 0;
		m_aaData = NULL; }

	tType** const Get( ) const {

		return m_aaData; }

	tType* Get( size_t iY ) const {

		return m_aaData[ iY ]; }

	tType& Get( size_t iY, size_t iX ) const {

		return m_aaData[ iY ][ iX ]; }

	void Set( size_t iY, size_t iX, const tType& Value ) {

		m_aaData[ iY ][ iX ] = Value; }

	void Initialize( size_t iR, size_t iC, const CFullMatrix& Mat ) {

		Initialize( iR, iC, Mat.m_aaData ); }

	virtual void Initialize( size_t iR, size_t iC, tType** aaData = NULL ) {
		size_t	i;

		Reset( );
		m_iR = iR;
		m_iC = iC;
		if( m_fMemory = !( m_aaData = aaData ) ) {
			m_aaData = new tType*[ m_iR ];
			for( i = 0; i < m_iR; ++i )
				m_aaData[ i ] = new tType[ m_iC ]; } }

	size_t GetRows( ) const {

		return m_iR; }

	size_t GetColumns( ) const {

		return m_iC; }

	bool Open( std::istream& istm, bool fBinary ) {

		return ( fBinary ? OpenBinary( istm ) : OpenText( istm ) ); }

	bool Save( std::ostream& ostm, bool fBinary ) const {

		return ( fBinary ? SaveBinary( ostm ) : SaveText( ostm ) ); }

	void Clear( ) {
		size_t	i;

		if( !m_aaData )
			return;

		for( i = 0; i < m_iR; ++i )
			memset( m_aaData[ i ], 0, sizeof(*m_aaData[ i ]) * m_iC ); }

	bool AddRows( size_t iR ) {
		tType**	m_aaNew;
		size_t	i;

		if( !m_fMemory )
			return false;

		m_aaNew = new tType*[ m_iR + iR ];
		memcpy( m_aaNew, m_aaData, m_iR * sizeof(*m_aaData) );
		delete[] m_aaData;
		m_aaData = m_aaNew;
		for( i = 0; i < iR; ++i )
			m_aaNew[ m_iR + i ] = new tType[ m_iC ];
		m_iR += iR;

		return true; }

protected:
	static const size_t	c_iBuffer	= 8192;

	bool OpenBinary( std::istream& istm ) {
		size_t		i, j;
		uint32_t	iSize;

		Reset( );
		istm.read( (char*)&iSize, sizeof(iSize) );
		m_iR = iSize;
		istm.read( (char*)&iSize, sizeof(iSize) );
		m_iC = iSize;
		m_fMemory = true;
		m_aaData = new tType*[ m_iR ];
		for( i = 0; i < m_iR; ++i ) {
			m_aaData[ i ] = new tType[ m_iC ];
			for( j = 0; j < m_iC; ++j )
				istm.read( (char*)&m_aaData[ i ][ j ], sizeof(m_aaData[ i ][ j ]) ); }

		return true; }

	bool OpenText( std::istream& istm ) {
		char				szLine[ c_iBuffer ];
		std::vector<tType>	vecRow;
		std::vector<tType*>	vecpRows;
		std::stringstream	sstm;
		tType				Cur;
		tType*				aCur;
		size_t				i;

		Reset( );
		m_fMemory = true;
		while( istm.peek( ) != EOF ) {
			istm.getline( szLine, c_iBuffer - 1 );
			if( !szLine[ 0 ] )
				break;
			sstm.clear( );
			sstm.str( szLine );
			vecRow.clear( );
			while( sstm.peek( ) != EOF ) {
				sstm >> Cur;
				vecRow.push_back( Cur ); }
			if( !m_iC )
				m_iC = vecRow.size( );
			else if( m_iC != vecRow.size( ) )
				return false;

			aCur = new tType[ m_iC ];
			for( i = 0; i < m_iC; ++i )
				aCur[ i ] = vecRow[ i ];
			vecpRows.push_back( aCur ); }

		m_aaData = new tType*[ m_iR = vecpRows.size( ) ];
		for( i = 0; i < m_iR; ++i )
			m_aaData[ i ] = vecpRows[ i ];

		return true; }

	bool SaveBinary( std::ostream& ostm ) const {
		size_t		i, j;
		uint32_t	iSize;

		iSize = (uint32_t)m_iR;
		ostm.write( (const char*)&iSize, sizeof(iSize) );
		iSize = (uint32_t)m_iC;
		ostm.write( (const char*)&iSize, sizeof(iSize) );
		for( i = 0; i < m_iR; ++i )
			for( j = 0; j < m_iC; ++j )
				ostm.write( (const char*)&m_aaData[ i ][ j ], sizeof(m_aaData[ i ][ j ]) );

		return true; }

	bool SaveText( std::ostream& ostm ) const {
		size_t	i, j;

		if( !( m_iR && m_iC ) )
			return false;

		for( i = 0; i < m_iR; ++i ) {
			ostm << m_aaData[ i ][ 0 ];
			for( j = 1; j < m_iC; ++j )
				ostm << '\t' << m_aaData[ i ][ j ];
			ostm << std::endl; }

		return true; }

	size_t	m_iR;
	size_t	m_iC;
	bool	m_fMemory;
	tType**	m_aaData;
};

typedef CFullMatrix<float>	CDataMatrix;

}

#endif // FULLMATRIX_H
