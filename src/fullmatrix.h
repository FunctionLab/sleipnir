/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef FULLMATRIX_H
#define FULLMATRIX_H

#include <sstream>
#include <vector>

#include "typesi.h"

namespace Sleipnir {

/*!
 * \brief
 * An asymmetric two-dimensional matrix.
 * 
 * \param tType
 * Type of element contained by the matrix.
 * 
 * \remarks
 * Essentially nothing except a lack of bounds checking for efficiency's sake separates this from every
 * other 2D matrix implementation known to man.
 * 
 * \see
 * CCompactMatrix | CHalfMatrix
 */
template<class tType>
class CFullMatrix {
public:
	CFullMatrix( ) {

		m_iR = m_iC = 0;
		m_aaData = NULL; }

	virtual ~CFullMatrix( ) {

		Reset( ); }

	/*!
	 * \brief
	 * Empties the matrix and deallocates all associated memory.
	 */
	void Reset( ) {
		size_t	i;

		if( m_aaData && m_fMemory ) {
			for( i = 0; i < GetRows( ); ++i )
				delete[] m_aaData[ i ];
			delete[] m_aaData; }

		m_iR = m_iC = 0;
		m_aaData = NULL; }

	/*!
	 * \brief
	 * Return the two-dimensional array backing the matrix.
	 * 
	 * \returns
	 * Two-dimensional array backing the matrix.
	 */
	tType** const Get( ) const {

		return m_aaData; }

	/*!
	 * \brief
	 * Return a single row of the matrix.
	 * 
	 * \param iY
	 * Matrix row.
	 * 
	 * \returns
	 * Requested matrix row.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than the size.  The
	 * returned array will contain a number of elements equal to the number of columns.
	 * 
	 * \see
	 * Set
	 */
	tType* Get( size_t iY ) const {

		return m_aaData[ iY ]; }

	/*!
	 * \brief
	 * Returns the value at the requested matrix position.
	 * 
	 * \param iY
	 * Matrix row.
	 * 
	 * \param iX
	 * Matrix column.
	 * 
	 * \returns
	 * Value at the requested matrix position.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row and column must be smaller than the
	 * size.
	 * 
	 * \see
	 * Set
	 */
	tType& Get( size_t iY, size_t iX ) const {

		return m_aaData[ iY ][ iX ]; }

	/*!
	 * \brief
	 * Set the value at the requested matrix position.
	 * 
	 * \param iY
	 * Matrix row.
	 * 
	 * \param iX
	 * Matrix column.
	 * 
	 * \param Value
	 * Value to store.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row and column must be smaller than the
	 * size.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, size_t iX, const tType& Value ) {

		m_aaData[ iY ][ iX ] = Value; }

	/*!
	 * \brief
	 * Set a single row of the matrix.
	 * 
	 * \param iY
	 * Matrix row.
	 * 
	 * \param aValues
	 * Data to be copied into the requested row.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than the size,
	 * and the given data array must be at least as long as the number of columns.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, const tType* aValues ) {

		memcpy( m_aaData[ iY ], aValues, GetColumns( ) ); }

	/*!
	 * \brief
	 * Create a new matrix using a reference to the given matrix.
	 * 
	 * \param Mat
	 * Matrix to duplicate in the created matrix.
	 * 
	 * \remarks
	 * A reference is held to the given matrix; it is not copied, so it should not be destroyed before
	 * the newly initialized matrix.
	 */
	void Initialize( const CFullMatrix& Mat ) {

		Initialize( Mat.GetRows( ), Mat.GetColumns( ), Mat.m_aaData ); }

	/*!
	 * \brief
	 * Create a new matrix of the requested size and, optionally, referencing the given data.
	 * 
	 * \param iR
	 * Matrix rows.
	 * 
	 * \param iC
	 * Matrix columns.
	 * 
	 * \param aaData
	 * If non-null, the memory that will back the newly created matrix.
	 * 
	 * \remarks
	 * A reference is held to the given memory; it is not copied, so it should not be destroyed before
	 * the newly initialized matrix.
	 */
	virtual void Initialize( size_t iR, size_t iC, tType** aaData = NULL ) {
		size_t	i;

		Reset( );
		m_iR = iR;
		m_iC = iC;
		if( m_fMemory = !( m_aaData = aaData ) ) {
			m_aaData = new tType*[ GetRows( ) ];
			for( i = 0; i < GetRows( ); ++i )
				m_aaData[ i ] = new tType[ GetColumns( ) ]; } }

	/*!
	 * \brief
	 * Return the number of rows in the matrix.
	 * 
	 * \returns
	 * Number of rows in the matrix.
	 * 
	 * \see
	 * GetColumns
	 */
	size_t GetRows( ) const {

		return m_iR; }

	/*!
	 * \brief
	 * Return the number of columns in the matrix.
	 * 
	 * \returns
	 * Number of columns in the matrix.
	 * 
	 * \see
	 * GetRows
	 */
	size_t GetColumns( ) const {

		return m_iC; }

	/*!
	 * \brief
	 * Loads a matrix from the given stream.
	 * 
	 * \param istm
	 * Stream from which matrix is loaded.
	 * 
	 * \param fBinary
	 * If true, matrix and stream are in binary format.
	 * 
	 * \returns
	 * True if matrix was loaded successfully.
	 * 
	 * \remarks
	 * A binary matrix is saved in row-major order and should have been produced using Save.  A text matrix
	 * is a simple tab-delimited file from which matrix elements are read using >>.
	 */
	bool Open( std::istream& istm, bool fBinary ) {

		return ( fBinary ? OpenBinary( istm ) : OpenText( istm ) ); }

	/*!
	 * \brief
	 * Saves a matrix to the given stream in either binary or tab-delimited text format.
	 * 
	 * \param ostm
	 * Stream to which matrix is saved.
	 * 
	 * \param fBinary
	 * If true, matrix is saved in binary format; otherwise, matrix is saved as tab-delimited text.
	 * 
	 * \returns
	 * True if matrix was saved successfully.
	 * 
	 * \remarks
	 * A binary matrix is saved in row-major order; a text matrix is a simple tab-delimited file from which
	 * matrix elements are read using >>.
	 * 
	 * \see
	 * Open
	 */
	bool Save( std::ostream& ostm, bool fBinary ) const {

		return ( fBinary ? SaveBinary( ostm ) : SaveText( ostm ) ); }

	/*!
	 * \brief
	 * Saves a matrix to the given stream in either binary or tab-delimited text format.
	 * 
	 * \param szFile
	 * File to which matrix is saved.
	 * 
	 * \param fBinary
	 * If true, matrix is saved in binary format; otherwise, matrix is saved as tab-delimited text.
	 * 
	 * \returns
	 * True if matrix was saved successfully.
	 * 
	 * \remarks
	 * A binary matrix is saved in row-major order; a text matrix is a simple tab-delimited file from which
	 * matrix elements are read using >>.
	 * 
	 * \see
	 * Open
	 */
	bool Save( const char* szFile, bool fBinary ) const {
		std::ofstream	ofsm;

		ofsm.open( szFile, ( fBinary ? ios_base::binary : ios_base::out ) | ios_base::out );
		return Save( ofsm, fBinary ); }

	/*!
	 * \brief
	 * Sets all entries of the matrix to 0 without changing its size.
	 */
	void Clear( ) {
		size_t	i;

		if( !m_aaData )
			return;

		for( i = 0; i < GetRows( ); ++i )
			memset( m_aaData[ i ], 0, sizeof(*m_aaData[ i ]) * GetColumns( ) ); }

	/*!
	 * \brief
	 * Increases the number of rows in the matrix by the requested amount.
	 * 
	 * \param iR
	 * Number of rows to append to the matrix.
	 * 
	 * \returns
	 * True if the rows were appended successfully.
	 * 
	 * \remarks
	 * On success, GetRows will be increased by the requested amount.  The news rows are not initialized
	 * to any particular value.
	 */
	bool AddRows( size_t iR ) {
		tType**	m_aaNew;
		size_t	i;

		if( !m_fMemory )
			return false;

		m_aaNew = new tType*[ GetRows( ) + iR ];
		memcpy( m_aaNew, m_aaData, GetRows( ) * sizeof(*m_aaData) );
		delete[] m_aaData;
		m_aaData = m_aaNew;
		for( i = 0; i < iR; ++i )
			m_aaNew[ GetRows( ) + i ] = new tType[ GetColumns( ) ];
		m_iR += iR;

		return true; }

private:
	static const size_t	c_iBuffer	= 131072;

	bool OpenBinary( std::istream& istm ) {
		size_t		i;
		uint32_t	iSize;

		Reset( );
		if( !istm.good( ) || istm.eof( ) )
			return false;
		istm.read( (char*)&iSize, sizeof(iSize) );
		m_iR = iSize;
		istm.read( (char*)&iSize, sizeof(iSize) );
		m_iC = iSize;
		m_fMemory = true;
		m_aaData = new tType*[ GetRows( ) ];
		for( i = 0; i < GetRows( ); ++i ) {
			m_aaData[ i ] = new tType[ GetColumns( ) ];
			istm.read( (char*)m_aaData[ i ], (std::streamsize)( sizeof(*m_aaData[ i ]) * GetColumns( ) ) ); }

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
			if( !GetColumns( ) )
				m_iC = vecRow.size( );
			else if( GetColumns( ) != vecRow.size( ) )
				return false;

			aCur = new tType[ GetColumns( ) ];
			for( i = 0; i < GetColumns( ); ++i )
				aCur[ i ] = vecRow[ i ];
			vecpRows.push_back( aCur ); }

		m_aaData = new tType*[ m_iR = vecpRows.size( ) ];
		for( i = 0; i < GetRows( ); ++i )
			m_aaData[ i ] = vecpRows[ i ];

		return true; }

	bool SaveBinary( std::ostream& ostm ) const {
		size_t		i;
		uint32_t	iSize;

		iSize = (uint32_t)GetRows( );
		ostm.write( (const char*)&iSize, sizeof(iSize) );
		iSize = (uint32_t)GetColumns( );
		ostm.write( (const char*)&iSize, sizeof(iSize) );
		for( i = 0; i < GetRows( ); ++i )
			ostm.write( (const char*)m_aaData[ i ], (std::streamsize)( sizeof(*m_aaData[ i ]) *
				GetColumns( ) ) );

		return true; }

	bool SaveText( std::ostream& ostm ) const {
		size_t	i, j;

		if( !( GetRows( ) && GetColumns( ) ) )
			return false;

		for( i = 0; i < GetRows( ); ++i ) {
			ostm << m_aaData[ i ][ 0 ];
			for( j = 1; j < GetColumns( ); ++j )
				ostm << '\t' << m_aaData[ i ][ j ];
			ostm << std::endl; }

		return true; }

	size_t	m_iR;
	size_t	m_iC;
	bool	m_fMemory;
	tType**	m_aaData;
};

/*!
 * \brief
 * A full matrix of four-byte floats, used for almost all Sleipnir continuously valued asymmetric matrices.
 */
typedef CFullMatrix<float>	CDataMatrix;

}

#endif // FULLMATRIX_H
