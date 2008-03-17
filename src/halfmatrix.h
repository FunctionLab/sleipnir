#ifndef HALFMATRIX_H
#define HALFMATRIX_H

#include "halfmatrixi.h"

namespace Sleipnir {

/*!
 * \brief
 * A symmetric two-dimensional matrix.
 * 
 * \param tType
 * Type of element contained by the matrix.
 * 
 * \remarks
 * Nothing is bounds checked to increase efficiency.  Symmetric matrices consume half the space of a full
 * matrix and guarantee that element i,j is always identical to element j,i.  Everything is stored (and thus
 * accessed most efficiently) in row-major order.  Accessing elements on the diagonal generally doesn't work.
 * 
 * \see
 * CCompactMatrix | CFullMatrix
 */
template<class tType>
class CHalfMatrix : protected CHalfMatrixBase {
public:
	/*!
	 * \brief
	 * Return the number of entries in a symmetric matrix of the given size.
	 * 
	 * \param iSize
	 * Number of elements in a symmetric matrix for which the entries are calculated.
	 * 
	 * \returns
	 * Number of entries in a symmetric matrix of the given size.
	 * 
	 * \remarks
	 * A symmetric matrix over N elements always contains N(N-1)/2 entries.
	 */
	static size_t GetSpace( size_t iSize ) {

		return ( ( ( iSize * ( iSize - 1 ) ) / 2 ) * sizeof(tType) ); }

	CHalfMatrix( ) : m_aaData(NULL) { }

	virtual ~CHalfMatrix( ) {

		Reset( ); }

	/*!
	 * \brief
	 * Empties the matrix and deallocates all associated memory.
	 */
	void Reset( ) {
		size_t	i;

		if( m_aaData && m_fMemory ) {
			for( i = 0; ( i + 1 ) < m_iSize; ++i )
				delete[] m_aaData[ i ];
			delete[] m_aaData; }

		m_iSize = 0;
		m_aaData = NULL; }

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
	 * returned array will contain a number of elements equal to GetSize - iY - 1.
	 * 
	 * \see
	 * Set
	 */
	const tType* Get( size_t iY ) const {

		return m_aaData[ iY ]; }

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
	 * returned array will contain a number of elements equal to GetSize - iY - 1.
	 * 
	 * \see
	 * Set
	 */
	tType* Get( size_t iY ) {

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
		static tType	c_Zero	= 0;

		if( iX == iY )
			return c_Zero;

		HalfIndex( iY, iX );
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

		if( iX == iY )
			return;

		HalfIndex( iY, iX );
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
	 * and the given data array must contain at least GetSize - iY - 1 elements.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, const tType* aValues ) {

		memcpy( m_aaData[ iY ], aValues, sizeof(*m_aaData[ iY ]) * ( m_iSize - iY - 1 ) ); }

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
	void Initialize( const CHalfMatrix& Mat ) {

		Initialize( Mat.GetSize( ), Mat.m_aaData ); }

	/*!
	 * \brief
	 * Create a new matrix of the requested size and, optionally, referencing the given data.
	 * 
	 * \param iSize
	 * Matrix elements.
	 * 
	 * \param aaData
	 * If non-null, the memory that will back the newly created matrix.
	 * 
	 * \remarks
	 * A reference is held to the given memory; it is not copied, so it should not be destroyed before
	 * the newly initialized matrix.
	 */
	virtual void Initialize( size_t iSize, tType** aaData = NULL ) {
		size_t	i;

		Reset( );
		m_iSize = iSize;
		if( m_fMemory = !( m_aaData = aaData ) ) {
			m_aaData = iSize ? new tType*[ m_iSize - 1 ] : NULL;
			for( i = 0; ( i + 1 ) < m_iSize; ++i )
				m_aaData[ i ] = new tType[ m_iSize - i - 1 ]; } }

	/*!
	 * \brief
	 * Return the number of elements (row/columns) of the matrix.
	 * 
	 * \returns
	 * Number of elements of the matrix.
	 * 
	 * \see
	 * GetSpace
	 */
	size_t GetSize( ) const {

		return m_iSize; }

	/*!
	 * \brief
	 * Sets all entries of the matrix to 0 without changing its size.
	 */
	void Clear( ) {
		size_t	i;

		for( i = 0; ( i + 1 ) < m_iSize; ++i )
			memset( m_aaData[ i ], 0, ( m_iSize - i - 1 ) * sizeof(*m_aaData[ i ]) ); }

protected:
	/*!
	 * \brief
	 * True if the matrix is responsible for disposing of the underlying memory.
	 */
	bool	m_fMemory;
	/*!
	 * \brief
	 * Two-dimensional array backing the symmetric matrix.
	 */
	tType**	m_aaData;
};

/*!
 * \brief
 * A half matrix of four-byte floats, used for almost all Sleipnir continuously valued symmetric matrices.
 */
typedef CHalfMatrix<float>	CDistanceMatrix;

/*!
 * \brief
 * A special symmetric matrix in which each entry consumes exactly one bit.
 * 
 * \remarks
 * Sort of the logical extreme of CCompactMatrix, except symmetric and much more simple to implement.
 */
class CBinaryMatrix : public CHalfMatrix<unsigned char> {
public:
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
	bool Get( size_t iY, size_t iX ) const {

		if( iX == iY )
			return false;
		HalfIndex( iY, iX );
		return ( ( m_aaData[ iY ][ iX / 8 ] >> ( iX % 8 ) ) & 1 ); }

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
	 * \param fValue
	 * Value to store.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row and column must be smaller than the
	 * size.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, size_t iX, bool fValue ) {
		unsigned char	c;

		if( iX == iY )
			return;
		HalfIndex( iY, iX );
		c = 1 << (unsigned char)( iX % 8 );
		m_aaData[ iY ][ iX / 8 ] = ( m_aaData[ iY ][ iX / 8 ] & ~c ) | ( fValue ? c : 0 ); }

	/*!
	 * \brief
	 * Create a new matrix of the requested size and, optionally, referencing the given data.
	 * 
	 * \param iSize
	 * Matrix elements.
	 * 
	 * \param fClear
	 * If true, initialize the new matrix to contain only false values.
	 * 
	 * \param aabData
	 * If non-null, the memory that will back the newly created matrix.
	 * 
	 * \remarks
	 * A reference is held to the given memory; it is not copied, so it should not be destroyed before
	 * the newly initialized matrix.
	 */
	void Initialize( size_t iSize, bool fClear = false, unsigned char** aabData = NULL ) {
		size_t	i, iCount;

		Reset( );
		m_iSize = iSize;
		if( m_fMemory = !( m_aaData = aabData ) ) {
			m_aaData = new unsigned char*[ m_iSize - 1 ];
			for( i = 0; ( i + 1 ) < m_iSize; ++i ) {
				m_aaData[ i ] = new unsigned char[ iCount = ( ( ( m_iSize - i - 1 ) / 8 ) + 1 ) ];
				if( fClear )
					memset( m_aaData[ i ], 0, iCount * sizeof(*m_aaData[ i ]) ); } } }
};

}

#endif // HALFMATRIX_H
