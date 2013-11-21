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
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef COMPACTMATRIX_H
#define COMPACTMATRIX_H

#include "compactmatrixi.h"

namespace Sleipnir {

/*!
 * \brief
 * Store a discrete symmetric matrix using the fewest possible bytes.
 * 
 * Store a discrete symmetric matrix (also referred to as a half or triangular matrix) using the fewest
 * possible bytes.  This is done using sub-byte alignment, with each value occupying the smallest number of
 * bits possible to store the total number of possible different values.
 * 
 * \see
 * CCompactFullMatrix | CHalfMatrix
 */
class CCompactMatrix : CCompactMatrixImpl {
public:
	bool Open( std::istream& istm );
	const unsigned char* Open( const unsigned char* pbData );
	void Save( std::ostream& ostm ) const;

	/*!
	 * \brief
	 * Pseudorandomize the compact half matrix.
	 * 
	 * \remarks
	 * This doesn't really do a great job randomizing the data, since it ignores sub-byte alignment and
	 * just shuffles bytes, but it's good enough in a pinch.
	 */
	void Randomize( ) {

		if( m_aiData )
			std::random_shuffle( m_aiData, m_aiData + CountWords( ) ); }

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
	 * \param cValue
	 * Value to store.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given value must be smaller than the maximum
	 * value with which the matrix was initialized, and the row and column must be smaller than the size.
	 * As a symmetric matrix, the value at position XY will always equal the value at position YX.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, size_t iX, unsigned char cValue ) {

		if( iX == iY )
			return;

		HalfIndex( iY, iX );
		CCompactMatrixBase::Set( iY, iX, cValue ); }

	/*!
	 * \brief
	 * Initialize a new compact half matrix with the requested numbers of elements and discrete values.
	 * 
	 * \param iSize
	 * Number of elements (rows/columns) in the matrix.  As a symmetric matrix, the number of rows and
	 * columns is equal, so a single size is provided.
	 * 
	 * \param cValues
	 * Number of different values that each entry in the matrix can take.
	 * 
	 * \param fClear
	 * If true, set each matrix entry to zero after allocation.
	 * 
	 * Allocates and, optionally, clears enough memory to store a symmetric matrix with the requested
	 * number of elements taking the requested number of different values.  A symmetric matrix with iSize
	 * elements will have iSize * (iSize - 1) / 2 entries, each of which takes ceil(log2(cValues)) bits.
	 * This means that roughly ceil(4 * log2(cValues) * iSize * (iSize - 1)) bytes will be allocated for
	 * the entire half matrix.
	 * 
	 * \remarks
	 * cValues should be greater than one.
	 */
	void Initialize( size_t iSize, unsigned char cValues, bool fClear ) {

		m_iSize = (uint32_t)iSize;
		CCompactMatrixBase::Initialize( cValues, fClear ); }

	/*!
	 * \brief
	 * Returns the number of elements in the matrix.
	 * 
	 * \returns
	 * Number of elements in the matrix.
	 * 
	 * \remarks
	 * Since a symmetric matrix must be square, the number of rows equals the number of columns and is thus
	 * referred to as the number of elements.
	 */
	size_t GetSize( ) const {

		return m_iSize; }

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
	 * size.  As a symmetric matrix, the value at position XY will always equal the value at position YX.
	 * 
	 * \see
	 * Set
	 */
	unsigned char Get( size_t iY, size_t iX ) const {

		if( iX == iY )
			return 0;

		HalfIndex( iY, iX );
		return CCompactMatrixBase::Get( iY, iX ); }
};

/*!
 * \brief
 * Store a discrete matrix using the fewest possible bytes.
 * 
 * Store a discrete matrix using the fewest possible bytes.  This is done using sub-byte alignment, with
 * each value occupying the smallest number of bits possible to store the total number of possible different values.
 * 
 * \see
 * CCompactMatrix | CFullMatrix
 */
class CCompactFullMatrix : CCompactFullMatrixImpl {
public:
	/*!
	 * \brief
	 * Initialize a new compact matrix with the requested numbers of rows, columns, and discrete values.
	 * 
	 * \param iRows
	 * Number of rows in the matrix.
	 * 
	 * \param iColumns
	 * Number of columns in the matrix.
	 * 
	 * \param cValues
	 * Number of different values that each entry in the matrix can take.
	 * 
	 * \param fClear
	 * If true, set each matrix entry to zero after allocation.
	 * 
	 * Allocates and, optionally, clears enough memory to store a matrix with the requested number of rows
	 * and columns, with each entry taking the requested number of different values.  A matrix with iRows
	 * rows and iColumns columns will have iRows * iColumns entries, each of which takes ceil(log2(cValues))
	 * bits.  This means that roughly ceil(8 * log2(cValues) * iRows * iColumns) bytes will be allocated for
	 * the entire matrix.
	 * 
	 * \remarks
	 * cValues should be greater than one.
	 */
	void Initialize( size_t iRows, size_t iColumns, unsigned char cValues, bool fClear ) {

		m_iRows = (uint32_t)iRows;
		m_iColumns = (uint32_t)iColumns;
		CCompactMatrixBase::Initialize( cValues, fClear ); }

	/*!
	 * \brief
	 * Returns the number of rows in the matrix.
	 * 
	 * \returns
	 * Number of rows in the matrix.
	 * 
	 * \see
	 * GetColumns
	 */
	size_t GetRows( ) const {

		return m_iRows; }

	/*!
	 * \brief
	 * Returns the number of columns in the matrix.
	 * 
	 * \returns
	 * Number of columns in the matrix.
	 * 
	 * \see
	 * GetRows
	 */
	size_t GetColumns( ) const {

		return m_iColumns; }

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
	unsigned char Get( size_t iY, size_t iX ) const {

		return CCompactMatrixBase::Get( iY, iX ); }

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
	 * \param cValue
	 * Value to store.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given value must be smaller than the maximum
	 * value with which the matrix was initialized, and the row and column must be smaller than the size.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, size_t iX, unsigned char cValue ) {

		CCompactMatrixBase::Set( iY, iX, cValue ); }
};

}

#endif // COMPACTMATRIX_H
