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
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <map>

#include "sparsematrixi.h"

namespace Sleipnir {

/*!
 * \brief
 * An asymmetric two-dimensional sparse matrix using maps for each row.
 * 
 * \param tType
 * Type of element contained by the matrix.
 * 
 * Implements a two-dimensional matrix assuming few non-default entries in each row.  These entries are
 * maintained by a map for each row, allowing fairly rapid lookup at the expense of moderate memory usage.
 * 
 * \see
 * CSparseListMatrix | CFullMatrix
 */
template<class tType>
class CSparseMapMatrix : protected CSparseMatrixImpl<tType> {
public:
	/*!
	 * \brief
	 * Create a new sparse matrix with the given default value.
	 * 
	 * \param Default
	 * Default value provided for entries not in the matrix.
	 */
	CSparseMapMatrix( const tType& Default ) : CSparseMatrixImpl( Default ), m_amapData(NULL) { }

	~CSparseMapMatrix( ) {

		Reset( ); }

	/*!
	 * \brief
	 * Empties the matrix and deallocates all associated memory.
	 */
	void Reset( ) {

		if( m_amapData )
			delete[] m_amapData;
		m_iR = 0; }

	/*!
	 * \brief
	 * Create a new matrix of the requested size.
	 * 
	 * \param iR
	 * Matrix rows.
	 * 
	 * \remarks
	 * Columns are not specified since they are allocated dynamically by each row's map.
	 */
	void Initialize( size_t iR ) {

		Reset( );
		m_amapData = new std::map<size_t,tType>[ m_iR = iR ]; }

	/*!
	 * \brief
	 * Return the number of rows in the matrix.
	 * 
	 * \returns
	 * Number of rows in the matrix.
	 */
	size_t GetRows( ) const {

		return CSparseMatrixImpl::GetRows( ); }

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
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than GetRows.  If
	 * no entry is found for the given position, the default value is returned.
	 * 
	 * \see
	 * Set
	 */
	tType Get( size_t iY, size_t iX ) const {
		typename std::map<size_t,tType>::const_iterator	iter;

		return ( ( ( iter = m_amapData[ iY ].find( iX ) ) == m_amapData[ iY ].end( ) ) ? m_Default :
			iter->second ); }

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
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than GetRows.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, size_t iX, const tType& Value ) {

		m_amapData[ iY ][ iX ] = Value; }

	/*!
	 * \brief
	 * Return the default value for entries not in the matrix.
	 * 
	 * \returns
	 * Default value for entries not in the matrix.
	 */
	const tType& GetDefault( ) const {

		return CSparseMatrixImpl::GetDefault( ); }

private:
	std::map<size_t,tType>*	m_amapData;
};

/*!
 * \brief
 * An asymmetric two-dimensional sparse matrix using linked lists for each row.
 * 
 * \param tType
 * Type of element contained by the matrix.
 * 
 * Implements a two-dimensional matrix assuming few non-default entries in each row.  These entries are
 * maintained by a linked list for each row, allowing very low memory usage (for sufficiently sparse
 * matrices) at the expense of potentially slow lookup.
 * 
 * \see
 * CSparseMapMatrix | CFullMatrix
 */
template<class tType>
class CSparseListMatrix : protected CSparseMatrixImpl<tType> {
private:
	struct SNode {
		size_t	m_iX;
		tType	m_Value;
		SNode*	m_pNext;

		SNode( size_t iX, tType Value, SNode* pNext ) : m_iX(iX), m_Value(Value), m_pNext(pNext) { }
	};

public:
	/*!
	 * \brief
	 * Create a new sparse matrix with the given default value.
	 * 
	 * \param Default
	 * Default value provided for entries not in the matrix.
	 */
	CSparseListMatrix( const tType& Default ) : CSparseMatrixImpl( Default ), m_apsData(NULL) { }

	~CSparseListMatrix( ) {

		Reset( ); }

	/*!
	 * \brief
	 * Empties the matrix and deallocates all associated memory.
	 */
	void Reset( ) {
		size_t	i;
		SNode*	pNode;
		SNode*	pNext;

		if( m_apsData ) {
			for( i = 0; i < m_iR; ++i )
				for( pNode = m_apsData[ i ]; pNode; pNode = pNext ) {
					pNext = pNode->m_pNext;
					delete pNode; }
			delete[] m_apsData; }
		m_iR = 0; }

	/*!
	 * \brief
	 * Create a new matrix of the requested size.
	 * 
	 * \param iR
	 * Matrix rows.
	 * 
	 * \remarks
	 * Columns are not specified since they are allocated dynamically by each row's map.
	 */
	void Initialize( size_t iR ) {

		Reset( );
		m_apsData = new SNode*[ m_iR = iR ];
		memset( m_apsData, 0, m_iR * sizeof(*m_apsData) ); }

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
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than GetRows.  If
	 * no entry is found for the given position, the default value is returned.
	 * 
	 * \see
	 * Set
	 */
	tType Get( size_t iY, size_t iX ) const {
		SNode*	pNode;

		for( pNode = m_apsData[ iY ]; pNode && ( pNode->m_iX >= iX ); pNode = pNode->m_pNext )
			if( pNode->m_iX == iX )
				return pNode->m_Value;

		return m_Default; }

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
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than GetRows.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, size_t iX, const tType& Value ) {
		SNode*	pNode;

		if( !( pNode = m_apsData[ iY ] ) || ( iX > pNode->m_iX ) ) {
			m_apsData[ iY ] = new SNode( iX, Value, pNode );
			return; }
		for( ; pNode->m_pNext; pNode = pNode->m_pNext ) {
			if( pNode->m_iX == iX ) {
				pNode->m_Value = Value;
				return; }
			if( iX > pNode->m_pNext->m_iX ) {
				pNode->m_pNext = new SNode( iX, Value, pNode->m_pNext );
				return; } }
		if( pNode->m_iX == iX )
			pNode->m_Value = Value;
		else
			pNode->m_pNext = new SNode( iX, Value, NULL ); }

	/*!
	 * \brief
	 * Return the default value for entries not in the matrix.
	 * 
	 * \returns
	 * Default value for entries not in the matrix.
	 */
	const tType& GetDefault( ) const {

		return CSparseMatrixImpl::GetDefault( ); }

	/*!
	 * \brief
	 * Return the number of rows in the matrix.
	 * 
	 * \returns
	 * Number of rows in the matrix.
	 */
	size_t GetRows( ) const {

		return CSparseMatrixImpl::GetRows( ); }

private:
	SNode**	m_apsData;
};

}

#endif // SPARSEMATRIX_H
