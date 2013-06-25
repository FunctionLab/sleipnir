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
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <map>

#include "sparsematrixi.h"

namespace Sleipnir {

template <typename tType>
struct CPair {
	utype i;
	tType v;
};

template <typename tType>
struct CAscendingKey{
	bool operator()(const CPair<tType>& lx, const CPair<tType>& rx) const{
		if(lx.i < rx.i) return true;
		if(lx.i > rx.i) return false;
		if(lx.v < rx.v) return true;
		return false;
	}
};

template <typename tType>
struct CAscendingValue{
	bool operator()(const CPair<tType> &lx, const CPair<tType> &rx) const{
		if(lx.v < rx.v) return true;
		if(lx.v > rx.v) return false;
		if(lx.i < rx.i) return true;
		return false;
	}
};

template <typename tType>
struct CDescendingKey{
	bool operator()(const CPair<tType>& lx, const CPair<tType>& rx) const{
		if(lx.i < rx.i) return false;
		if(lx.i > rx.i) return true;
		if(lx.v < rx.v) return false;
		return true;
	}
};

template <typename tType>
struct CDescendingValue{
	bool operator()(const CPair<tType> &lx, const CPair<tType> &rx) const{
		if(lx.v < rx.v) return false;
		if(lx.v > rx.v) return true;
		if(lx.i < rx.i) return false;
		return true;
	}
};


//A full matrix
template<class tType>
class CSparseFlatMatrix : protected CSparseMatrixImpl<tType> {
public:
	CSparseFlatMatrix(const tType& Default): CSparseMatrixImpl<tType>(Default){}
	~CSparseFlatMatrix(){
		Reset();
	}
	void Initialize(size_t iR){
		Reset();
		CSparseMatrixImpl<tType>::m_iR = iR;
		m_vecData.resize(iR);
		m_currentIndex.resize(iR);
		m_bSorted.resize(iR);
	}
	//num is initial capacity
	void InitializeRow(size_t rowID, size_t num){
		m_vecData[rowID] = vector<CPair<tType> >();
		m_vecData[rowID].reserve(num);
		m_currentIndex[rowID] = 0;
		m_bSorted[rowID] = false;
	}
	void Reset() {
		size_t i;
		for(i=0; i<CSparseMatrixImpl<tType>::m_iR; i++)
			m_vecData[i].clear();
		m_vecData.clear();
		CSparseMatrixImpl<tType>::m_iR = 0; 
	}
	const tType& GetDefault() const {
		return CSparseMatrixImpl<tType>::GetDefault();
	}
	size_t GetRows() const {
		return CSparseMatrixImpl<tType>::GetRows(); 
	}
	//does not check if [iY][iX] already exists
	//nor if the row m_vecData[iY] is already full
	void Add(size_t iY, size_t iX, tType v){
		CPair<tType> cp;
		cp.i = (utype) iX;
		cp.v = v;
		m_vecData[iY].push_back(cp);
		m_bSorted[iY] = false;
		m_currentIndex[iY]++;
	}
	const vector<CPair<tType> >& GetRow(size_t iY) const{
		return m_vecData[iY];
	}
	typename vector<CPair<tType> >::iterator RowBegin(size_t iY){
		return m_vecData[iY].begin();
	}
	typename vector<CPair<tType> >::iterator RowEnd(size_t iY){
		return m_vecData[iY].end();
	}
	void Shrink(){
		size_t i;
		for(i=0; i<m_vecData.size(); i++)
			m_vecData[i].resize(m_currentIndex[i]);
	}
	void Organize(){
		size_t i;
		for(i=0; i<m_vecData.size(); i++){
			sort(m_vecData[i].begin(), m_vecData[i].end(), CAscendingKey<tType>());	
			m_bSorted[i] = true;
		}
	}
	bool Check(size_t iY, size_t iX){
		if(!m_bSorted[iY])
			SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		if(ind!=(size_t) -1) 
			return true;
		return false;
	}
	//assume element at this coordinate exists
	tType Get(size_t iY, size_t iX) const {
		if(!m_bSorted[iY])
			SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		return m_vecData[iY][ind].v;
	}
	//assume element at this coordinate exists
	void Set(size_t iY, size_t iX, tType v) {
		if(!m_bSorted[iY])
			SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		m_vecData[iY][ind].v = v;
	}
	void Increment(size_t iY, size_t iX, tType v){
		if(!m_bSorted[iY])
			SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		m_vecData[iY][ind].v += v;
	}

private:
	vector<vector<CPair<tType> > > m_vecData;
	//parent class contains m_iR and m_Default
	vector<bool> m_bSorted;
	vector<size_t> m_currentIndex;

	void SortRow(size_t iY){
		sort(m_vecData[iY].begin(), m_vecData[iY].end(), CAscendingKey<tType>());
		m_bSorted[iY] = true;
	}
	size_t GetIndex(size_t iY, size_t iX) {
		//suppose m_bSorted[iY] is true
		return BinarySearch(m_vecData[iY], iX);
	}
	size_t BinarySearch(vector<CPair<tType> > &A, size_t iX){
		int iMax = A.size()-1;
		int iMin = 0;
		while(iMax>=iMin){
			int iMid = (iMin + iMax) / 2;
			if(A[iMid].i < iX)
				iMin = iMid + 1;
			else if(A[iMid].i > iX)
				iMax = iMid - 1;
			else
				return (size_t) iMid;
		}
		//fprintf(stderr, "Element %d is not found!\n", iX);
		return (size_t) -1;
	}
};

//A half-matrix
template<class tType>
class CSparseFlatHalfMatrix : protected CSparseMatrixImpl<tType> {
public:
	CSparseFlatHalfMatrix(const tType& Default): CSparseMatrixImpl<tType>(Default){}
	~CSparseFlatHalfMatrix(){
		Reset();
	}
	void Copy(const CSparseFlatMatrix<tType> &cf){ //a full matrix (symmetric)
		CSparseMatrixImpl<tType>::m_Default = cf.GetDefault();
		Initialize(cf.GetRows());
		size_t i,j;
		for(i=0; i<CSparseMatrixImpl<tType>::m_iR; i++){
			const vector<CPair<tType> > &allR = cf.GetRow(i);
			InitializeRow(i, allR.size());
			for(j=0; j<allR.size(); j++)
				if(allR[j].i > i)
					Add(i, allR[j].i, allR[j].v);
		}
		Organize();
	}
	void Initialize(size_t iR){
		Reset();
		CSparseMatrixImpl<tType>::m_iR = iR;
		m_vecData.resize(iR);
		m_currentIndex.resize(iR);
		m_bSorted.resize(iR);
	}
	//num is initial capacity
	void InitializeRow(size_t rowID, size_t num){
		m_vecData[rowID] = vector<CPair<tType> >();
		m_vecData[rowID].reserve(num);
		m_currentIndex[rowID] = 0;
		m_bSorted[rowID] = false;
	}
	void Reset() {
		size_t i;
		for(i=0; i<CSparseMatrixImpl<tType>::m_iR; i++)
			m_vecData[i].clear();
		m_vecData.clear();
		CSparseMatrixImpl<tType>::m_iR = 0; 
	}
	const tType& GetDefault() const {
		return CSparseMatrixImpl<tType>::GetDefault();
	}
	size_t GetRows() const {
		return CSparseMatrixImpl<tType>::GetRows(); 
	}
	//does not check if [iY][iX] already exists
	//nor if the row m_vecData[iY] is already full
	void Add(size_t iY, size_t iX, tType v){
		AdjustCoord(iY, iX);
		CPair<tType> cp;
		cp.i = (utype) iX;
		cp.v = v;
		m_vecData[iY].push_back(cp);
		m_bSorted[iY] = false;
		m_currentIndex[iY]++;
	}
	const vector<CPair<tType> >& GetRow(size_t iY) const{
		return m_vecData[iY];
	}
	typename vector<CPair<tType> >::iterator RowBegin(size_t iY){
		return m_vecData[iY].begin();
	}
	typename vector<CPair<tType> >::iterator RowEnd(size_t iY){
		return m_vecData[iY].end();
	}
	void Shrink(){
		size_t i;
		for(i=0; i<m_vecData.size(); i++)
			m_vecData[i].resize(m_currentIndex[i]);
	}
	void SortRow(size_t iY){
		sort(m_vecData[iY].begin(), m_vecData[iY].end(), CAscendingKey<tType>());
		m_bSorted[iY] = true;
	}
	void Organize(){
		size_t i;
		for(i=0; i<m_vecData.size(); i++)
			SortRow(i);
	}
	void AdjustCoord(size_t &iY, size_t &iX){
		if(iY>=iX){ //second must be the greater
			size_t tmp = iY;
			iY = iX;
			iX = tmp;
		}
	}
	bool Check(size_t iY, size_t iX){
		AdjustCoord(iY, iX);
		if(!m_bSorted[iY]) SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		if(ind!=(size_t) -1) return true;
		return false;
	}
	CPair<tType>* GetElement(size_t iY, size_t iX){
		AdjustCoord(iY, iX);
		if(!m_bSorted[iY]) SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		if(ind==(size_t)-1) return NULL;
		return &m_vecData[iY][ind];
	}
	//assume element at this coordinate exists
	tType Get(size_t iY, size_t iX) const {
		AdjustCoord(iY, iX);
		if(!m_bSorted[iY]) SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		return m_vecData[iY][ind].v;
	}
	//assume element at this coordinate exists
	void Set(size_t iY, size_t iX, tType v) {
		AdjustCoord(iY, iX);
		if(!m_bSorted[iY]) SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		m_vecData[iY][ind].v = v;
	}
	void Increment(size_t iY, size_t iX, tType v){
		AdjustCoord(iY, iX);
		if(!m_bSorted[iY]) SortRow(iY);
		size_t ind = GetIndex(iY, iX);
		m_vecData[iY][ind].v += v;
	}

private:
	//parent class contains m_iR and m_Default
	vector<vector<CPair<tType> > > m_vecData;
	vector<bool> m_bSorted;
	vector<size_t> m_currentIndex;

	size_t GetIndex(size_t iY, size_t iX) {
		return BinarySearch(m_vecData[iY], iX); //suppose m_bSorted[iY] = true
	}
	size_t BinarySearch(vector<CPair<tType> > &A, size_t iX){
		int iMax = A.size()-1;
		int iMin = 0;
		while(iMax>=iMin){
			int iMid = (iMin + iMax) / 2;
			if(A[iMid].i < iX)
				iMin = iMid + 1;
			else if(A[iMid].i > iX)
				iMax = iMid - 1;
			else
				return (size_t) iMid;
		}
		return (size_t) -1;
	}
};
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
	CSparseMapMatrix( const tType& Default ) : CSparseMatrixImpl<tType>( Default ), m_amapData(NULL) { }

	~CSparseMapMatrix( ) {

		Reset( ); }

	/*!
	 * \brief
	 * Empties the matrix and deallocates all associated memory.
	 */
	void Reset( ) {

		if( m_amapData )
			delete[] m_amapData;
		CSparseMatrixImpl<tType>::m_iR = 0; }

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
		m_amapData = new std::map<size_t,tType>[ CSparseMatrixImpl<tType>::m_iR = iR ]; }

	/*!
	 * \brief
	 * Return the number of rows in the matrix.
	 * 
	 * \returns
	 * Number of rows in the matrix.
	 */
	size_t GetRows( ) const {

		return CSparseMatrixImpl<tType>::GetRows( ); }

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

		return ( ( ( iter = m_amapData[ iY ].find( iX ) ) == m_amapData[ iY ].end( ) ) ?
			CSparseMatrixImpl<tType>::m_Default : iter->second ); }

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

		return CSparseMatrixImpl<tType>::GetDefault( ); }

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
	CSparseListMatrix( const tType& Default ) : CSparseMatrixImpl<tType>( Default ), m_apsData(NULL) { }

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
			for( i = 0; i < CSparseMatrixImpl<tType>::m_iR; ++i )
				for( pNode = m_apsData[ i ]; pNode; pNode = pNext ) {
					pNext = pNode->m_pNext;
					delete pNode; }
			delete[] m_apsData; }
		CSparseMatrixImpl<tType>::m_iR = 0; }

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
		m_apsData = new SNode*[ CSparseMatrixImpl<tType>::m_iR = iR ];
		memset( m_apsData, 0, CSparseMatrixImpl<tType>::m_iR * sizeof(*m_apsData) ); }

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

		return CSparseMatrixImpl<tType>::m_Default; }

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

		return CSparseMatrixImpl<tType>::GetDefault( ); }

	/*!
	 * \brief
	 * Return the number of rows in the matrix.
	 * 
	 * \returns
	 * Number of rows in the matrix.
	 */
	size_t GetRows( ) const {

		return CSparseMatrixImpl<tType>::GetRows( ); }

private:
	SNode**	m_apsData;
};

}

#endif // SPARSEMATRIX_H
