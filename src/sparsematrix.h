#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <map>

namespace Sleipnir {

template<class tType>
class CSparseMapMatrix {
public:
	CSparseMapMatrix( tType Default ) : m_iR(0), m_amapData(NULL), m_Default(Default) { }

	~CSparseMapMatrix( ) {

		Reset( ); }

	void Reset( ) {

		if( m_amapData )
			delete[] m_amapData;
		m_iR = 0; }

	void Initialize( size_t iR ) {

		Reset( );
		m_amapData = new std::map<size_t,tType>[ m_iR = iR ]; }

	tType Get( size_t iY, size_t iX ) const {
		typename std::map<size_t,tType>::const_iterator	iter;

		return ( ( ( iter = m_amapData[ iY ].find( iX ) ) == m_amapData[ iY ].end( ) ) ? m_Default : iter->second ); }

	void Set( size_t iY, size_t iX, const tType& Value ) {

		m_amapData[ iY ][ iX ] = Value; }

protected:
	size_t					m_iR;
	std::map<size_t,tType>*	m_amapData;
	tType					m_Default;
};

template<class tType>
class CSparseListMatrix {
private:
	struct SNode {
		size_t	m_iX;
		tType	m_Value;
		SNode*	m_pNext;

		SNode( size_t iX, tType Value, SNode* pNext ) : m_iX(iX), m_Value(Value), m_pNext(pNext) { }
	};

public:
	CSparseListMatrix( tType Default ) : m_iR(0), m_apsData(NULL), m_Default(Default) { }

	~CSparseListMatrix( ) {

		Reset( ); }

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

	void Initialize( size_t iR ) {

		Reset( );
		m_apsData = new SNode*[ m_iR = iR ];
		memset( m_apsData, 0, m_iR * sizeof(*m_apsData) ); }

	tType Get( size_t iY, size_t iX ) const {
		SNode*	pNode;

		for( pNode = m_apsData[ iY ]; pNode && ( pNode->m_iX >= iX ); pNode = pNode->m_pNext )
			if( pNode->m_iX == iX )
				return pNode->m_Value;

		return m_Default; }

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

protected:
	size_t	m_iR;
	SNode**	m_apsData;
	tType	m_Default;
};

}

#endif // SPARSEMATRIX_H
