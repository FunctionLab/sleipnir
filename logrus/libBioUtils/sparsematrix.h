#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <map>

namespace libBioUtils {

template<class tType>
class CSparseMatrix {
public:
	CSparseMatrix( tType Default ) : m_iR(0), m_amapData(NULL), m_Default(Default) { }

	~CSparseMatrix( ) {

		Reset( ); }

	void Reset( ) {

		if( m_amapData )
			delete[] m_amapData;
		m_iR = 0; }

	void Initialize( size_t iR ) {

		Reset( );
		m_amapData = new std::map<size_t,tType>[ m_iR = iR ]; }

	tType Get( size_t iY, size_t iX ) const {
		std::map<size_t,tType>::const_iterator	iter;

		return ( ( ( iter = m_amapData[ iY ].find( iX ) ) == m_amapData[ iY ].end( ) ) ? m_Default : iter->second ); }

	void Set( size_t iY, size_t iX, const tType& Value ) {

		m_amapData[ iY ][ iX ] = Value; }

protected:
	size_t					m_iR;
	std::map<size_t,tType>*	m_amapData;
	tType					m_Default;
};

}

#endif // SPARSEMATRIX_H
