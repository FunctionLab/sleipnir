#ifndef SPARSEMATRIXI_H
#define SPARSEMATRIXI_H

#include <map>

namespace Sleipnir {

template<class tType>
class CSparseMatrixImpl {
protected:
	CSparseMatrixImpl( const tType& Default ) : m_iR(0), m_Default(Default) { }

	size_t GetRows( ) const {

		return m_iR; }

	const tType& GetDefault( ) const {

		return m_Default; }

	size_t	m_iR;
	tType	m_Default;
};

}

#endif // SPARSEMATRIXI_H
