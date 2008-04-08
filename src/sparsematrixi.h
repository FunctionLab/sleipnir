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
