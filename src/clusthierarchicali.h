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
#ifndef CLUSHIERARCHICALI_H
#define CLUSHIERARCHICALI_H

#include <string>
#include <vector>

#include "fullmatrix.h"
#include "halfmatrix.h"

namespace Sleipnir {

class CHierarchy;

class CHierarchyImpl {
protected:
	~CHierarchyImpl( );

	std::string GetSave( ) const;
	bool Save( std::ostream&, size_t ) const;

	bool IsGene( ) const {

		return !( m_pLeft && m_pRight ); }

	size_t				m_iID;
	size_t				m_iWeight;
	float				m_dScore;
	const CHierarchy*	m_pLeft;
	const CHierarchy*	m_pRight;
};

class CClustHierarchicalImpl {
protected:

	static CHierarchy* Cluster( const CDistanceMatrix&, const std::vector<bool>* = NULL );
	static void AssertParentage( std::vector<size_t>&, std::vector<size_t>&, std::vector<size_t>&, size_t, size_t );
	static CHierarchy* ConstructHierarchy( const std::vector<size_t>&, const std::vector<size_t>&,
		const std::vector<float>&, size_t );
	static void UpdateDistances( size_t, size_t, CDistanceMatrix&, size_t, size_t, std::vector<float>&,
		std::vector<size_t>& );
};

}

#endif // CLUSHIERARCHICALI_H
