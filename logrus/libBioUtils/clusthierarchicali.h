#ifndef CLUSHIERARCHICALI_H
#define CLUSHIERARCHICALI_H

#include <string>
#include <vector>

#include "fullmatrix.h"
#include "halfmatrix.h"

namespace libBioUtils {

class CHierarchy;

class CHierarchyImpl {
protected:
	~CHierarchyImpl( );

	bool IsGene( ) const;
	std::string GetSave( ) const;
	bool Save( std::ostream&, size_t ) const;

	size_t				m_iID;
	size_t				m_iWeight;
	float				m_dScore;
	const CHierarchy*	m_pLeft;
	const CHierarchy*	m_pRight;
};

class CClustHierarchicalImpl {
protected:

	static void AssertParentage( std::vector<size_t>&, std::vector<size_t>&, std::vector<size_t>&, size_t, size_t );
	static CHierarchy* ConstructHierarchy( const std::vector<size_t>&, const std::vector<size_t>&,
		const std::vector<float>&, size_t );
};

}

#endif // CLUSHIERARCHICALI_H
