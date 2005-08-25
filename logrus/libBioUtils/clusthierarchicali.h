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
	struct SHierarchy {
		CHierarchy*	m_pHier;
		bool		m_fUsed;
	};

	typedef std::vector<std::vector<float> >	TVecVecD;
	typedef std::pair<CHierarchy*,CHierarchy*>	TPrHier;

	static void FindHighestSimilarity( const CDistanceMatrix&, const vector<SHierarchy>&,
		const vector<SHierarchy>&, const TVecVecD&, const CDistanceMatrix&, size_t, float&,
		CHierarchy*&, CHierarchy*& );
	static void RefreshDistances( const CDistanceMatrix&, const vector<SHierarchy>&,
		const vector<SHierarchy>&, TVecVecD&, CDistanceMatrix&, size_t );
	static float CalculateDistance( const CHierarchy&, const CHierarchy&, const CDistanceMatrix&,
		const TVecVecD&, const CDistanceMatrix& );
	static void CalculateDistance( const CHierarchy&, const CHierarchy&, const CDistanceMatrix&,
		const TVecVecD&, const CDistanceMatrix&, size_t&, float& );
	static void CleanDistance( const CHierarchy&, vector<SHierarchy>&, vector<SHierarchy>&,
		TVecVecD& );
};

}

#endif // CLUSHIERARCHICALI_H
