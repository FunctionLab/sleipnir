#ifndef CLUSHIERARCHICAL_H
#define CLUSHIERARCHICAL_H

#include "clusthierarchicali.h"

namespace Sleipnir {

class CPCL;
class IMeasure;

class CHierarchy : public CHierarchyImpl {
public:
	CHierarchy( size_t, float, const CHierarchy*, const CHierarchy* );

	const CHierarchy& Get( bool ) const;
	size_t GetID( ) const;
	float GetSimilarity( ) const;
	bool IsGene( ) const;
	void Save( ostream&, size_t ) const;
	void GetGenes( std::vector<size_t>& ) const;
	float SortChildren( const std::vector<float>& );
	void Destroy( );
};

class CClustHierarchical : CClustHierarchicalImpl {
public:
	static CHierarchy* Cluster( const CDistanceMatrix& );
	static CHierarchy* Cluster( const CDistanceMatrix&, const std::vector<bool>& );
};

}

#endif // CLUSHIERARCHICAL_H
