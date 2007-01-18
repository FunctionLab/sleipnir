#ifndef ORTHOLOGY_H
#define ORTHOLOGY_H

#include "orthologyi.h"

namespace libBioUtils {

class COrthology : COrthologyImpl {
public:
	bool Open( std::istream& );

	size_t GetClusters( ) const {

		return m_vecvecpGenes.size( ); }

	size_t GetGenes( size_t iCluster ) const {

		return m_vecvecpGenes[ iCluster ].size( ); }

	CGene& GetGene( size_t iCluster, size_t iGene ) const {

		return *m_vecvecpGenes[ iCluster ][ iGene ]; }
};

}

#endif // ORTHOLOGY_H
