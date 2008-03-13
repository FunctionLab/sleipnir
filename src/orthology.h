#ifndef ORTHOLOGY_H
#define ORTHOLOGY_H

#include "orthologyi.h"

namespace Sleipnir {

class COrthology : COrthologyImpl {
public:
	bool Open( std::istream& );
	void Save( std::ostream& ) const;

	size_t GetClusters( ) const {

		return m_vecvecpGenes.size( ); }

	size_t GetGenes( size_t iCluster ) const {

		return m_vecvecpGenes[ iCluster ].size( ); }

	CGene& GetGene( size_t iCluster, size_t iGene ) const {

		return *m_vecvecpGenes[ iCluster ][ iGene ]; }

	size_t GetGenomes( ) const {

		return m_vecpGenomes.size( ); }

	CGenome* GetGenome( size_t iOrganism ) const {

		return m_vecpGenomes[ iOrganism ]; }

	const std::string& GetOrganism( size_t iOrganism ) const {

		return m_vecstrOrganisms[ iOrganism ]; }

	size_t GetOrganism( const CGene& Gene ) const {
		TMapGeneI::const_iterator	iterGenome;

		return ( ( ( iterGenome = m_mapGenes.find( (CGene*)&Gene ) ) == m_mapGenes.end( ) ) ? NULL :
			iterGenome->second ); }
};

}

#endif // ORTHOLOGY_H
