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

	CGenome* GetGenome( const CGene& Gene ) const {
		TMapGeneGenome::const_iterator	iterGenome;

		return ( ( ( iterGenome = m_mapGenes.find( (CGene*)&Gene ) ) == m_mapGenes.end( ) ) ? NULL :
			iterGenome->second ); }

	size_t GetGenomes( ) const {

		return m_vecpGenomes.size( ); }

	CGenome* GetGenome( size_t iOrganism ) const {

		return m_vecpGenomes[ iOrganism ]; }

	const std::string& GetOrganism( size_t iOrganism ) const {

		return m_vecstrOrganisms[ iOrganism ]; }
};

}

#endif // ORTHOLOGY_H
