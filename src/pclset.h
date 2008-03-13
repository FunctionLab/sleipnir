#ifndef PCLSET_H
#define PCLSET_H

#include "pclseti.h"
#include "meta.h"
#include "pcl.h"

namespace Sleipnir {

class CPCLSet : CPCLSetImpl {
public:
	bool Open( const std::vector<std::string>&, size_t = 2, bool = false );

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	size_t GetPCLs( ) const {

		return m_iPCLs; }

	float Get( size_t iPCL, size_t iGene, size_t iExp ) const {
		size_t	iMap;

		if( ( iMap = m_Genes.Get( iPCL, iGene ) ) == -1 )
			return CMeta::GetNaN( );

		return m_aPCLs[ iPCL ].Get( iMap, iExp ); }

	const float* Get( size_t iPCL, size_t iGene ) const {
		size_t	iMap;

		if( ( iMap = m_Genes.Get( iPCL, iGene ) ) == -1 )
			return NULL;

		return m_aPCLs[ iPCL ].Get( iMap ); }

	size_t GetGene( const std::string& strGene ) const {
		TMapStrI::const_iterator	iterGene;

		return ( ( ( iterGene = m_mapGenes.find( strGene ) ) == m_mapGenes.end( ) ) ? -1 :
			iterGene->second ); }

	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	const std::vector<std::string>& GetGeneNames( ) const {

		return m_vecstrGenes; }

	const CPCL& Get( size_t iPCL ) const {

		return m_aPCLs[ iPCL ]; }
};

}

#endif // PCLSET_H
