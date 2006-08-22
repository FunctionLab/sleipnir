#ifndef PCL_H
#define PCL_H

#include <iostream>
#include <string>

#include "pcli.h"

namespace libBioUtils {

class CPCL : CPCLImpl {
public:
	static size_t GetSkip( );

	void Open( const CPCL& );
	bool Open( std::istream& );
	bool Open( std::istream&, size_t );
	void Open( const std::vector<size_t>&, const std::vector<string>&,
		const std::vector<string>& );
	void Open( const std::vector<string>&, const std::vector<string>& );
	void Save( std::ostream&, const std::vector<size_t>* = NULL ) const;
	void SaveGene( std::ostream&, size_t, size_t = -1 ) const;
	void SaveHeader( std::ostream&, bool = false ) const;
	void SortGenes( const std::vector<size_t>& );
	void RankTransform( );
	bool AddGenes( const std::vector<std::string>& );
	void Normalize( );
	void Reset( );

	float Get( size_t iGene, size_t iExp ) const {

		return m_Data.Get( iGene, iExp ); }

	const float* Get( size_t iGene ) const {

		return m_Data.Get( iGene ); }

	const CDataMatrix& Get( ) const {

		return m_Data; }

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	const std::vector<std::string>& GetGeneNames( ) const {

		return m_vecstrGenes; }

	size_t GetExperiments( ) const {

		return m_vecstrExperiments.size( ); }

	const std::string& GetFeature( size_t iGene, size_t iFeature ) const {

		return m_vecvecstrFeatures[ iFeature ][ iGene ]; }

	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	const std::string& GetExperiment( size_t iExp ) const {

		return m_vecstrExperiments[ iExp ]; }

	void MaskGene( size_t iGene ) {

	//	m_setiGenes.push_back( iGene ); }
		m_setiGenes.insert( iGene ); }

	void Set( size_t iX, size_t iY, float dValue ) {

		m_Data.Set( iX, iY, dValue ); }

	size_t GetGene( const std::string& strGene ) const {
		size_t	i;

		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			if( m_vecstrGenes[ i ] == strGene )
				return i;

		return -1; }

	std::string GetFeature( size_t iGene, const char* szFeature ) const {
		size_t	 i;

		for( i = 0; i < m_vecstrFeatures.size( ); ++i )
			if( m_vecstrFeatures[ i ] == szFeature )
				return GetFeature( iGene, i );

		return ""; }
};

}

#endif // PCL_H
