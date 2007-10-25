#ifndef DATABASEI_H
#define DATABASEI_H

#define DATABASE_NIBBLES

#include <fstream>
#include <vector>

#include "typesi.h"

namespace libBioUtils {

class CDatabaselet {
public:
	CDatabaselet( );
	~CDatabaselet( );

	bool Open( const std::string&, const std::vector<std::string>&, uint32_t, uint32_t );
	bool Open( const std::string& );
	bool Open( const std::vector<std::string>&, const std::vector<std::string>&, bool );
	void Write( size_t, size_t, size_t, unsigned char );
	bool Get( size_t, size_t, std::vector<unsigned char>& ) const;
	bool Get( size_t, std::vector<unsigned char>& ) const;

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

private:
	size_t GetDatasets( ) const {

		return ( ( m_iDatasets
#ifdef DATABASE_NIBBLES
			+ 1 ) / 2
#else // DATABASE_NIBBLES
			)
#endif // DATABASE_NIBBLES
			); }

	size_t GetDataset( size_t iDataset ) const {

		return ( iDataset
#ifdef DATABASE_NIBBLES
			/ 2
#endif // DATABASE_NIBBLES
			); }

	size_t GetOffset( ) const {

		return ( GetDatasets( ) * m_iGenes ); }

	size_t GetOffset( size_t iGene ) const {

		return ( m_iHeader + ( GetOffset( ) * iGene ) ); }

	size_t GetOffset( size_t iOne, size_t iTwo ) const {

		return ( GetOffset( iOne ) + ( GetDatasets( ) * iTwo ) ); }

	size_t GetOffset( size_t iOne, size_t iTwo, size_t iDataset ) const {

		return ( GetOffset( iOne, iTwo ) + GetDataset( iDataset ) ); }

	uint32_t					m_iHeader;
	uint32_t					m_iGenes;
	uint32_t					m_iDatasets;
	std::vector<std::string>	m_vecstrGenes;
	mutable std::fstream		m_fstm;
	mutable pthread_mutex_t*	m_pmutx;
};

class CDatabaseImpl {
protected:
	static const char	c_acDAB[];
	static const char	c_acExtension[];

	~CDatabaseImpl( ) {

		Clear( ); }

	bool Open( const std::vector<std::string>&, const std::string&, bool, const std::vector<std::string>& );

	void Clear( ) {
		size_t	i;

		m_mapstriGenes.clear( );
		for( i = 0; i < m_vecpDBs.size( ); ++i )
			delete m_vecpDBs[ i ];
		m_vecpDBs.clear( ); }

	size_t GetGene( const std::string& strGene ) const {
		std::map<std::string, size_t>::const_iterator	iterGene;

		return ( ( ( iterGene = m_mapstriGenes.find( strGene ) ) == m_mapstriGenes.end( ) ) ? -1 :
			iterGene->second ); }

	std::vector<CDatabaselet*>		m_vecpDBs;
	std::map<std::string, size_t>	m_mapstriGenes;
};

}

#endif // DATABASEI_H
