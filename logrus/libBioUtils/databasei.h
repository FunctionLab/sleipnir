#ifndef DATABASEI_H
#define DATABASEI_H

#include <fstream>
#include <vector>

#include "typesi.h"

namespace libBioUtils {

class CDatabaselet {
public:
	CDatabaselet( ) : m_pfstm(NULL) { }

	~CDatabaselet( ) {

		if( m_pfstm )
			delete m_pfstm; }

	bool Open( const std::string&, const std::vector<std::string>&, uint32_t, uint32_t );
	bool Open( const std::string& );
	void Write( size_t, size_t, size_t, unsigned char );
	bool Get( size_t, size_t, std::vector<unsigned char>& ) const;

	size_t GetGene( const std::string& strGene ) const {
		size_t	i;

		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			if( strGene == m_vecstrGenes[ i ] )
				return i;

		return -1; }

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

private:
	size_t GetOffset( ) const {

		return ( m_iDatasets * m_iGenes ); }

	size_t GetOffset( size_t iGene ) const {

		return ( m_iHeader + ( GetOffset( ) * iGene ) ); }

	size_t GetOffset( size_t iOne, size_t iTwo ) const {

		return ( GetOffset( iOne ) + ( m_iDatasets * iTwo ) ); }

	size_t GetOffset( size_t iOne, size_t iTwo, size_t iDataset ) const {

		return ( GetOffset( iOne, iTwo ) + iDataset ); }

	uint32_t					m_iHeader;
	uint32_t					m_iGenes;
	uint32_t					m_iDatasets;
	std::vector<std::string>	m_vecstrGenes;
	std::fstream*				m_pfstm;
};

class CDatabaseImpl {
protected:
	static const char	c_acDAB[];
	static const char	c_acExtension[];

	std::vector<CDatabaselet>	m_vecDBs;
	size_t						m_iGenes;
};

}

#endif // DATABASEI_H
