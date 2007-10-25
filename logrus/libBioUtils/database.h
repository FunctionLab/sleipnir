#ifndef DATABASE_H
#define DATABASE_H

#include "databasei.h"

namespace libBioUtils {

class IBayesNet;

class CDatabase : CDatabaseImpl {
public:
	bool Open( const std::vector<std::string>&, const std::string&, const IBayesNet*, const std::string&,
		size_t, bool );
	bool Open( const std::string& );

	bool Get( size_t iOne, size_t iTwo, std::vector<unsigned char>& vecbData ) const {

		return m_vecpDBs[ iOne % m_vecpDBs.size( ) ]->Get( iOne / m_vecpDBs.size( ), iTwo, vecbData ); }

	bool Get( size_t iGene, std::vector<unsigned char>& vecbData ) const {

		return m_vecpDBs[ iGene % m_vecpDBs.size( ) ]->Get( iGene / m_vecpDBs.size( ), vecbData ); }

	size_t GetGenes( ) const {

		return m_mapstriGenes.size( ); }

	size_t GetGene( const std::string& strGene ) const {

		return CDatabaseImpl::GetGene( strGene ); }
};

}

#endif // DATABASE_H
