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

		return m_vecDBs[ iOne % m_vecDBs.size( ) ].Get( iOne / m_vecDBs.size( ), iTwo, vecbData ); }

	bool Get( size_t iGene, std::vector<unsigned char>& vecbData ) const {

		return m_vecDBs[ iGene % m_vecDBs.size( ) ].Get( iGene / m_vecDBs.size( ), vecbData ); }

	size_t GetGenes( ) const {

		return m_iGenes; }
};

}

#endif // DATABASE_H
