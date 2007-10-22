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

		return m_vecDBs[ iOne % m_vecDBs.size( ) ].Get( m_veciGenes[ iOne ], iTwo, vecbData ); }
};

}

#endif // DATABASE_H
