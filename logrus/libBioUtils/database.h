#ifndef DATABASE_H
#define DATABASE_H

#include "databasei.h"

namespace libBioUtils {

class IBayesNet;

class CDatabase : CDatabaseImpl {
public:
	bool Open( const std::vector<std::string>&, const std::string&, const IBayesNet*, const std::string&,
		bool );
	bool Open( const std::vector<std::string>&, size_t, const std::string& );
	bool Get( size_t, size_t, std::vector<unsigned char>& ) const;
};

}

#endif // DATABASE_H
