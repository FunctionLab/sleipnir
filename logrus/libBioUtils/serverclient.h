#ifndef SERVERCLIENT_H
#define SERVERCLIENT_H

#include "typesi.h"

namespace libBioUtils {

class CPropertyFile;
class CServer;

class IServerClient {
public:
	virtual IServerClient* NewInstance( SOCKET, uint32_t, uint16_t, const CPropertyFile* ) = 0;
	virtual bool ProcessMessage(
#ifdef XML_ENABLED
		xercesc::DOMDocument*
#else // XML_ENABLED
		const std::vector<unsigned char>&
#endif // XML_ENABLED
		) = 0;
	virtual void Destroy( ) = 0;
};

}

#endif // SERVERCLIENT_H
