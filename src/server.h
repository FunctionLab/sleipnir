#ifndef SERVER_H
#define SERVER_H

#include "serveri.h"

namespace Sleipnir {

class CServer : CServerImpl {
public:
	bool Initialize( const char*, IServerClient* );
	bool Initialize( size_t, size_t, IServerClient* );
	bool Start( );
	void Listen( );
	void Stop( );
};

}

#endif // SERVER_H
