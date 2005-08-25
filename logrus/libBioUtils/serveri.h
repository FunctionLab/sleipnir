#ifndef SERVERI_H
#define SERVERI_H

#include "propertyfile.h"
#include "thread.h"

namespace libBioUtils {

class CServer;
class IServerClient;

class CServerImpl {
protected:
	typedef std::map<IServerClient*, CThread>	TMapThreads;
	typedef std::vector<IServerClient*>			TVecClients;

	static const char*	c_szPort;
	static const char*	c_szTimeout;
	static CServer*		s_pServer;

#ifndef WIN32
	static void Alarm( int );
#endif // WIN32

	CPropertyFile	m_Config;
	int				m_iPort;
	int				m_iTimeout;
	SOCKET			m_iSocket;
	bool			m_fStop;
	IServerClient*	m_pClient;
};

}

#endif // SERVERI_H
