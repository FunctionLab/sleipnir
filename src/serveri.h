#ifndef SERVERI_H
#define SERVERI_H

namespace Sleipnir {

class CServer;
class IServerClient;

class CServerImpl {
protected:
	typedef std::map<IServerClient*, pthread_t>	TMapThreads;
	typedef std::vector<IServerClient*>			TVecClients;

	static const char*	c_szPort;
	static const char*	c_szTimeout;
	static CServer*		s_pServer;

#ifndef _MSC_VER
	static void Alarm( int );
#endif // _MSC_VER

	void Listen( );

	int				m_iPort;
	int				m_iTimeout;
	SOCKET			m_iSocket;
	bool			m_fStop;
	IServerClient*	m_pClient;
};

class CServerClientImpl {
public:
	static void* StartRoutine( void* );

	CServerClientImpl( SOCKET, IServerClient* );

private:
	static const size_t	c_iBuffer	= 131072;

	~CServerClientImpl( );

	void StartRoutine( );
	void ReadMessage( );

	SOCKET						m_iSocket;
	IServerClient*				m_pClient;
	std::vector<unsigned char>	m_vecbMessage;
};

}

#endif // SERVERI_H
