#ifndef SERVERCLIENTI_H
#define SERVERCLIENTI_H

namespace Sleipnir {

class CServer;

class CServerClientImpl {
protected:
	static void* StartRoutine( void* );

	CServerClientImpl( );
	~CServerClientImpl( );

	virtual void ProcessMessage( ) = 0;

	SOCKET		m_iSocket;
	CServer*	m_pParent;
	bool		m_fOpen;

private:
	void ReadMessage( );
};

}

#endif // SERVERCLIENTI_H
