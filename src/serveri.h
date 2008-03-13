#ifndef SERVERI_H
#define SERVERI_H

#include "propertyfile.h"

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

	CPropertyFile	m_Config;
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
#ifdef XML_ENABLED
	static const XMLCh	c_pxszLS[];

	xercesc::DOMDocument*		m_pxdocMessage;
	xercesc::DOMBuilder*		m_pxbld;
#else // XML_ENABLED
	static const size_t	c_iBuffer	= 131072;

	std::vector<unsigned char>	m_vecbMessage;
#endif // XML_ENABLED

	~CServerClientImpl( );

	void StartRoutine( );
	void ReadMessage( );

	SOCKET			m_iSocket;
	IServerClient*	m_pClient;
};

#ifdef XML_ENABLED
class CStringInputSource : public xercesc::DOMInputSource {
public:
	CStringInputSource( const std::string& );
	
	const XMLCh* getEncoding( ) const;
	const XMLCh* getPublicId( ) const;
	const XMLCh* getSystemId( ) const;
	const XMLCh* getBaseURI( ) const;
	void setEncoding( const XMLCh* const );
	void setPublicId( const XMLCh* const );
	void setSystemId( const XMLCh* const );
	void setBaseURI( const XMLCh* const );
	void setIssueFatalErrorIfNotFound( const bool );
	bool getIssueFatalErrorIfNotFound( ) const;
	xercesc::BinInputStream* makeStream( ) const;
	void release( );

private:
	const std::string&	m_str;
};

class CBinStringInputStream : public xercesc::BinInputStream {
public:
	CBinStringInputStream( const std::string& );
	
	unsigned int curPos( ) const;
	unsigned int readBytes( XMLByte* const, const unsigned int );

private:
	const std::string&	m_str;
	int					m_iPos;
	bool				m_fHalf;
};
#endif // XML_ENABLED

}

#endif // SERVERI_H
