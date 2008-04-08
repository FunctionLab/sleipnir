/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "server.h"
#include "serverclient.h"

namespace Sleipnir {

CServer*	CServerImpl::s_pServer		= NULL;
const char*	CServerImpl::c_szPort		= "port";
const char*	CServerImpl::c_szTimeout	= "timeout";

/*!
 * \brief
 * Prepare a server object to listen for incoming connections on the specified port.
 * 
 * \param iPort
 * TCP/IP port on which the server will listen.
 * 
 * \param iTimeout
 * Timeout interval for socket listening.
 * 
 * \param pServerClient
 * Pointer to a template server client object from which new clients will be created to handle incoming
 * requests.
 * 
 * \returns
 * True if the server was initialized successfully.
 * 
 * Prepares the server to listen for connection requests on the given port.  No actual socket manipulation
 * happens until the server is started.
 * 
 * \remarks
 * On platforms where listening on a socket must be interrupted periodically, the interrupt will occur
 * with the given timeout frequency.  This is usually a nonissue, and a value of ~100ms is adequate.
 * 
 * \see
 * Start
 */
bool CServer::Initialize( size_t iPort, size_t iTimeout, IServerClient* pServerClient ) {
#ifndef _MSC_VER
	struct sigaction	Sigact;
#endif // _MSC_VER

	m_pClient = pServerClient;
	m_iPort = iPort;
	m_iTimeout = iTimeout;
#ifndef _MSC_VER
	Sigact.sa_handler = Alarm;
	memset( &Sigact.sa_mask, 0, sizeof(Sigact.sa_mask) );
	Sigact.sa_flags = 0;
	sigaction( SIGALRM, &Sigact, NULL );
#endif // _MSC_VER

	return true; }

#ifndef _MSC_VER
void CServerImpl::Alarm( int iSig ) { }
#endif // _MSC_VER

/*!
 * \brief
 * Opens a server socket and blocks, listening for incoming connections to which server client threads are
 * attached.
 * 
 * \returns
 * True if the server was bound and closed cleanly, false otherwise.
 * 
 * Begins listening on a TCP/IP server socket for incoming connections.  When such a connection is
 * detected, a new IServerClient object as provided to Initialize is created and given control of a new thread
 * in which to handle the request.  The server object then continues to listen for additional connections on
 * the main thread.  This method will block and not return until the server is closed, e.g. by Stop.
 * 
 * \remarks
 * Initialize must be called before the server is started.
 */
bool CServer::Start( ) {
	sockaddr_in	Addr;
	char		cOn;

#ifdef _MSC_VER
	WSADATA		sWSA;
	WSAStartup( MAKEWORD(2, 0), &sWSA );
#endif // _MSC_VER

	m_fStop = false;
	m_iSocket = socket( PF_INET, SOCK_STREAM, 0 );
	cOn = true;
	setsockopt( m_iSocket, SOL_SOCKET, SO_REUSEADDR, &cOn, sizeof(cOn) );

	Addr.sin_family = AF_INET;
	Addr.sin_port = htons( m_iPort );
	Addr.sin_addr.s_addr = INADDR_ANY;
	s_pServer = this;
	if( bind( m_iSocket, (const sockaddr*)&Addr, sizeof(Addr) ) ||
		listen( m_iSocket, INT_MAX ) ) {
#ifdef _MSC_VER
		{
			char	acError[ 1024 ];

			strerror_s( acError, ARRAYSIZE(acError) - 1, errno );
			g_CatSleipnir.error( "CServer::Start( ) bind failed: %s", acError );
		}
#else // _MSC_VER
		g_CatSleipnir.error( "CServer::Start( ) bind failed: %s", strerror( errno ) );
#endif // _MSC_VER
		return false; }

	g_CatSleipnir.notice( "CServer::Start( ) bound to port %d", m_iPort );
	Listen( );
	g_CatSleipnir.info( "CServer::Start( ) preparing to shutdown..." );

#ifdef _MSC_VER
	WSACleanup( );
#endif // _MSC_VER

	return true; }

void CServerImpl::Listen( ) {
	SOCKET				iClient;
	socklen_t			iSize;
	CServerClientImpl*	pClient;
	pthread_t			thrdClient;
	sockaddr_in			Addr;
#ifndef _MSC_VER
	itimerval			Time;

	memset( &Time, 0, sizeof(Time) );
#endif // _MSC_VER
	while( !m_fStop ) {
#ifndef _MSC_VER
		Time.it_value.tv_usec = 1000 * m_iTimeout;
		setitimer( ITIMER_REAL, &Time, NULL );
		Time.it_value.tv_usec = 0;
#endif // _MSC_VER
		iSize = sizeof(Addr);
		if( ( iClient = accept( m_iSocket, (sockaddr*)&Addr, &iSize ) ) == -1 )
			continue;
#ifndef _MSC_VER
		setitimer( ITIMER_REAL, &Time, NULL );
#endif // _MSC_VER

		pClient = new CServerClientImpl( iClient, m_pClient->NewInstance( iClient,
			iSize = ntohl( Addr.sin_addr.s_addr ), ntohs( Addr.sin_port ) ) );
		g_CatSleipnir.info( "CServer::Listen( ) client 0x%08x connected from %d.%d.%d.%d",
			pClient, ( iSize >> 24 ) & 0xFF, ( iSize >> 16 ) & 0xFF, ( iSize >> 8 ) & 0xFF,
			iSize & 0xFF );
		pthread_create( &thrdClient, NULL, CServerClientImpl::StartRoutine, pClient );
		pthread_detach( thrdClient ); } }

}
