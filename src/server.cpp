#include "stdafx.h"
#include "server.h"
#include "serverclient.h"

namespace Sleipnir {

CServer*	CServerImpl::s_pServer		= NULL;
const char*	CServerImpl::c_szPort		= "port";
const char*	CServerImpl::c_szTimeout	= "timeout";

bool CServer::Initialize( const char* szConfig, IServerClient* pClient ) {
	SVariant			VarPort, VarTimeout;

	if( !m_Config.Open( szConfig ) )
		return false;
	VarPort = m_Config.Get( c_szPort );
	if( VarPort.m_eType != SVariant::eInt )
		return false;
	VarTimeout = m_Config.Get( c_szTimeout );
	if( VarTimeout.m_eType != SVariant::eInt )
		return false;

	return Initialize( VarPort.m_i, VarTimeout.m_i, pClient ); }

bool CServer::Initialize( size_t iPort, size_t iTimeout, IServerClient* pClient ) {
#ifndef _MSC_VER
	struct sigaction	Sigact;
#endif // _MSC_VER

	m_pClient = pClient;
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

void CServer::Listen( ) {
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
			iSize = ntohl( Addr.sin_addr.s_addr ), ntohs( Addr.sin_port ), &m_Config ) );
		g_CatSleipnir.info( "CServer::Listen( ) client 0x%08x connected from %d.%d.%d.%d",
			pClient, ( iSize >> 24 ) & 0xFF, ( iSize >> 16 ) & 0xFF, ( iSize >> 8 ) & 0xFF,
			iSize & 0xFF );
		pthread_create( &thrdClient, NULL, CServerClientImpl::StartRoutine, pClient );
		pthread_detach( thrdClient ); } }

void CServer::Stop( ) {

	m_fStop = true; }

}
