#include "stdafx.h"
#include "server.h"
#include "serverclient.h"

namespace libBioUtils {

const char*	CServerImpl::c_szPort		= "port";
const char*	CServerImpl::c_szTimeout	= "timeout";
CServer*	CServerImpl::s_pServer		= NULL;

bool CServer::Initialize( const char* szConfig, IServerClient* pClient ) {
	SVariant			Var;
#ifndef _MSC_VER
	struct sigaction	Sigact;
#endif // _MSC_VER

	m_pClient = pClient;
	if( !m_Config.Open( szConfig ) )
		return false;
	
	Var = m_Config.Get( c_szPort );
	if( Var.m_eType != SVariant::eInt )
		return false;
	m_iPort = Var.m_i;

	Var = m_Config.Get( c_szTimeout );
	if( Var.m_eType != SVariant::eInt )
		return false;
	m_iTimeout = Var.m_i;

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
			g_CatBioUtils.error( "CServer::Start( ) bind failed: %s", acError );
		}
#else // _MSC_VER
		g_CatBioUtils.error( "CServer::Start( ) bind failed: %s", strerror( errno ) );
#endif // _MSC_VER
		return false; }

	g_CatBioUtils.notice( "CServer::Start( ) bound to port %d", m_iPort );
	Listen( );
	g_CatBioUtils.info( "CServer::Start( ) preparing to shutdown..." );

#ifdef _MSC_VER
	WSACleanup( );
#endif // _MSC_VER

	return true; }

void CServer::Listen( ) {
	SOCKET			iClient;
	socklen_t		iSize;
	IServerClient*	pClient;
	pthread_t		thrdClient;
	sockaddr_in		Addr;
#ifndef _MSC_VER
	itimerval		Time;

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

		pClient = m_pClient->NewInstance( );
		pClient->SetParent( this );
		pClient->SetSocket( iClient );
		pClient->SetConfig( &m_Config );

		iSize = ntohl( Addr.sin_addr.s_addr );
		g_CatBioUtils.info( "CServer::Listen( ) client 0x%08x connected from %d.%d.%d.%d",
			pClient, ( iSize >> 24 ) & 0xFF, ( iSize >> 16 ) & 0xFF, ( iSize >> 8 ) & 0xFF,
			iSize & 0xFF );
		pthread_create( &thrdClient, NULL, pClient->GetStartRoutine( ), pClient );
		pthread_detach( thrdClient ); } }

void CServer::Stop( ) {

	m_fStop = true; }

}
